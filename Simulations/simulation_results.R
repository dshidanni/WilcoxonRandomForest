rm(list = ls())

sim_dir <- "Simulations"
dat_test_file  <- file.path(sim_dir, "dat_test.RData")
dat_test2_file <- file.path(sim_dir, "dat_test2.RData")

s1_tree_dir   <- file.path(sim_dir, "s1_trees")
s1_forest_dir <- file.path(sim_dir, "s1_forests")
s2_forest_dir <- file.path(sim_dir, "s2")
cpm_forest_dir <- file.path(sim_dir, "s1_CPMforests")

RMSE <- function(est, trueval) {
  sqrt(mean((est - trueval)^2, na.rm = TRUE))
}

meanbias <- function(est, trueval) {
  mean(est - trueval, na.rm = TRUE)
}

pinball <- function(est, obsval, q) {
  mean((obsval - est) * ifelse(obsval >= est, q, q - 1), na.rm = TRUE)
}

brier <- function(thres, estprob, obsval) {
  mean((estprob - as.integer(obsval <= thres))^2, na.rm = TRUE)
}

load_one_object <- function(file, object_name) {
  if (!file.exists(file)) stop("File not found: ", file)
  envn <- new.env(parent = emptyenv())
  loaded <- load(file, envir = envn)
  if (!(object_name %in% loaded)) {
    stop("Object '", object_name, "' was not found in ", file)
  }
  envn[[object_name]]
}

load_rdata_dir <- function(path) {
  if (!dir.exists(path)) stop("Directory not found: ", path)

  files <- list.files(path, pattern = "\\.RData$", full.names = TRUE)
  if (length(files) == 0L) stop("No .RData files found in: ", path)

  envn <- new.env(parent = emptyenv())
  for (f in files) load(f, envir = envn)
  envn
}

get_obj <- function(envn, name) {
  if (!exists(name, envir = envn, inherits = FALSE)) {
    stop("Required object '", name, "' is missing from the loaded directory.")
  }
  get(name, envir = envn, inherits = FALSE)
}

column_truth_metric <- function(est, trueval, metric_fun) {
  est <- as.matrix(est)
  mean(vapply(
    seq_len(ncol(est)),
    function(i) metric_fun(est[, i], trueval),
    numeric(1)))
}

paired_replicate_metric <- function(est, obs, metric_fun, ...) {
  est <- as.matrix(est)
  obs <- as.matrix(obs)

  if (!identical(dim(est), dim(obs))) {
    stop("Prediction and observed outcome matrices must have the same dimensions.")
  }

  mean(vapply(
    seq_len(ncol(est)),
    function(i) metric_fun(est[, i], obs[, i], ...),
    numeric(1)))
}

subset_prediction_matrix <- function(x, sim_idx) {
  x <- as.matrix(x)
  if (max(sim_idx) > ncol(x)) stop("Requested simulation column is unavailable.")
  x[, sim_idx, drop = FALSE]
}

subset_scalar_matrix <- function(x, sim_idx) {
  x <- as.matrix(x)
  if (max(sim_idx) > nrow(x)) stop("Requested simulation row is unavailable.")
  x[sim_idx, , drop = FALSE]
}


# Evaluation metrics:
# (mean) RMSE, bias, mean & median of CRPS, coverage/non-exceedance rate
summarize_method <- function(results_env, prefix,
                             truth, test_y, sim_idx,
                             threshold = NULL, base_time = NULL) {
  q10 <- subset_prediction_matrix(
    get_obj(results_env, paste0(prefix, "_qt_10")), sim_idx)
  q50 <- subset_prediction_matrix(
    get_obj(results_env, paste0(prefix, "_qt_50")), sim_idx)
  q90 <- subset_prediction_matrix(
    get_obj(results_env, paste0(prefix, "_qt_90")), sim_idx)
  
  mean_est <- subset_prediction_matrix(
    get_obj(results_env, paste0(prefix, "_mean")), sim_idx)
  crps <- subset_prediction_matrix(
    get_obj(results_env, paste0(prefix, "_crps")), sim_idx)
  cvm <- subset_scalar_matrix(
    get_obj(results_env, paste0(prefix, "_cvm")), sim_idx)
  elapsed <- subset_scalar_matrix(
    get_obj(results_env, paste0(prefix, "_time")), sim_idx)[, 1]

  obs <- subset_prediction_matrix(test_y, sim_idx)

  out <- c(
    qt10_RMSE = column_truth_metric(q10, truth$q10, RMSE),
    qt10_bias = column_truth_metric(q10, truth$q10, meanbias),
    qt10_pinball = paired_replicate_metric(q10, obs, pinball, q = 0.1),

    qt50_RMSE = column_truth_metric(q50, truth$q50, RMSE),
    qt50_bias = column_truth_metric(q50, truth$q50, meanbias),
    qt50_pinball = paired_replicate_metric(q50, obs, pinball, q = 0.5),

    qt90_RMSE = column_truth_metric(q90, truth$q90, RMSE),
    qt90_bias = column_truth_metric(q90, truth$q90, meanbias),
    qt90_pinball = paired_replicate_metric(q90, obs, pinball, q = 0.9),

    mean_RMSE = column_truth_metric(mean_est, truth$mean, RMSE),
    mean_bias = column_truth_metric(mean_est, truth$mean, meanbias))

  if (!is.null(threshold)) {
    prob <- subset_prediction_matrix(
      get_obj(results_env, paste0(prefix, "_prob")), sim_idx)
    if (is.null(truth$prob) || length(truth$prob) != nrow(prob)) {
      stop("Probability truth vector is missing or has an incompatible length.")
    }

    out <- c(out,
      prob_RMSE = column_truth_metric(prob, truth$prob, RMSE),
      prob_bias = column_truth_metric(prob, truth$prob, meanbias),
      prob_brier = paired_replicate_metric(prob, obs,
        function(est, y) brier(threshold, est, y)))
  }

  # Match original aggregation exactly.
  out <- c(out,
    CRPS_mean = mean(colMeans(crps)),
    CRPS_med = mean(apply(crps, 2, median)),
    CvM_mean = mean(colMeans(cvm)),
    covg_rate = mean(colMeans(q10 <= obs & q90 >= obs)),
    covgwid_med = mean(apply(q90 - q10, 2, median)),
    nonexc10 = mean(colMeans(obs <= q10)),
    nonexc90 = mean(colMeans(obs <= q90)))

  # For CPM-refined tree, it's Wilcoxon tree + CPM refinement time
  if (!is.null(base_time)) {
    base_elapsed <- subset_scalar_matrix(base_time, sim_idx)[, 1]
    elapsed <- elapsed + base_elapsed
  }
  c(out, time = mean(elapsed))
}

# build the summary data table
build_summary <- function(result_dir,
                          ntrains, nodesizes, model_map, truth,
                          transformed = FALSE, nsim, threshold = NULL) {
  rows <- vector("list", length(ntrains) * length(model_map))
  k <- 0

  for (i in seq_along(ntrains)) {
    nt <- ntrains[i]
    envn <- load_rdata_dir(file.path(result_dir, as.character(nt)))

    obs_name <- if (transformed) "test_Ytr" else "test_Y"
    test_y <- get_obj(envn, obs_name)
    sim_idx <- seq_len(nsim)

    for (model_name in names(model_map)) {
      k <- k + 1L
      prefix <- unname(model_map[[model_name]])
      if (transformed) prefix <- paste0("tr", prefix)

      metrics <- summarize_method(
        results_env = envn, prefix = prefix,
        truth = truth, test_y = test_y,
        sim_idx = sim_idx, threshold = threshold)

      rows[[k]] <- data.frame(
        Model = model_name, ntrain = nt, node = nodesizes[i],
        as.list(metrics), check.names = FALSE, row.names = NULL)
    }

    rm(envn)
    gc(verbose = FALSE)
  }

  do.call(rbind, rows)
}




##### S1.1: load true values #####
dat_test <- load_one_object(dat_test_file, "dat_test")

truth_s11 <- list(
  q10 = dat_test$true10, q50 = dat_test$true50, q90 = dat_test$true90,
  mean = dat_test$truemean, prob = dat_test$trueprob)

truth_s12 <- list(
  q10 = dat_test$trtrue10, q50 = dat_test$trtrue50, q90 = dat_test$trtrue90,
  mean = dat_test$trtruemean,
  # c_Ytr = h(0), so this is the same latent threshold probability.
  prob = dat_test$trueprob)

c_Y <- 0
c_Ytr <- qchisq(pnorm(0), df = 15)

ntrains <- c(1000, 2000, 5000)



##### S1.2: Full single-tree results #####
tree_models <- c(
  "SSE" = "sse", "Wilcoxon" = "wil", "Pinball" = "pbl",
  "CPM-probit" = "cpmp", "CPM-logit" = "cpml")

s11_tree_res <- build_summary(
  result_dir = s1_tree_dir, ntrains = ntrains, nodesizes = c(60, 120, 300),
  model_map = tree_models, truth = truth_s11,
  transformed = FALSE, nsim = 1000, threshold = c_Y)
s11_tree_res3 <- s11_tree_res; s11_tree_res3[, 4:25] <- round(s11_tree_res3[, 4:25], 3)
write.csv(s11_tree_res, 
          file = file.path(sim_dir, "s11_tree_res.csv"), row.names = FALSE)
write.csv(s11_tree_res3, 
          file = file.path(sim_dir, "s11_tree_res3.csv"), row.names = FALSE)


s12_tree_res <- build_summary(
  result_dir = s1_tree_dir, ntrains = ntrains, nodesizes = c(60, 120, 300),
  model_map = tree_models, truth = truth_s12,
  transformed = TRUE, nsim = 1000, threshold = c_Ytr)
s12_tree_res3 <- s12_tree_res; s12_tree_res3[, 4:25] <- round(s12_tree_res3[, 4:25], 3)
write.csv(s12_tree_res, 
          file = file.path(sim_dir, "s12_tree_res.csv"), row.names = FALSE)
write.csv(s12_tree_res3, 
          file = file.path(sim_dir, "s12_tree_res3.csv"), row.names = FALSE)




##### S2.3: Full WRF/QRF forest results for Scenarios 1.1 and 1.2 #####
forest_models <- c("QRF" = "sse", "WRF" = "wil")

s11_forest_res <- build_summary(
  result_dir = s1_forest_dir, ntrains = ntrains, nodesizes = c(10, 20, 50),
  model_map = forest_models, truth = truth_s11,
  transformed = FALSE, nsim = 250, threshold = c_Y)
s11_forest_res3 <- s11_forest_res; s11_forest_res3[, 4:25] <- round(s11_forest_res3[, 4:25], 3)


s12_forest_res <- build_summary(
  result_dir = s1_forest_dir, ntrains = ntrains, nodesizes = c(10, 20, 50),
  model_map = forest_models, truth = truth_s12,
  transformed = TRUE, nsim = 250, threshold = c_Ytr)
s12_forest_res3 <- s12_forest_res; s12_forest_res3[, 4:25] <- round(s12_forest_res3[, 4:25], 3)
write.csv(s12_forest_res, 
          file = file.path(sim_dir, "s12_forest_res.csv"), row.names = FALSE)
write.csv(s12_forest_res3, 
          file = file.path(sim_dir, "s12_forest_res3.csv"), row.names = FALSE)




##### S2.3: Scenario 2 forest results
dat_test2 <- load_one_object(dat_test2_file, "dat_test2")

truth_s2 <- list(q10 = dat_test2$true10, q50 = dat_test2$true50,
                 q90 = dat_test2$true90, mean = dat_test2$truemean)

s2_forest_res <- build_summary(
  result_dir = s2_forest_dir, ntrains = ntrains, nodesizes = c(5, 10, 25),
  model_map = forest_models, truth = truth_s2, transformed = FALSE,
  nsim = 250, threshold = NULL)
s2_forest_res3 <- s2_forest_res; s2_forest_res3[, 4:22] <- round(s2_forest_res3[, 4:22], 3)
write.csv(s2_forest_res, 
          file = file.path(sim_dir, "s2_forest_res.csv"), row.names = FALSE)
write.csv(s2_forest_res3, 
          file = file.path(sim_dir, "s2_forest_res3.csv"), row.names = FALSE)



##### S2.4: Exploratory CPM-refined forest analysis #####
# Compare the CPM-refined WRF against the corresponding ECDF-based WRF using
# the same first 100 Scenario 1 forest replicates (seeds 43,...,142).
base_e <- load_rdata_dir(file.path(s1_forest_dir, "2000"))
cpm_e <- load_rdata_dir(cpm_forest_dir)
sim_idx_cpm <- seq_len(100)

test_Y_2000 <- get_obj(base_e, "test_Y")
test_Ytr_2000 <- get_obj(base_e, "test_Ytr")

# Scenario 1.1
wrf_s11_metrics <- summarize_method(
  results_env = base_e, prefix = "wil", truth = truth_s11,
  test_y = test_Y_2000, sim_idx = sim_idx_cpm, threshold = c_Y)

cpm_s11_metrics <- summarize_method(
  results_env = cpm_e, prefix = "cpmp", truth = truth_s11,
  test_y = test_Y_2000, sim_idx = sim_idx_cpm, threshold = c_Y,
  base_time = get_obj(base_e, "wil_time"))

s11_cpmforest_res <- data.frame(
  Model = c("WRF", "CPM-refined WRF"), ntrain = 2000, node = 20,
  rbind(wrf_s11_metrics, cpm_s11_metrics), check.names = FALSE, row.names = NULL)
s11_cpmforest_res3 <- s11_cpmforest_res; s11_cpmforest_res3[, 4:25] <- round(s11_cpmforest_res3[, 4:25], 3)
write.csv(s11_cpmforest_res, 
          file = file.path(sim_dir, "s11_cpmforest_res.csv"), row.names = FALSE)
write.csv(s11_cpmforest_res3, 
          file = file.path(sim_dir, "s11_cpmforest_res3.csv"), row.names = FALSE)


# Scenario 1.2
wrf_s12_metrics <- summarize_method(
  results_env = base_e, prefix = "trwil", truth = truth_s12,
  test_y = test_Ytr_2000, sim_idx = sim_idx_cpm, threshold = c_Ytr)

cpm_s12_metrics <- summarize_method(
  results_env = cpm_e, prefix = "trcpmp", truth = truth_s12,
  test_y = test_Ytr_2000, sim_idx = sim_idx_cpm, threshold = c_Ytr,
  base_time = get_obj(base_e, "trwil_time"))

s12_cpmforest_res <- data.frame(
  Model = c("WRF", "CPM-refined WRF"), ntrain = 2000, node = 20,
  rbind(wrf_s12_metrics, cpm_s12_metrics), check.names = FALSE, row.names = NULL)
s12_cpmforest_res3 <- s12_cpmforest_res; s12_cpmforest_res3[,4:25] <- round(s12_cpmforest_res[,4:25],3)
write.csv(s12_cpmforest_res, 
          file = file.path(sim_dir, "s12_cpmforest_res.csv"), row.names = FALSE)
write.csv(s12_cpmforest_res3, 
          file = file.path(sim_dir, "s12_cpmforest_res3.csv"), row.names = FALSE)


# Combined file for convenient construction of the Supplementary S2.4 table.
cpmforest_res <- rbind(
  data.frame(Scenario = "1.1", s11_cpmforest_res, check.names = FALSE),
  data.frame(Scenario = "1.2", s12_cpmforest_res, check.names = FALSE))
cpmforest_res3 <- cpmforest_res; cpmforest_res3[, 5:26] <- round(cpmforest_res3[,5:26],3)
write.csv(cpmforest_res, 
          file = file.path(sim_dir, "cpmforest_res.csv"), row.names = FALSE)
write.csv(cpmforest_res3, 
          file = file.path(sim_dir, "cpmforest_res3.csv"), row.names = FALSE)

