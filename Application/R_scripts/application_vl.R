# load WRF functions
source("R/myRF.R")
app_path <- "Application"

synth_dat <- read.csv(file.path(app_path, "data", "synth_dat1.csv"))
head(synth_dat)
syn_dat <- synth_dat[, names(synth_dat) != "out_cd4"]

N <- nrow(syn_dat)
n_fold <- 10

brier <- function(thres, estprob, obsval) {
  mean((estprob-as.integer(obsval<thres))^2,na.rm=TRUE) }

# 10-fold
set.seed(42)
cvfold_idx <- sample(1:n_fold, N, replace=TRUE)
table(cvfold_idx)
#    1    2    3    4    5    6    7    8    9   10 
# 1140 1089 1079 1112 1124 1055 1103 1057 1081 1230 


##### Predict 6-month viral load (VL) #####
nodesize <- 20
vl_threshold <- 80
vl_cv_wrf <- replicate(n_fold, list(), simplify = FALSE)
vl_cv_wrf$brier <- numeric(n_fold)
for (i in seq_len(n_fold)) {
  dat_train <- syn_dat[cvfold_idx!=i,]
  dat_test <- syn_dat[cvfold_idx==i,]
  
  set.seed(42)
  fit_wrf <- myRF(out_VL~., data=dat_train, ntree=1000, nodesize=nodesize, split_rule="Wilcoxon")
  vl_cv_wrf[[i]]$predvl <- predict(fit_wrf, newdata=dat_test, what="probability", outcomes=vl_threshold,
                                   linear_intpl=TRUE, include=FALSE, lower.tail=TRUE)
  vl_cv_wrf$brier[i] <- brier(vl_threshold, vl_cv_wrf[[i]]$predvl[,1], dat_test$out_VL)
}

# saveRDS(vl_cv_wrf, file=file.path(app_path, "results", "vl_cv_wrf.rds"))
print(mean(vl_cv_wrf$brier)) # 0.1271001




##### Compare with logistic regression & classification forest #####
### Logistic regression with restricted cubic splines for continuous covariates
library(rms)
packageVersion("rms") # 8.1-0
vl_lr_brier <- numeric(10)
for (i in seq_len(n_fold)) {
  fit_lr <- lrm(I(out_VL < 80) ~ rcs(age_art, 3) + male_y + site + infmode +
                  prior_aids + rcs(sqrt(base_cd4), 5) + rcs(log10(base_VL), 5) +
                  art_regimen + rcs(months_to_vl, 3) + rcs(year_art, 3),
                data = syn_dat[cvfold_idx!=i,])
  vl_lr_predprob <- predict(fit_lr, newdata=syn_dat[cvfold_idx==i,], type="fitted")
  vl_lr_brier[i] <- brier(80, vl_lr_predprob, syn_dat[cvfold_idx==i,]$out_VL)
}
print(mean(vl_lr_brier)) # 0.1309335



### Classification forest
library(ranger)
packageVersion("ranger") # 0.18.0
syn_dat$vl_supp <- factor(ifelse(syn_dat$out_VL < 80, "yes", "no"),
                         levels = c("no", "yes"))
vl_rf_brier <- numeric(10)
vl_rf_predprob <- numeric(nrow(syn_dat))

for (i in seq_len(n_fold)) {
  dat_train <- syn_dat[cvfold_idx != i, ]
  dat_test  <- syn_dat[cvfold_idx == i, ]
  
  fit_rf <- ranger(
    vl_supp ~ age_art + male_y + site + infmode +
      prior_aids + base_cd4 + base_VL +
      art_regimen + months_to_vl + year_art,
    data = dat_train, probability = TRUE, num.trees = 1000,
    mtry = 3, min.node.size = 20, sample.fraction = 1-exp(-1),
    replace = FALSE, seed = 42 + i)
  
  rf_predprob <- predict(fit_rf, data = dat_test)$predictions[, "yes"]
  
  vl_rf_predprob[cvfold_idx == i] <- rf_predprob
  vl_rf_brier[i] <- brier(80, rf_predprob, dat_test$out_VL)
}
print(mean(vl_rf_brier)) # 0.1228797

p_vl_below80 <- mean(syn_dat$out_VL<80)
vl_brier_comparison <- data.frame(
  Method = c("Uninformative model",
             "Logistic regression with restricted cubic splines",
             "WRF",
             "Classification random forest"),
  Mean_Brier = c(p_vl_below80 * (1-p_vl_below80),
                 mean(vl_lr_brier),
                 mean(vl_cv_wrf$brier),
                 mean(vl_rf_brier)))
write.csv(vl_brier_comparison, file=file.path(app_path, "results", "vl_brier_comparison.csv"), row.names=FALSE)



##### Figure: ECDF for VL #####
two_profiles <- read.csv(file.path(app_path, "results", "two_profiles.csv"))
syn_dat <- syn_dat[, names(syn_dat) != "vl_supp"]
set.seed(42)
vl_wrf_fit <- myRF(out_VL~., data=syn_dat, ntree=1000, mtry=3, nodesize=nodesize, split_rule="Wilcoxon")
vl_ecdf_wrf <- list()
vl_ecdf_wrf$subj1 <- predict(vl_wrf_fit, newdata=two_profiles[1,], what="ECDF")
vl_ecdf_wrf$subj2 <- predict(vl_wrf_fit, newdata=two_profiles[2,], what="ECDF")
vl_ecdf_wrf$subj1$y.sorted <- log(vl_ecdf_wrf$subj1$y.sorted)
vl_ecdf_wrf$subj2$y.sorted <- log(vl_ecdf_wrf$subj2$y.sorted)

dat_vl <- cn_dat[,names(cn_dat) %in% c("art_regimen","site","male_y","age_art","year_art","base_VL",
                                       "base_cd4","months_to_vl","out_VL","prior_aids","infmode")]

dat_vl <- dat_vl[, c("art_regimen","site","male_y","age_art","year_art","base_VL",
                     "base_cd4","months_to_vl","prior_aids","infmode","out_VL")]
dat_vl$out_VL[dat_vl$out_VL<80] <- 40
dat_vl <- as.data.frame(dat_vl)

set.seed(42)
vl_wrf_fit <- myRF(out_VL~., data=dat_vl, ntree=1000, mtry=3, nodesize=20, split_rule="Wilcoxon")
vl_ecdf_wrf <- list()
vl_ecdf_wrf$subj1 <- predict(vl_wrf_fit, newdata=cd4_2subj[1,1:10], what="ECDF")
vl_ecdf_wrf$subj2 <- predict(vl_wrf_fit, newdata=cd4_2subj[2,1:10], what="ECDF")
vl_ecdf_wrf$subj1$y.sorted <- log(vl_ecdf_wrf$subj1$y.sorted)
vl_ecdf_wrf$subj2$y.sorted <- log(vl_ecdf_wrf$subj2$y.sorted)

library(ggplot2)
make_ecdf_df <- function(ecdf_list) {
  subjects <- c("Baseline CD4=200, VL=1e6",
                "Baseline CD4=500, VL=1e4")
  do.call(rbind, lapply(1:2, function(i) {
    subj <- ecdf_list[[i]]
    data.frame(y = subj$y.sorted,
               ECDF = subj$ECDF,
               Subject = subjects[i]) }))
}

vlcdf_pl_d <- make_ecdf_df(vl_ecdf_wrf)
vlcdf_pl_d$Subject <- factor(vlcdf_pl_d$Subject,
                             levels = c("Baseline CD4=200, VL=1e6",
                                        "Baseline CD4=500, VL=1e4"))

vl_80_x <- log(79)
vl_axis_values_all <- c(79, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7)
vl_axis_breaks_all <- log(vl_axis_values_all)
vl_axis_labels_all <- c("<80", "10^2", "10^3", "10^4", "10^5", "10^6", "10^7")
vl_x_min <- min(vlcdf_pl_d$y, na.rm = TRUE)
vl_x_max <- max(vlcdf_pl_d$y, na.rm = TRUE)
keep_ticks <- vl_axis_breaks_all >= vl_x_min - 1e-8 &
  vl_axis_breaks_all <= vl_x_max + 1e-8

vl_axis_breaks <- vl_axis_breaks_all[keep_ticks]
vl_axis_labels <- vl_axis_labels_all[keep_ticks]
line_types <- c(
  "Baseline CD4=200, VL=1e6" = "solid",
  "Baseline CD4=500, VL=1e4" = "dashed")

fig1_vl <- ggplot(vlcdf_pl_d, aes(x = y, y = ECDF, linetype = Subject)) +
  geom_step(linewidth = 0.7) +
  scale_linetype_manual(
    values = line_types, breaks = names(line_types),
    labels = legend_labels) +
  scale_x_continuous(
    breaks = vl_axis_breaks,
    labels = vl_axis_labels,
    expand = expansion(mult = c(0.02, 0.03))) +
  scale_y_continuous(
    limits = c(0, 1), breaks = seq(0, 1, by = 0.25),
    labels = sprintf("%.2f", seq(0, 1, by = 0.25))) +
  labs(
    x = "6-month HIV-1 RNA viral load (copies/mL; log scale above 80)",
    y = NULL, linetype = NULL) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    legend.key.width = unit(2, "cm"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()) + guides(linetype = guide_legend(ncol = 2))

ggsave(filename = file.path(app_path, "figures", "Figure_vl_ecdf.pdf"),
       plot = fig1_vl, width = 10, height = 6)



##### Variable importance based on OOB Brier score #####
n_rep <- 10
vl_threshold <- 80
covar_names <- vl_wrf_fit$xnames
n_covar <- length(covar_names)
y <- vl_wrf_fit$y
X <- vl_wrf_fit$X
ntree <- vl_wrf_fit$ntree
vl_80 <- as.integer(y < vl_threshold)

# Store the leaf membership of all training observations under each fitted tree,
# identify which variables are used by each tree.
all_leaf_ids <- vector("list", ntree)
used_variables <- vector("list", ntree)
for (t in seq_len(ntree)) {
  tree_t <- vl_wrf_fit$RF[[t]]
  all_leaf_ids[[t]] <- if (!is.null(tree_t$leaf_id_by_obs)) {
    tree_t$leaf_id_by_obs
  } else {
    .get_leaf_ids_vectorized(tree_t, X)
  }
  used_variables[[t]] <- .get_tree_used_vars(tree_t)
}


# Obtain baseline tree-level OOB Brier scores
baseline_tree_brier <- rep(NA_real_, ntree)

for (t in seq_len(ntree)) {
  tree_weights <- .tree_oob_weights_for_vimp(
    object = vl_wrf_fit, tree_index = t,
    response_members = "all",
    lids_base = all_leaf_ids[[t]],
    used_vars = used_variables[[t]])
  oob_index <- tree_weights$oob_idx
  
  if (length(oob_index) == 0L) {
    next
  }
  
  # weighted P(VL<80|X) for each OOB observation.
  predicted_probability <- as.numeric(tree_weights$W %*% vl_80)
  baseline_tree_brier[t] <- mean((predicted_probability - vl_80[oob_index])^2)
}

baseline_mean_brier <- mean(baseline_tree_brier) # 0.1349043


## Permutation OOB Brier scores - variable importance
permuted_mean_brier <- matrix(NA_real_, nrow = n_rep, ncol = n_covar,
                              dimnames = list(paste0("rep", seq_len(n_rep)), covar_names))
set.seed(42)
for (r in seq_len(n_rep)) {
  cat("Permutation repetition", r, "of", n_rep, "\n")
  
  for (j in seq_len(n_covar)) {
    permuted_tree_brier <- baseline_tree_brier
    for (t in which(valid_trees)) {
      if (!(j %in% used_variables[[t]])) {
        next
      }
      tree_weights_permuted <- .tree_oob_weights_for_vimp(
        object = vl_wrf_fit, tree_index = t, response_members = "all",
        lids_base = all_leaf_ids[[t]], used_vars = used_variables[[t]],
        permute_var = j)
      oob_index <- tree_weights_permuted$oob_idx
      predicted_probability_permuted <- as.numeric(
        tree_weights_permuted$W %*% vl_80)
      permuted_tree_brier[t] <- mean((predicted_probability_permuted-vl_80[oob_index])^2)
    }
    permuted_mean_brier[r, j] <- mean(permuted_tree_brier[valid_trees])
    cat("  ",covar_names[j],": mean OOB Brier =",signif(permuted_mean_brier[r, j], 6),"\n")
  }
}


## summarize variable importance
mean_permuted_brier <- colMeans(permuted_mean_brier, na.rm = TRUE)

vl_variable_importance <- data.frame(
  Variable = covar_names,
  Baseline_mean_Brier = baseline_mean_brier,
  Permuted_mean_Brier = as.numeric(mean_permuted_brier),
  VI_Brier = as.numeric(mean_permuted_brier - baseline_mean_brier),
  VI_percent = as.numeric(100*(mean_permuted_brier-baseline_mean_brier)/baseline_mean_brier),
  row.names = NULL)

# Sort from most to least important.
vl_variable_importance <- vl_variable_importance[
  order(vl_variable_importance$VI_percent, decreasing = TRUE),]
rownames(vl_variable_importance) <- NULL
print(vl_variable_importance)
#        Variable Baseline_mean_Brier Permuted_mean_Brier     VI_Brier VI_percent
# 1       base_VL           0.1349043           0.1476429 0.0127386171  9.4427085
# 2      base_cd4           0.1349043           0.1411816 0.0062773018  4.6531527
# 3  months_to_vl           0.1349043           0.1397068 0.0048024925  3.5599262
# 4          site           0.1349043           0.1371961 0.0022918737  1.6988889
# 5    prior_aids           0.1349043           0.1363634 0.0014591044  1.0815850
# 6   art_regimen           0.1349043           0.1362825 0.0013781821  1.0216001
# 7       age_art           0.1349043           0.1361764 0.0012721232  0.9429822
# 8       infmode           0.1349043           0.1361732 0.0012689539  0.9406329
# 9      year_art           0.1349043           0.1355839 0.0006796223  0.5037811
# 10       male_y           0.1349043           0.1354731 0.0005688666  0.4216817


write.csv(vl_variable_importance, 
          file = file.path(app_path, "results", "vl_variable_importance.csv"), row.names = FALSE)
vl_wrf_vimp <- list(
  SummaryTab = vl_variable_importance,
  permuted_mean_Brier = permuted_mean_brier,
  baseline_tree_Brier = baseline_tree_brier,
  baseline_mean_Brier = baseline_mean_brier,
  response_members = "all", threshold = vl_threshold,
  repetitions = n_rep, seed = 42)
# saveRDS(vl_wrf_vimp, file=file.path(app_path,"results","vl_wrf_vimp.rds"))
