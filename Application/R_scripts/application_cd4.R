# load WRF functions
source("R/myRF.R")
app_path <- "Application"

synth_dat <- read.csv(file.path(app_path, "data", "synth_dat1.csv"))
head(synth_dat)
syn_dat <- synth_dat[, names(synth_dat) != "out_VL"]

N <- nrow(syn_dat)
n_fold <- 10
cd4vals <- c(200, 350, 500)
nodesizes <- c(5, 10, 20, 50, 100, 200)

brier <- function(thres, estprob, obsval) {
  mean((estprob-as.integer(obsval<=thres))^2,na.rm=TRUE) }

# 10-fold
set.seed(42)
cvfold_idx <- sample(1:n_fold, N, replace=TRUE)
table(cvfold_idx)
#    1    2    3    4    5    6    7    8    9   10 
# 1140 1089 1079 1112 1124 1055 1103 1057 1081 1230 


##### Predict 6-month CD4 threshold probabilities #####
nodesize <- 5
cd4_cv_qrf <- cd4_cv_wrf <- replicate(10, list(), simplify = FALSE)
cd4_cv_qrf$meanCRPS <- cd4_cv_wrf$meanCRPS <- 
  cd4_cv_qrf$medCRPS <- cd4_cv_wrf$medCRPS <- 
  cd4_cv_qrf$cvm <- cd4_cv_wrf$cvm <- numeric(10)

for (i in seq_len(n_fold)) {
  dat_train <- syn_dat[cvfold_idx!=i,]
  dat_test <- syn_dat[cvfold_idx==i,]
  
  set.seed(42)
  fit_qrf <- myRF(out_cd4~., data=dat_train, ntree=1000, nodesize=nodesize, split_rule="SSE")
  cd4_cv_qrf[[i]]$predcd4 <- predict(fit_qrf, newdata=dat_test, what="probability", 
                                     outcomes=cd4vals, linear_intpl=TRUE, 
                                     include=TRUE, lower.tail=TRUE)
  cd4_cv_qrf[[i]]$crps <- predict(fit_qrf, newdata=dat_test, what="CRPS", y_obsvd=dat_test$out_cd4)
  cd4_cv_qrf$cvm[i] <- predict(fit_qrf, newdata=dat_test, what="CvM", y_obsvd=dat_test$out_cd4)
  cd4_cv_qrf[[i]]$brier <- data.frame(CD4=cd4vals, 
                                      brierscore=c(brier(200, cd4_cv_qrf[[i]]$predcd4[,1], dat_test$out_cd4),
                                                   brier(350, cd4_cv_qrf[[i]]$predcd4[,2], dat_test$out_cd4),
                                                   brier(500, cd4_cv_qrf[[i]]$predcd4[,3], dat_test$out_cd4)))
  cd4_cv_qrf$meanCRPS[i] <- mean(cd4_cv_qrf[[i]]$crps)
  cd4_cv_qrf$medCRPS[i] <- median(cd4_cv_qrf[[i]]$crps)
  
  set.seed(42)
  fit_wrf <- myRF(out_cd4~., data=dat_train, ntree=1000, nodesize=nodesize, split_rule="Wilcoxon")
  cd4_cv_wrf[[i]]$predcd4 <- predict(fit_wrf, newdata=dat_test, what="probability", 
                                     outcomes=cd4vals, linear_intpl=TRUE, 
                                     include=TRUE, lower.tail=TRUE)
  cd4_cv_wrf[[i]]$crps <- predict(fit_wrf, newdata=dat_test, what="CRPS", y_obsvd=dat_test$out_cd4)
  cd4_cv_wrf$cvm[i] <- predict(fit_wrf, newdata=dat_test, what="CvM", y_obsvd=dat_test$out_cd4)
  cd4_cv_wrf[[i]]$brier <- data.frame(CD4=cd4vals, 
                                      brierscore=c(brier(200, cd4_cv_wrf[[i]]$predcd4[,1], dat_test$out_cd4),
                                                   brier(350, cd4_cv_wrf[[i]]$predcd4[,2], dat_test$out_cd4),
                                                   brier(500, cd4_cv_wrf[[i]]$predcd4[,3], dat_test$out_cd4)))
  cd4_cv_wrf$meanCRPS[i] <- mean(cd4_cv_wrf[[i]]$crps)
  cd4_cv_wrf$medCRPS[i] <- median(cd4_cv_wrf[[i]]$crps)
}
# saveRDS(cd4_cv_qrf, file=file.path(app_path, "results", "cd4_cv_qrf.rds"))
# saveRDS(cd4_cv_wrf, file=file.path(app_path, "results", "cd4_cv_wrf.rds"))


##### Summarize the 10-fold out-of-fold predictions

dir.create(file.path(app_path, "results"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(app_path, "figures"),
           recursive = TRUE, showWarnings = FALSE)

# Store each subject's out-of-fold predictions in the original row order.
qrf_prob <- matrix(NA_real_, nrow = N, ncol = length(cd4vals),
                   dimnames = list(NULL, paste0("P_CD4_le_", cd4vals)))
wrf_prob <- matrix(NA_real_, nrow = N, ncol = length(cd4vals),
                   dimnames = list(NULL, paste0("P_CD4_le_", cd4vals)))
qrf_crps <- rep(NA_real_, N)
wrf_crps <- rep(NA_real_, N)

for (i in seq_len(n_fold)) {
  test_idx <- which(cvfold_idx == i)
  qrf_pred_i <- as.matrix(cd4_cv_qrf[[i]]$predcd4)
  wrf_pred_i <- as.matrix(cd4_cv_wrf[[i]]$predcd4)
  qrf_prob[test_idx, ] <- qrf_pred_i
  wrf_prob[test_idx, ] <- wrf_pred_i
  qrf_crps[test_idx] <- cd4_cv_qrf[[i]]$crps
  wrf_crps[test_idx] <- cd4_cv_wrf[[i]]$crps
}


# Calculate the Brier scores from all pooled out-of-fold predictions.
# This appropriately accounts for the unequal fold sizes.
qrf_brier <- vapply(seq_along(cd4vals),
                    function(j) { brier(thres = cd4vals[j], estprob = qrf_prob[, j], 
                                        obsval = syn_dat$out_cd4)},
                    numeric(1))
# 0.08055242 0.10011121 0.11071733
wrf_brier <- vapply(seq_along(cd4vals),
                    function(j) { brier(thres = cd4vals[j], estprob = wrf_prob[, j],
                                        obsval = syn_dat$out_cd4)}, 
                    numeric(1))
# 0.07908923 0.09895795 0.11095063

# Your stored CvM values are fold-specific statistics.
# The table below reports their sum across the 10 validation folds.
cd4_cv_results <- data.frame(Method = c("QRF", "WRF"),
                             Mean_CRPS = c(mean(qrf_crps), mean(wrf_crps)),
                             Median_CRPS = c(median(qrf_crps), median(wrf_crps)),
                             Mean_CvM = c(mean(cd4_cv_qrf$cvm), mean(cd4_cv_wrf$cvm)),
                             Brier_CD4_200 = c(qrf_brier[1], wrf_brier[1]),
                             Brier_CD4_350 = c(qrf_brier[2], wrf_brier[2]),
                             Brier_CD4_500 = c(qrf_brier[3], wrf_brier[3]))
#   Method Mean_CRPS Median_CRPS Mean_CvM Brier_CD4_200 Brier_CD4_350 Brier_CD4_500
# 1    QRF  83.89909    57.99335 1.295340    0.08055242    0.10011121     0.1107173
# 2    WRF  84.28489    57.93124 1.210231    0.07908923    0.09895795     0.1109506

write.csv(cd4_cv_results, file = file.path(app_path, "results", "cd4_cv_results.csv"), row.names = FALSE)




##### Variable importance based on OOB CRPS #####
set.seed(42)
qrf_fit <- myRF(out_cd4~., data=syn_dat, ntree=1000, mtry=3, nodesize=5,
                split_rule="SSE", replace=FALSE)
cd4_qrf_vimp <- variable_importance_crps.myRF(qrf_fit, R=10, response_members="all",
                                              seed=42, sort=TRUE, verbose=TRUE)
# saveRDS(cd4_qrf_vimp, file=file.path(app_path, "results", "cd4_qrf_vimp.rds"))
cd4_qrf_vimp_summary <- cd4_qrf_vimp$summary
write.csv(cd4_qrf_vimp_summary, file=file.path(app_path, "results", "cd4_qrf_vimp_summary.csv"), row.names=FALSE)
print(cd4_qrf_vimp_summary)
#        variable baseline_mean_CRPS permuted_mean_CRPS     VI_CRPS  VI_percent
# 1      base_cd4           97.26415          223.62852 126.3643609 129.9187370
# 2       base_VL           97.26415          103.98973   6.7255745   6.9147515
# 3          site           97.26415          102.00125   4.7370932   4.8703381
# 4    prior_aids           97.26415          101.38897   4.1248115   4.2408342
# 5       age_art           97.26415          101.20368   3.9395226   4.0503335
# 6  months_to_vl           97.26415          100.19294   2.9287895   3.0111705
# 7      year_art           97.26415           99.16897   1.9048150   1.9583937
# 8       infmode           97.26415           99.06771   1.8035522   1.8542825
# 9   art_regimen           97.26415           98.35433   1.0901721   1.1208364
# 10       male_y           97.26415           97.65571   0.3915533   0.4025669

set.seed(42)
wrf_fit <- myRF(out_cd4~., data=syn_dat, ntree=1000, mtry=3, nodesize=5,
                split_rule="Wilcoxon", replace=FALSE)
cd4_wrf_vimp <- variable_importance_crps.myRF(wrf_fit, R=10, response_members="all",
                                              seed=42, sort=TRUE, verbose=TRUE)
# saveRDS(cd4_wrf_vimp, file=file.path(app_path, "results", "cd4_wrf_vimp.rds"))
cd4_wrf_vimp_summary <- cd4_wrf_vimp$summary
write.csv(cd4_wrf_vimp_summary, file=file.path(app_path, "results", "cd4_wrf_vimp_summary.csv"), row.names=FALSE)
print(cd4_wrf_vimp_summary)
#        variable baseline_mean_CRPS permuted_mean_CRPS     VI_CRPS  VI_percent
# 1      base_cd4           97.40172          222.87584 125.4741223 128.8212585
# 2       base_VL           97.40172          104.21620   6.8144786   6.9962610
# 3          site           97.40172          101.90048   4.4987615   4.6187701
# 4    prior_aids           97.40172          101.82830   4.4265758   4.5446587
# 5       age_art           97.40172          101.35281   3.9510924   4.0564914
# 6  months_to_vl           97.40172          100.25462   2.8528960   2.9289996
# 7       infmode           97.40172           99.35011   1.9483886   2.0003636
# 8      year_art           97.40172           99.34599   1.9442737   1.9961390
# 9   art_regimen           97.40172           98.43401   1.0322852   1.0598223
# 10       male_y           97.40172           97.86653   0.4648118   0.4772111





##### Figure: Predict ECDF for CD4 #####
library(ggplot2)

find_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

two_profiles <- syn_dat[1:2,]
for (i in seq_len(ncol(two_profiles))) {
  if (is.numeric(dat[,i]) & length(unique(dat[,i]))>2) {
    two_profiles[,i] <- median(syn_dat[,i])
  } else {
    two_profiles[,i] <- find_mode(syn_dat[,i])
  }
}
two_profiles$base_cd4 <- c(200, 500)
two_profiles$base_VL <- c(1e6, 1e4)
write.csv(two_profiles, file=file.path(app_path, "results", "two_profiles.csv"), row.names=FALSE)

cd4_ecdf_wrf <- list()

set.seed(42)
cd4wrf_fit <- myRF(out_cd4~., data=syn_dat, ntree=1000, mtry=3, nodesize=5, split_rule="Wilcoxon")
cd4_ecdf_wrf$subj1 <- predict(cd4wrf_fit, newdata=two_profiles[1,], what="ECDF")
cd4_ecdf_wrf$subj2 <- predict(cd4wrf_fit, newdata=two_profiles[2,], what="ECDF")



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

cdf_pl_d <- make_ecdf_df(cd4_ecdf_wrf)
cdf_pl_d$Subject <- factor(cdf_pl_d$Subject,
                           levels = c("Baseline CD4=200, VL=1e6",
                                      "Baseline CD4=500, VL=1e4"))

fig1_cd4 <- ggplot(cdf_pl_d, aes(x = y, y = ECDF, linetype = Subject)) +
  geom_step(linewidth = 0.7) +
  geom_vline(xintercept = cd4vals, color = "gray60", linetype = "longdash", linewidth = 0.4) +
  scale_x_continuous(breaks = c(0, cd4vals, 1000, 1500)) +
  coord_cartesian(xlim = c(0, 1500)) + 
  labs(x = expression(paste("6-month CD4 (cells/", mu, "L)")),
       y = "Cumulative probability", linetype = NULL) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9),
        legend.key.width = unit(2, "cm"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
  guides(linetype = guide_legend(ncol = 2))
ggsave(filename = file.path(app_path, "figures", "Figure_cd4_ecdf.pdf"),
       plot = fig1_cd4, width = 6, height = 5)



##### Figure: calibration plot #####
make_calibration_data <- function(pred, observed_y, threshold,
                                  method, n_groups = 10L) {
  prediction_order <- order(pred, seq_along(pred))
  group <- integer(length(pred))
  group[prediction_order] <- pmin(n_groups,
                                  ceiling(seq_along(prediction_order) *
                                            n_groups/length(prediction_order)))
  do.call(rbind,
    lapply(sort(unique(group)), function(g) {
      idx <- group == g
      data.frame(Method = method, Threshold = threshold, Group = g,
                 N = sum(idx), Mean_predicted_probability = mean(pred[idx]),
                 Observed_proportion = mean(observed_y[idx] <= threshold))
    }))
}

calibration_list <- list()
counter <- 0
for (j in seq_len(3)) {
  counter <- counter + 1
  calibration_list[[counter]] <- make_calibration_data(
    pred = qrf_prob[, j], observed_y = syn_dat$out_cd4,
    threshold = cd4vals[j], method = "QRF")
  counter <- counter + 1
  calibration_list[[counter]] <- make_calibration_data(
    pred = wrf_prob[, j], observed_y = syn_dat$out_cd4,
    threshold = cd4vals[j], method = "WRF")
}
cd4_calibration <- do.call(rbind, calibration_list)
cd4_calibration$Threshold_label <- factor(
  cd4_calibration$Threshold, levels = cd4vals,
  labels = paste0("P(CD4 <= ", cd4vals, ")"))

write.csv(cd4_calibration,
          file = file.path(app_path, "results", "cd4_calibration.csv"), row.names = FALSE)

calibration_plot <- ggplot(cd4_calibration,
                           aes(x = Mean_predicted_probability,
                               y = Observed_proportion,
                               group = Method, linetype = Method, 
                               shape = Method, colour = Method)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") + geom_line() +
  geom_point(size = 2) + facet_wrap(~ Threshold_label, nrow = 1) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Mean predicted probability", y = "Observed proportion",
       title = "Calibration plots for synthetic 6-month CD4 data") +
  theme_bw() +
  theme(legend.position = "bottom", strip.background = element_blank())
  
ggsave(filename = file.path(app_path, "figures", "Figure_cd4_calibration.pdf"),
       plot = calibration_plot, width = 10, height = 3.8)




