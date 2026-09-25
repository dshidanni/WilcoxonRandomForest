dir_path <- "Application"
load("ori_dat.RData") # not provided publicly
library(synthpop)
set.seed(42)
syn_obj <- syn(ori_dat, method="cart", k=nrow(ori_dat), seed=42, 
               m=5) # generate 5 synthetic data sets
compare(syn_obj, ori_dat, stat="counts") # compare with marginal distributions

syn_safe <- sdc(syn_obj, ori_dat, rm.replicated.uniques = TRUE)
# no. of replicated uniques removed from each synthetic data set are all 0.
syn_safedat <- syn_safe$syn
write.csv(syn_safedat[[1]], 
          file=file.path(dir_path, "data", "synth_dat1.csv"), row.names=FALSE)
write.csv(syn_safedat[[2]], 
          file=file.path(dir_path, "data", "synth_dat2.csv"), row.names=FALSE)
write.csv(syn_safedat[[3]], 
          file=file.path(dir_path, "data", "synth_dat3.csv"), row.names=FALSE)
write.csv(syn_safedat[[4]], 
          file=file.path(dir_path, "data", "synth_dat4.csv"), row.names=FALSE)
write.csv(syn_safedat[[5]], 
          file=file.path(dir_path, "data", "synth_dat5.csv"), row.names=FALSE)
