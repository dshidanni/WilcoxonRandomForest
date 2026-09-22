rm(list=ls())
source(file.path("R", "myRF.R"))
source(file.path("R", "CPMtree.R"))

probs <- c(0.1, 0.5, 0.9)
c4 <- 1/7
c11 <- qnorm(1/6)
c12 <- qnorm(7/12)
c3 <- qt(1/5,df=10)
h_func <- function(y) qchisq(pnorm(y), df = 15)

datgen <- function(n) {
  X1 <- rnorm(n, 0, 1)
  X2 <- rbinom(n, 1, 0.5)
  X3 <- rt(n, df = 10)
  X4 <- runif(n, 0, 1)
  X5 <- rbeta(n, 3, 2)
  X6 <- rpois(n, 1)
  eps <- rnorm(n, 0, 1)
  
  group <- rep(7, n)
  group[X1 > c12 & X2 == 1] <- 6
  group[X1 <= c12 & X2 == 0] <- 5
  group[X1 <= c12 & X2 == 1] <- 4
  group[X3 <= c3] <- 3
  group[X1 <= c11] <- 2
  group[X4 <= c4] <- 1
  
  base_Y <- (group-4)/2
  Y <- base_Y + eps
  Ytr <- h_func(Y)
  
  data.frame(X1, X2, X3, X4, X5, X6, Y, Ytr, group)
}

dir_path <- "Simulations"
load(file.path(dir_path, "dat_test.RData"))




##### Scenario 1: Trees #####
ntrains <- c(1000, 2000, 5000)
nsim_tree <- 1000

# convert probability-scale residual to CvM statistic for calibration
PSR_to_CvM <- function(PSR, num) {
  U_vec <- 0.5*(PSR+1)
  U_sorted <- sort(U_vec)
    i <- seq_len(num)
    term_sum <- sum((U_sorted - (2 * i - 1) / (2 * num))^2)
    W2 <- term_sum + 1 / (12 * num)
    return(W2)
}

c_Y <- 0; c_Ytr <- qchisq(pnorm(0), df=15)
nodesizes_tree <- c(60, 120, 300)

for (n in 1:3) {
  nt <- ntrains[n]
  node_size <- nodesizes_tree[n]
  
  sse_qt_10 <- sse_qt_50 <- sse_qt_90 <- sse_prob <- sse_mean <- sse_crps <- sse_psr <- 
    trsse_qt_10 <- trsse_qt_50 <- trsse_qt_90 <- trsse_prob <- trsse_mean <- trsse_crps <- trsse_psr <- 
    wil_qt_10 <- wil_qt_50 <- wil_qt_90 <- wil_prob <- wil_mean <- wil_crps <- wil_psr <- 
    trwil_qt_10 <- trwil_qt_50 <- trwil_qt_90 <- trwil_prob <- trwil_mean <- trwil_crps <- trwil_psr <- 
    pbl_qt_10 <- pbl_qt_50 <- pbl_qt_90 <- pbl_prob <- pbl_mean <- pbl_crps <- pbl_psr <- 
    trpbl_qt_10 <- trpbl_qt_50 <- trpbl_qt_90 <- trpbl_prob <- trpbl_mean <- trpbl_crps <- trpbl_psr <- 
    cpmp_qt_10 <- cpmp_qt_50 <- cpmp_qt_90 <- cpmp_prob <- cpmp_mean <- cpmp_crps <- cpmp_psr <- 
    trcpmp_qt_10 <- trcpmp_qt_50 <- trcpmp_qt_90 <- trcpmp_prob <- trcpmp_mean <- trcpmp_crps <- trcpmp_psr <- 
    cpml_qt_10 <- cpml_qt_50 <- cpml_qt_90 <- cpml_prob <- cpml_mean <- cpml_crps <- cpml_psr <- 
    trcpml_qt_10 <- trcpml_qt_50 <- trcpml_qt_90 <- trcpml_prob <- trcpml_mean <- trcpml_crps <- trcpml_psr <-
    test_Y <- test_Ytr <- data.frame(matrix(0, nrow=nrow(dat_test), ncol=nsim_tree))
  sse_cvm <- sse_time <- trsse_cvm <- trsse_time <- 
    wil_time <- trwil_time <- wil_cvm <- trwil_cvm <- 
    pbl_cvm <- pbl_time <- trpbl_cvm <- trpbl_time <- 
    cpmp_cvm <- cpmp_time <- trcpmp_cvm <- trcpmp_time <- 
    cpml_cvm <- cpml_time <- trcpml_cvm <- trcpml_time <- 
    data.frame(matrix(0, ncol=1, nrow=nsim_tree))
  
  for (sim in 1:nsim_tree) {
    set.seed(42+sim)
    dat_train <- datgen(nt)
    
    Y <- with(dat_test, truemean + rnorm(nrow(dat_test))) # untransformed
    Ytr <- qchisq(pnorm(Y), df = 15) # transformed
    test_Y[, sim] <- Y
    test_Ytr[, sim] <- Ytr

    form_Y <- Y ~ X1 + X2 + X3 + X4 + X5 + X6
    t_Y_sse <- system.time({ m_Y_sse <- myRF(form_Y, data=dat_train, ntree=1, mtry=6, 
                                              replace=FALSE, subsample_prop=1, 
                                             nodesize=node_size, split_rule="SSE") })[["elapsed"]]
    t_Y_wil <- system.time({ m_Y_wil <- myRF(form_Y, data=dat_train, ntree=1, mtry=6, 
                                              replace=FALSE, subsample_prop=1, 
                                             nodesize=node_size, split_rule="Wilcoxon") })[["elapsed"]]
    t_Y_pbl <- system.time({ m_Y_pbl <- myRF(form_Y, data=dat_train, ntree=1, mtry=6, 
                                              replace=FALSE, subsample_prop=1, 
                                             nodesize=node_size, split_rule="pinball", 
                                             pinball_taus=c(0.1,0.5,0.9)) })[["elapsed"]]
    t_Y_cpmp <- system.time({ m_Y_cpmp <- myCPMForest(m_Y_wil, family="probit") })[["elapsed"]]
    t_Y_cpmp <- t_Y_cpmp + t_Y_wil
    t_Y_cpml <- system.time({ m_Y_cpml <- myCPMForest(m_Y_wil, family="logistic") })[["elapsed"]]
    t_Y_cpml <- t_Y_cpml + t_Y_wil
    
    sse_qt <- predict(m_Y_sse, what="quantiles", newdata=dat_test, probs=probs)
    sse_qt_10[,sim] <- sse_qt[,1]; sse_qt_50[,sim] <- sse_qt[,2]; sse_qt_90[,sim] <- sse_qt[,3]
    sse_pr <- predict(m_Y_sse, what="probability", newdata=dat_test, outcomes=c_Y)
    sse_prob[,sim] <- sse_pr[,1]
    sse_mean[,sim] <- predict(m_Y_sse, what="mean", newdata=dat_test)
    sse_crps[,sim] <- predict(m_Y_sse, what="CRPS", newdata=dat_test, y_obsvd=Y)
    sse_psr[,sim] <- predict(m_Y_sse, what="PSR", newdata=dat_test, y_obsvd=Y)
    sse_cvm[sim,1] <- PSR_to_CvM(as.vector(sse_psr[,sim]), 1000)
    sse_time[sim,1] <- t_Y_sse
    print(paste0(nt,", ",sim,", SSE done"))
    
    wil_qt <- predict(m_Y_wil, what="quantiles", newdata=dat_test, probs=probs)
    wil_qt_10[,sim] <- wil_qt[,1]; wil_qt_50[,sim] <- wil_qt[,2]; wil_qt_90[,sim] <- wil_qt[,3]
    wil_pr <- predict(m_Y_wil, what="probability", newdata=dat_test, outcomes=c_Y)
    wil_prob[,sim] <- wil_pr[, 1]
    wil_mean[,sim] <- predict(m_Y_wil, what="mean", newdata=dat_test)
    wil_crps[,sim] <- predict(m_Y_wil, what="CRPS", newdata=dat_test, y_obsvd=Y)
    wil_psr[,sim] <- predict(m_Y_wil, what="PSR", newdata=dat_test, y_obsvd=Y)
    wil_cvm[sim,1] <- PSR_to_CvM(as.vector(wil_psr[,sim]), 1000)
    wil_time[sim,1] <- t_Y_wil
    print(paste0(nt,", ",sim,", Wilcoxon done"))
    
    pbl_qt <- predict(m_Y_pbl, what="quantiles", newdata=dat_test, probs=probs)
    pbl_qt_10[,sim] <- pbl_qt[,1]; pbl_qt_50[,sim] <- pbl_qt[,2]; pbl_qt_90[,sim] <- pbl_qt[,3]
    pbl_pr <- predict(m_Y_pbl, what="probability", newdata=dat_test, outcomes=c_Y)
    pbl_prob[,sim] <- pbl_pr[, 1]
    pbl_mean[,sim] <- predict(m_Y_pbl, what="mean", newdata=dat_test)
    pbl_crps[,sim] <- predict(m_Y_pbl, what="CRPS", newdata=dat_test, y_obsvd=Y)
    pbl_psr[,sim] <- predict(m_Y_pbl, what="PSR", newdata=dat_test, y_obsvd=Y)
    pbl_cvm[sim,1] <- PSR_to_CvM(as.vector(pbl_psr[,sim]), 1000)
    pbl_time[sim,1] <- t_Y_pbl
    print(paste0(nt,", ",sim,", Pinball done"))
    
    cpmp_qt <- predict(m_Y_cpmp, what="quantiles", newdata=dat_test, probs=probs)
    cpmp_qt_10[,sim] <- cpmp_qt[,1]; cpmp_qt_50[,sim] <- cpmp_qt[,2]; cpmp_qt_90[,sim] <- cpmp_qt[,3]
    cpmp_pr <- predict(m_Y_cpmp, newdata=dat_test, what="probability", outcomes=c_Y)
    cpmp_prob[,sim] <- cpmp_pr[, 1]
    cpmp_mean[,sim] <- predict(m_Y_cpmp, what="mean", newdata=dat_test)
    cpmp_crps[,sim] <- predict(m_Y_cpmp, what="CRPS", newdata=dat_test, y_obsvd=Y)
    cpmp_psr[,sim] <- predict(m_Y_cpmp, what="PSR", newdata=dat_test, y_obsvd=Y)
    cpmp_cvm[sim,1] <- PSR_to_CvM(as.vector(cpmp_psr[,sim]), 1000)
    cpmp_time[sim,1] <- t_Y_cpmp
    print(paste0(nt,", ",sim,", CPM-probit done"))
    
    cpml_qt <- predict(m_Y_cpml, what="quantiles", newdata=dat_test, probs=probs)
    cpml_qt_10[,sim] <- cpml_qt[,1]; cpml_qt_50[,sim] <- cpml_qt[,2]; cpml_qt_90[,sim] <- cpml_qt[,3]
    cpml_pr <- predict(m_Y_cpml, newdata=dat_test, what="probability", outcomes=c_Y)
    cpml_prob[,sim] <- cpml_pr[, 1]
    cpml_mean[,sim] <- predict(m_Y_cpml, what="mean", newdata=dat_test)
    cpml_crps[,sim] <- predict(m_Y_cpml, what="CRPS", newdata=dat_test, y_obsvd=Y)
    cpml_psr[,sim] <- predict(m_Y_cpml, what="PSR", newdata=dat_test, y_obsvd=Y)
    cpml_cvm[sim,1] <- PSR_to_CvM(as.vector(cpml_psr[,sim]), 1000)
    cpml_time[sim,1] <- t_Y_cpml
    print(paste0(nt,", ",sim,", CPM-logit done"))
    
    save(sse_qt_10, file=paste0(dir_path,"/s1_trees/",nt,"/sse_qt_10.RData"))
    save(sse_qt_50, file=paste0(dir_path,"/s1_trees/",nt,"/sse_qt_50.RData"))
    save(sse_qt_90, file=paste0(dir_path,"/s1_trees/",nt,"/sse_qt_90.RData"))
    save(sse_prob, file=paste0(dir_path,"/s1_trees/",nt,"/sse_prob.RData"))
    save(sse_mean, file=paste0(dir_path,"/s1_trees/",nt,"/sse_mean.RData"))
    save(sse_crps, file=paste0(dir_path,"/s1_trees/",nt,"/sse_crps.RData"))
    save(sse_psr, file=paste0(dir_path,"/s1_trees/",nt,"/sse_psr.RData"))
    save(sse_cvm, file=paste0(dir_path,"/s1_trees/",nt,"/sse_cvm.RData"))
    save(sse_time, file=paste0(dir_path,"/s1_trees/",nt,"/sse_time.RData"))
    
    save(wil_qt_10, file=paste0(dir_path,"/s1_trees/",nt,"/wil_qt_10.RData"))
    save(wil_qt_50, file=paste0(dir_path,"/s1_trees/",nt,"/wil_qt_50.RData"))
    save(wil_qt_90, file=paste0(dir_path,"/s1_trees/",nt,"/wil_qt_90.RData"))
    save(wil_prob, file=paste0(dir_path,"/s1_trees/",nt,"/wil_prob.RData"))
    save(wil_mean, file=paste0(dir_path,"/s1_trees/",nt,"/wil_mean.RData"))
    save(wil_crps, file=paste0(dir_path,"/s1_trees/",nt,"/wil_crps.RData"))
    save(wil_psr, file=paste0(dir_path,"/s1_trees/",nt,"/wil_psr.RData"))
    save(wil_cvm, file=paste0(dir_path,"/s1_trees/",nt,"/wil_cvm.RData"))
    save(wil_time, file=paste0(dir_path,"/s1_trees/",nt,"/wil_time.RData"))
    
    save(pbl_qt_10, file=paste0(dir_path,"/s1_trees/",nt,"/pbl_qt_10.RData"))
    save(pbl_qt_50, file=paste0(dir_path,"/s1_trees/",nt,"/pbl_qt_50.RData"))
    save(pbl_qt_90, file=paste0(dir_path,"/s1_trees/",nt,"/pbl_qt_90.RData"))
    save(pbl_prob, file=paste0(dir_path,"/s1_trees/",nt,"/pbl_prob.RData"))
    save(pbl_mean, file=paste0(dir_path,"/s1_trees/",nt,"/pbl_mean.RData"))
    save(pbl_crps, file=paste0(dir_path,"/s1_trees/",nt,"/pbl_crps.RData"))
    save(pbl_psr, file=paste0(dir_path,"/s1_trees/",nt,"/pbl_psr.RData"))
    save(pbl_cvm, file=paste0(dir_path,"/s1_trees/",nt,"/pbl_cvm.RData"))
    save(pbl_time, file=paste0(dir_path,"/s1_trees/",nt,"/pbl_time.RData"))
    
    save(cpmp_qt_10, file=paste0(dir_path,"/s1_trees/",nt,"/cpmp_qt_10.RData"))
    save(cpmp_qt_50, file=paste0(dir_path,"/s1_trees/",nt,"/cpmp_qt_50.RData"))
    save(cpmp_qt_90, file=paste0(dir_path,"/s1_trees/",nt,"/cpmp_qt_90.RData"))
    save(cpmp_prob, file=paste0(dir_path,"/s1_trees/",nt,"/cpmp_prob.RData"))
    save(cpmp_mean, file=paste0(dir_path,"/s1_trees/",nt,"/cpmp_mean.RData"))
    save(cpmp_crps, file=paste0(dir_path,"/s1_trees/",nt,"/cpmp_crps.RData"))
    save(cpmp_psr, file=paste0(dir_path,"/s1_trees/",nt,"/cpmp_psr.RData"))
    save(cpmp_cvm, file=paste0(dir_path,"/s1_trees/",nt,"/cpmp_cvm.RData"))
    save(cpmp_time, file=paste0(dir_path,"/s1_trees/",nt,"/cpmp_time.RData"))
    
    save(cpml_qt_10, file=paste0(dir_path,"/s1_trees/",nt,"/cpml_qt_10.RData"))
    save(cpml_qt_50, file=paste0(dir_path,"/s1_trees/",nt,"/cpml_qt_50.RData"))
    save(cpml_qt_90, file=paste0(dir_path,"/s1_trees/",nt,"/cpml_qt_90.RData"))
    save(cpml_prob, file=paste0(dir_path,"/s1_trees/",nt,"/cpml_prob.RData"))
    save(cpml_mean, file=paste0(dir_path,"/s1_trees/",nt,"/cpml_mean.RData"))
    save(cpml_crps, file=paste0(dir_path,"/s1_trees/",nt,"/cpml_crps.RData"))
    save(cpml_psr, file=paste0(dir_path,"/s1_trees/",nt,"/cpml_psr.RData"))
    save(cpml_cvm, file=paste0(dir_path,"/s1_trees/",nt,"/cpml_cvm.RData"))
    save(cpml_time, file=paste0(dir_path,"/s1_trees/",nt,"/cpml_time.RData"))
    
    save(test_Y, file=paste0(dir_path,"/s1_trees/",nt,"/test_Y.RData"))
    save(test_Ytr, file=paste0(dir_path,"/s1_trees/",nt,"/test_Ytr.RData"))
    
    
    # transformed Y
    form_Ytr <- Ytr ~ X1 + X2 + X3 + X4 + X5 + X6
    t_Ytr_sse <- system.time({ m_Ytr_sse <- myRF(form_Ytr, data=dat_train, ntree=1, mtry=6, 
                                              replace=FALSE, subsample_prop=1, 
                                              nodesize=node_size, split_rule="SSE") })[["elapsed"]]
    t_Ytr_wil <- system.time({ m_Ytr_wil <- myRF(form_Ytr, data=dat_train, ntree=1, mtry=6, 
                                              replace=FALSE, subsample_prop=1, 
                                              nodesize=node_size, split_rule="Wilcoxon") })[["elapsed"]]
    t_Ytr_pbl <- system.time({ m_Ytr_pbl <- myRF(form_Ytr, data=dat_train, ntree=1, mtry=6, 
                                              replace=FALSE, subsample_prop=1, 
                                              nodesize=node_size, split_rule="pinball", 
                                              pinball_taus=c(0.1,0.5,0.9)) })[["elapsed"]]
    t_Ytr_cpmp <- system.time({ m_Ytr_cpmp <- myCPMForest(m_Ytr_wil, family="probit") })[["elapsed"]]
    t_Ytr_cpmp <- t_Ytr_cpmp + t_Ytr_wil
    t_Ytr_cpml <- system.time({ m_Ytr_cpml <- myCPMForest(m_Ytr_wil, family="logistic") })[["elapsed"]]
    t_Ytr_cpml <- t_Ytr_cpml + t_Ytr_wil
    
    trsse_qt <- predict(m_Ytr_sse, what="quantiles", newdata=dat_test, probs=probs)
    trsse_qt_10[,sim] <- trsse_qt[,1]; trsse_qt_50[,sim] <- trsse_qt[,2]; trsse_qt_90[,sim] <- trsse_qt[,3]
    trsse_pr <- predict(m_Ytr_sse, what="probability", newdata=dat_test, outcomes=c_Ytr)
    trsse_prob[,sim] <- trsse_pr[,1]
    trsse_mean[,sim] <- predict(m_Ytr_sse, what="mean", newdata=dat_test)
    trsse_crps[,sim] <- predict(m_Ytr_sse, what="CRPS", newdata=dat_test, y_obsvd=Ytr)
    trsse_psr[,sim] <- predict(m_Ytr_sse, what="PSR", newdata=dat_test, y_obsvd=Ytr)
    trsse_cvm[sim,1] <- PSR_to_CvM(as.vector(trsse_psr[,sim]), 1000)
    trsse_time[sim,1] <- t_Ytr_sse

    trwil_qt <- predict(m_Ytr_wil, what="quantiles", newdata=dat_test, probs=probs)
    trwil_qt_10[,sim] <- trwil_qt[,1]; trwil_qt_50[,sim] <- trwil_qt[,2]; trwil_qt_90[,sim] <- trwil_qt[,3]
    trwil_pr <- predict(m_Ytr_wil, what="probability", newdata=dat_test, outcomes=c_Ytr)
    trwil_prob[,sim] <- trwil_pr[,1]
    trwil_mean[,sim] <- predict(m_Ytr_wil, what="mean", newdata=dat_test)
    trwil_crps[,sim] <- predict(m_Ytr_wil, what="CRPS", newdata=dat_test, y_obsvd=Ytr)
    trwil_psr[,sim] <- predict(m_Ytr_wil, what="PSR", newdata=dat_test, y_obsvd=Ytr)
    trwil_cvm[sim,1] <- PSR_to_CvM(as.vector(trwil_psr[,sim]), 1000)
    trwil_time[sim,1] <- t_Ytr_wil

    trpbl_qt <- predict(m_Ytr_pbl, what="quantiles", newdata=dat_test, probs=probs)
    trpbl_qt_10[,sim] <- trpbl_qt[,1]; trpbl_qt_50[,sim] <- trpbl_qt[,2]; trpbl_qt_90[,sim] <- trpbl_qt[,3]
    trpbl_pr <- predict(m_Ytr_pbl, what="probability", newdata=dat_test, outcomes=c_Ytr)
    trpbl_prob[,sim] <- trpbl_pr[,1]
    trpbl_mean[,sim] <- predict(m_Ytr_pbl, what="mean", newdata=dat_test)
    trpbl_crps[,sim] <- predict(m_Ytr_pbl, what="CRPS", newdata=dat_test, y_obsvd=Ytr)
    trpbl_psr[,sim] <- predict(m_Ytr_pbl, what="PSR", newdata=dat_test, y_obsvd=Ytr)
    trpbl_cvm[sim,1] <- PSR_to_CvM(as.vector(trpbl_psr[,sim]), 1000)
    trpbl_time[sim,1] <- t_Ytr_pbl

    trcpmp_qt <- predict(m_Ytr_cpmp, what="quantiles", newdata=dat_test, probs=probs)
    trcpmp_qt_10[,sim] <- trcpmp_qt[,1]; trcpmp_qt_50[,sim] <- trcpmp_qt[,2]; trcpmp_qt_90[,sim] <- trcpmp_qt[,3]
    trcpmp_pr <- predict(m_Ytr_cpmp, newdata=dat_test, what="probability", outcomes=c_Ytr)
    trcpmp_prob[,sim] <- trcpmp_pr[,1]
    trcpmp_mean[,sim] <- predict(m_Ytr_cpmp, what="mean", newdata=dat_test)
    trcpmp_crps[,sim] <- predict(m_Ytr_cpmp, what="CRPS", newdata=dat_test, y_obsvd=Ytr)
    trcpmp_psr[,sim] <- predict(m_Ytr_cpmp, what="PSR", newdata=dat_test, y_obsvd=Ytr)
    trcpmp_cvm[sim,1] <- PSR_to_CvM(as.vector(trcpmp_psr[,sim]), 1000)
    trcpmp_time[sim,1] <- t_Ytr_cpmp

    trcpml_qt <- predict(m_Ytr_cpml, what="quantiles", newdata=dat_test, probs=probs)
    trcpml_qt_10[,sim] <- trcpml_qt[,1]; trcpml_qt_50[,sim] <- trcpml_qt[,2]; trcpml_qt_90[,sim] <- trcpml_qt[,3]
    trcpml_pr <- predict(m_Ytr_cpml, newdata=dat_test, what="probability", outcomes=c_Ytr)
    trcpml_prob[,sim] <- trcpml_pr[,1]
    trcpml_mean[,sim] <- predict(m_Ytr_cpml, what="mean", newdata=dat_test)
    trcpml_crps[,sim] <- predict(m_Ytr_cpml, what="CRPS", newdata=dat_test, y_obsvd=Ytr)
    trcpml_psr[,sim] <- predict(m_Ytr_cpml, what="PSR", newdata=dat_test, y_obsvd=Ytr)
    trcpml_cvm[sim,1] <- PSR_to_CvM(as.vector(trcpml_psr[,sim]), 1000)
    trcpml_time[sim,1] <- t_Ytr_cpml

    save(trsse_qt_10, file=paste0(dir_path,"/s1_trees/",nt,"/trsse_qt_10.RData"))
    save(trsse_qt_50, file=paste0(dir_path,"/s1_trees/",nt,"/trsse_qt_50.RData"))
    save(trsse_qt_90, file=paste0(dir_path,"/s1_trees/",nt,"/trsse_qt_90.RData"))
    save(trsse_prob, file=paste0(dir_path,"/s1_trees/",nt,"/trsse_prob.RData"))
    save(trsse_mean, file=paste0(dir_path,"/s1_trees/",nt,"/trsse_mean.RData"))
    save(trsse_time, file=paste0(dir_path,"/s1_trees/",nt,"/trsse_time.RData"))
    save(trsse_crps, file=paste0(dir_path,"/s1_trees/",nt,"/trsse_crps.RData"))
    save(trsse_psr, file=paste0(dir_path,"/s1_trees/",nt,"/trsse_psr.RData"))
    save(trsse_cvm, file=paste0(dir_path,"/s1_trees/",nt,"/trsse_cvm.RData"))
    
    save(trwil_qt_10, file=paste0(dir_path,"/s1_trees/",nt,"/trwil_qt_10.RData"))
    save(trwil_qt_50, file=paste0(dir_path,"/s1_trees/",nt,"/trwil_qt_50.RData"))
    save(trwil_qt_90, file=paste0(dir_path,"/s1_trees/",nt,"/trwil_qt_90.RData"))
    save(trwil_prob, file=paste0(dir_path,"/s1_trees/",nt,"/trwil_prob.RData"))
    save(trwil_mean, file=paste0(dir_path,"/s1_trees/",nt,"/trwil_mean.RData"))
    save(trwil_time, file=paste0(dir_path,"/s1_trees/",nt,"/trwil_time.RData"))
    save(trwil_crps, file=paste0(dir_path,"/s1_trees/",nt,"/trwil_crps.RData"))
    save(trwil_psr, file=paste0(dir_path,"/s1_trees/",nt,"/trwil_psr.RData"))
    save(trwil_cvm, file=paste0(dir_path,"/s1_trees/",nt,"/trwil_cvm.RData"))
    
    save(trpbl_qt_10, file=paste0(dir_path,"/s1_trees/",nt,"/trpbl_qt_10.RData"))
    save(trpbl_qt_50, file=paste0(dir_path,"/s1_trees/",nt,"/trpbl_qt_50.RData"))
    save(trpbl_qt_90, file=paste0(dir_path,"/s1_trees/",nt,"/trpbl_qt_90.RData"))
    save(trpbl_prob, file=paste0(dir_path,"/s1_trees/",nt,"/trpbl_prob.RData"))
    save(trpbl_mean, file=paste0(dir_path,"/s1_trees/",nt,"/trpbl_mean.RData"))
    save(trpbl_time, file=paste0(dir_path,"/s1_trees/",nt,"/trpbl_time.RData"))
    save(trpbl_crps, file=paste0(dir_path,"/s1_trees/",nt,"/trpbl_crps.RData"))
    save(trpbl_psr, file=paste0(dir_path,"/s1_trees/",nt,"/trpbl_psr.RData"))
    save(trpbl_cvm, file=paste0(dir_path,"/s1_trees/",nt,"/trpbl_cvm.RData"))
    
    save(trcpmp_qt_10, file=paste0(dir_path,"/s1_trees/",nt,"/trcpmp_qt_10.RData"))
    save(trcpmp_qt_50, file=paste0(dir_path,"/s1_trees/",nt,"/trcpmp_qt_50.RData"))
    save(trcpmp_qt_90, file=paste0(dir_path,"/s1_trees/",nt,"/trcpmp_qt_90.RData"))
    save(trcpmp_prob, file=paste0(dir_path,"/s1_trees/",nt,"/trcpmp_prob.RData"))
    save(trcpmp_mean, file=paste0(dir_path,"/s1_trees/",nt,"/trcpmp_mean.RData"))
    save(trcpmp_time, file=paste0(dir_path,"/s1_trees/",nt,"/trcpmp_time.RData"))
    save(trcpmp_crps, file=paste0(dir_path,"/s1_trees/",nt,"/trcpmp_crps.RData"))
    save(trcpmp_psr, file=paste0(dir_path,"/s1_trees/",nt,"/trcpmp_psr.RData"))
    save(trcpmp_cvm, file=paste0(dir_path,"/s1_trees/",nt,"/trcpmp_cvm.RData"))
    
    save(trcpml_qt_10, file=paste0(dir_path,"/s1_trees/",nt,"/trcpml_qt_10.RData"))
    save(trcpml_qt_50, file=paste0(dir_path,"/s1_trees/",nt,"/trcpml_qt_50.RData"))
    save(trcpml_qt_90, file=paste0(dir_path,"/s1_trees/",nt,"/trcpml_qt_90.RData"))
    save(trcpml_prob, file=paste0(dir_path,"/s1_trees/",nt,"/trcpml_prob.RData"))
    save(trcpml_mean, file=paste0(dir_path,"/s1_trees/",nt,"/trcpml_mean.RData"))
    save(trcpml_time, file=paste0(dir_path,"/s1_trees/",nt,"/trcpml_time.RData"))
    save(trcpml_crps, file=paste0(dir_path,"/s1_trees/",nt,"/trcpml_crps.RData"))
    save(trcpml_psr, file=paste0(dir_path,"/s1_trees/",nt,"/trcpml_psr.RData"))
    save(trcpml_cvm, file=paste0(dir_path,"/s1_trees/",nt,"/trcpml_cvm.RData"))
  }
}







##### Scenario 1: Forests #####
nsim_forest <- 250
forest_path <- file.path(dir_path, "s1_forests")
form_Y <- Y ~ X1 + X2 + X3 + X4 + X5 + X6
form_Ytr <- Ytr ~ X1 + X2 + X3 + X4 + X5 + X6
nodesizes1 <- c(10, 20, 50)
test_Y <- load(file.path(dir_path, ))

for (n in 1:3) {
  nt <- ntrains[n]
  sse_qt_10 <- sse_qt_50 <- sse_qt_90 <- sse_prob <- sse_mean <- sse_crps <- sse_psr <- 
    trsse_qt_10 <- trsse_qt_50 <- trsse_qt_90 <- trsse_prob <- trsse_mean <- trsse_crps <- trsse_psr <- 
    wil_qt_10 <- wil_qt_50 <- wil_qt_90 <- wil_prob <- wil_mean <- wil_crps <- wil_psr <- 
    trwil_qt_10 <- trwil_qt_50 <- trwil_qt_90 <- trwil_prob <- trwil_mean <- trwil_crps <- trwil_psr <- 
    test_Y <- test_Ytr <- data.frame(matrix(0, nrow=nrow(dat_test), ncol=nsim_forest))
  sse_cvm <- sse_time <- trsse_cvm <- trsse_time <- 
    wil_time <- trwil_time <- wil_cvm <- trwil_cvm <- 
    data.frame(matrix(0, ncol=1, nrow=nsim_forest))
  node_size <- nodesizes1[n]
  
  for (sim in 1:nsim_forest) {
    set.seed(42+sim)
    dat_train <- datgen(nt)
    Y <- with(dat_test, truemean + rnorm(nrow(dat_test)))
    Ytr <- qchisq(pnorm(Y), df = 15) 
    test_Y[, sim] <- Y
    test_Ytr[, sim] <- Ytr
    
    set.seed(42+sim)
    t_Y_sse <- system.time({ m_Y_sse <- myRF(form_Y, data=dat_train, ntree=1000, mtry=2, 
                                             replace=FALSE, nodesize=node_size, 
                                             split_rule="SSE") })[["elapsed"]]
    set.seed(42+sim)
    t_Y_wil <- system.time({ m_Y_wil <- myRF(form_Y, data=dat_train, ntree=1000, mtry=2, 
                                             replace=FALSE, nodesize=node_size,
                                             split_rule="Wilcoxon") })[["elapsed"]]
    
    sse_qt <- predict(m_Y_sse, what="quantiles", newdata=dat_test, probs=probs)
    sse_qt_10[,sim] <- sse_qt[,1]; sse_qt_50[,sim] <- sse_qt[,2]; sse_qt_90[,sim] <- sse_qt[,3]
    sse_pr <- predict(m_Y_sse, what="probability", newdata=dat_test, outcomes=c_Y)
    sse_prob[,sim] <- sse_pr[,1]
    sse_mean[,sim] <- predict(m_Y_sse, what="mean", newdata=dat_test)
    sse_crps[,sim] <- predict(m_Y_sse, what="CRPS", newdata=dat_test, y_obsvd=Y)
    sse_psr[,sim] <- predict(m_Y_sse, what="PSR", newdata=dat_test, y_obsvd=Y)
    sse_cvm[sim,1] <- PSR_to_CvM(as.vector(sse_psr[,sim]), 1000)
    sse_time[sim,1] <- t_Y_sse
    print(paste0(nt,", ",sim,", SSE done"))
    
    wil_qt <- predict(m_Y_wil, what="quantiles", newdata=dat_test, probs=probs)
    wil_qt_10[,sim] <- wil_qt[,1]; wil_qt_50[,sim] <- wil_qt[,2]; wil_qt_90[,sim] <- wil_qt[,3]
    wil_pr <- predict(m_Y_wil, what="probability", newdata=dat_test, outcomes=c_Y)
    wil_prob[,sim] <- wil_pr[,1]
    wil_mean[,sim] <- predict(m_Y_wil, what="mean", newdata=dat_test)
    wil_crps[,sim] <- predict(m_Y_wil, what="CRPS", newdata=dat_test, y_obsvd=Y)
    wil_psr[,sim] <- predict(m_Y_wil, what="PSR", newdata=dat_test, y_obsvd=Y)
    wil_cvm[sim,1] <- PSR_to_CvM(as.vector(wil_psr[,sim]), 1000)
    wil_time[sim,1] <- t_Y_wil
    print(paste0(nt,", ",sim,", Wilcoxon done"))
    
    save(sse_qt_10, file=paste0(forest_path,"/",nt,"/sse_qt_10.RData"))
    save(sse_qt_50, file=paste0(forest_path,"/",nt,"/sse_qt_50.RData"))
    save(sse_qt_90, file=paste0(forest_path,"/",nt,"/sse_qt_90.RData"))
    save(sse_prob, file=paste0(forest_path,"/",nt,"/sse_prob.RData"))
    save(sse_mean, file=paste0(forest_path,"/",nt,"/sse_mean.RData"))
    save(sse_crps, file=paste0(forest_path,"/",nt,"/sse_crps.RData"))
    save(sse_psr, file=paste0(forest_path,"/",nt,"/sse_psr.RData"))
    save(sse_cvm, file=paste0(forest_path,"/",nt,"/sse_cvm.RData"))
    save(sse_time, file=paste0(forest_path,"/",nt,"/sse_time.RData"))
    
    save(wil_qt_10, file=paste0(forest_path,"/",nt,"/wil_qt_10.RData"))
    save(wil_qt_50, file=paste0(forest_path,"/",nt,"/wil_qt_50.RData"))
    save(wil_qt_90, file=paste0(forest_path,"/",nt,"/wil_qt_90.RData"))
    save(wil_prob, file=paste0(forest_path,"/",nt,"/wil_prob.RData"))
    save(wil_mean, file=paste0(forest_path,"/",nt,"/wil_mean.RData"))
    save(wil_crps, file=paste0(forest_path,"/",nt,"/wil_crps.RData"))
    save(wil_psr, file=paste0(forest_path,"/",nt,"/wil_psr.RData"))
    save(wil_cvm, file=paste0(forest_path,"/",nt,"/wil_cvm.RData"))
    save(wil_time, file=paste0(forest_path,"/",nt,"/wil_time.RData"))
    
    save(test_Y, file=paste0(forest_path,"/",nt,"/test_Y.RData"))
    save(test_Ytr, file=paste0(forest_path,"/",nt,"/test_Ytr.RData"))
    
    
    set.seed(42+sim)
    t_Ytr_sse <- system.time({ m_Ytr_sse <- myRF(form_Ytr, data=dat_train, ntree=1000, mtry=2, 
                                                  replace=FALSE, nodesize=node_size, split_rule="SSE") })[["elapsed"]]
    set.seed(42+sim)
    t_Ytr_wil <- system.time({ m_Ytr_wil <- myRF(form_Ytr, data=dat_train, ntree=1000, mtry=2, 
                                                  replace=FALSE, nodesize=node_size, split_rule="Wilcoxon") })[["elapsed"]]

    trsse_qt <- predict(m_Ytr_sse, what="quantiles", newdata=dat_test, probs=probs)
    trsse_qt_10[,sim] <- trsse_qt[,1]; trsse_qt_50[,sim] <- trsse_qt[,2]; trsse_qt_90[,sim] <- trsse_qt[,3]
    trsse_pr <- predict(m_Ytr_sse, what="probability", newdata=dat_test, outcomes=c_Ytr)
    trsse_prob[,sim] <- trsse_pr[, 1]
    trsse_mean[,sim] <- predict(m_Ytr_sse, what="mean", newdata=dat_test)
    trsse_crps[,sim] <- predict(m_Ytr_sse, what="CRPS", newdata=dat_test, y_obsvd=Ytr)
    trsse_psr[,sim] <- predict(m_Ytr_sse, what="PSR", newdata=dat_test, y_obsvd=Ytr)
    trsse_cvm[sim,1] <- PSR_to_CvM(as.vector(trsse_psr[,sim]), 1000)
    trsse_time[sim,1] <- t_Ytr_sse
    print(paste0(nt,", ",sim,", trSSE done"))
    
    trwil_qt <- predict(m_Ytr_wil, what="quantiles", newdata=dat_test, probs=probs)
    trwil_qt_10[,sim] <- trwil_qt[,1]; trwil_qt_50[,sim] <- trwil_qt[,2]; trwil_qt_90[,sim] <- trwil_qt[,3]
    trwil_pr <- predict(m_Ytr_wil, what="probability", newdata=dat_test, outcomes=c_Ytr)
    trwil_prob[,sim] <- trwil_pr[, 1]
    trwil_mean[,sim] <- predict(m_Ytr_wil, what="mean", newdata=dat_test)
    trwil_crps[,sim] <- predict(m_Ytr_wil, what="CRPS", newdata=dat_test, y_obsvd=Ytr)
    trwil_psr[,sim] <- predict(m_Ytr_wil, what="PSR", newdata=dat_test, y_obsvd=Ytr)
    trwil_cvm[sim,1] <- PSR_to_CvM(as.vector(trwil_psr[,sim]), 1000)
    trwil_time[sim,1] <- t_Ytr_wil
    print(paste0(nt,", ",sim,", trWilcoxon done"))
    
    save(trsse_qt_10, file=paste0(forest_path,"/",nt,"/trsse_qt_10.RData"))
    save(trsse_qt_50, file=paste0(forest_path,"/",nt,"/trsse_qt_50.RData"))
    save(trsse_qt_90, file=paste0(forest_path,"/",nt,"/trsse_qt_90.RData"))
    save(trsse_prob, file=paste0(forest_path,"/",nt,"/trsse_prob.RData"))
    save(trsse_mean, file=paste0(forest_path,"/",nt,"/trsse_mean.RData"))
    save(trsse_time, file=paste0(forest_path,"/",nt,"/trsse_time.RData"))
    save(trsse_crps, file=paste0(forest_path,"/",nt,"/trsse_crps.RData"))
    save(trsse_psr, file=paste0(forest_path,"/",nt,"/trsse_psr.RData"))
    save(trsse_cvm, file=paste0(forest_path,"/",nt,"/trsse_cvm.RData"))
    
    save(trwil_qt_10, file=paste0(forest_path,"/",nt,"/trwil_qt_10.RData"))
    save(trwil_qt_50, file=paste0(forest_path,"/",nt,"/trwil_qt_50.RData"))
    save(trwil_qt_90, file=paste0(forest_path,"/",nt,"/trwil_qt_90.RData"))
    save(trwil_prob, file=paste0(forest_path,"/",nt,"/trwil_prob.RData"))
    save(trwil_mean, file=paste0(forest_path,"/",nt,"/trwil_mean.RData"))
    save(trwil_time, file=paste0(forest_path,"/",nt,"/trwil_time.RData"))
    save(trwil_crps, file=paste0(forest_path,"/",nt,"/trwil_crps.RData"))
    save(trwil_psr, file=paste0(forest_path,"/",nt,"/trwil_psr.RData"))
    save(trwil_cvm, file=paste0(forest_path,"/",nt,"/trwil_cvm.RData"))
  }
}



##### Scenario 1: CPM Forests #####
# CPM-refined WRF
nt <- 2000
node_size <- 20
cpmf_path <- file.path(dir_path, "s1_CPMforests")

cpmp_qt_10 <- cpmp_qt_50 <- cpmp_qt_90 <- cpmp_prob <- cpmp_mean <- cpmp_crps <- cpmp_psr <- 
  trcpmp_qt_10 <- trcpmp_qt_50 <- trcpmp_qt_90 <- trcpmp_prob <- trcpmp_mean <- trcpmp_crps <- trcpmp_psr <- 
  data.frame(matrix(0, nrow=1000, ncol=100))
cpmp_cvm <- cpmp_time <- trcpmp_cvm <- trcpmp_time <- data.frame(nrow = 100, ncol = 1)


for (sim in 1:100) {
  set.seed(42+sim)
  dat_train <- datgen(nt)
  Y <- with(dat_test, mu_X+rnorm(1000))
  Ytr <- qchisq(pnorm(Y), df = 15) 
  
  set.seed(42+sim)
  m_Y_wil <- myRF(form_Y, data=dat_train, ntree=1000, mtry=2, 
                  replace=FALSE, nodesize=node_size, split_rule="Wilcox")
  t_Y_cpmp  <- system.time({ m_Y_cpmp <- myCPMForest(m_Y_wil, family="probit") })[["elapsed"]]
  cpmp_qt <- predict(m_Y_cpmp, what="quantiles", newdata=dat_test, probs=probs)
  cpmp_qt_10[,sim] <- cpmp_qt[,1]; cpmp_qt_50[,sim] <- cpmp_qt[,2]; cpmp_qt_90[,sim] <- cpmp_qt[,3]
  cpmp_pr <- predict(m_Y_cpmp, newdata=dat_test, what="probability", outcomes=c_Y)
  cpmp_prob[,sim] <- cpmp_pr[,1]
  cpmp_mean[,sim] <- predict(m_Y_cpmp, what="mean", newdata=dat_test)
  cpmp_crps[,sim] <- predict(m_Y_cpmp, what="CRPS", newdata=dat_test, y_obsvd=Y)
  cpmp_psr[,sim] <- predict(m_Y_cpmp, what="PSR", newdata=dat_test, y_obsvd=Y)
  cpmp_cvm[sim,1] <- PSR_to_CvM(as.vector(cpmp_psr[,sim]), 1000)
  cpmp_time[sim,1] <- t_Y_cpmp
  
  save(cpmp_qt_10, file=paste0(cpmf_path,"/cpmp_qt_10.RData"))
  save(cpmp_qt_50, file=paste0(cpmf_path,"/cpmp_qt_50.RData"))
  save(cpmp_qt_90, file=paste0(cpmf_path,"/cpmp_qt_90.RData"))
  save(cpmp_prob, file=paste0(cpmf_path,"/cpmp_prob.RData"))
  save(cpmp_mean, file=paste0(cpmf_path,"/cpmp_mean.RData"))
  save(cpmp_crps, file=paste0(cpmf_path,"/cpmp_crps.RData"))
  save(cpmp_psr, file=paste0(cpmf_path,"/cpmp_psr.RData"))
  save(cpmp_cvm, file=paste0(cpmf_path,"/cpmp_cvm.RData"))
  save(cpmp_time, file=paste0(cpmf_path,"/cpmp_time.RData"))
  
  
  m_Ytr_wil <- myRF(form_Ytr, data=dat_train, ntree=1000, mtry=2, 
                    replace=FALSE, nodesize=node_size, split_rule="Wilcox")
  t_Ytr_cpmp  <- system.time({ m_Ytr_cpmp <- myCPMForest(m_Ytr_wil, family="probit") })[["elapsed"]]
  trcpmp_qt <- predict(m_Ytr_cpmp, what="quantiles", newdata=dat_test, probs=probs)
  trcpmp_qt_10[,sim] <- trcpmp_qt[,1]; trcpmp_qt_50[,sim] <- trcpmp_qt[,2]; trcpmp_qt_90[,sim] <- trcpmp_qt[,3]
  trcpmp_pr <- predict(m_Ytr_cpmp, newdata=dat_test, what="probability", outcomes=c_Ytr)
  trcpmp_prob[,sim] <- unlist(trcpmp_pr[1, ])
  trcpmp_mean[,sim] <- predict(m_Ytr_cpmp, what="mean", newdata=dat_test)
  trcpmp_crps[,sim] <- predict(m_Ytr_cpmp, what="CRPS", newdata=dat_test, y_obsvd=Ytr)
  trcpmp_psr[,sim] <- predict(m_Ytr_cpmp, what="PSR", newdata=dat_test, y_obsvd=Ytr)
  trcpmp_cvm[sim,1] <- PSR_to_CvM(as.vector(trcpmp_psr[,sim]), 1000)
  trcpmp_time[sim,1] <- t_Ytr_cpmp
  save(trcpmp_qt_10, file=paste0(cpmf_path,"/trcpmp_qt_10.RData"))
  save(trcpmp_qt_50, file=paste0(cpmf_path,"/trcpmp_qt_50.RData"))
  save(trcpmp_qt_90, file=paste0(cpmf_path,"/trcpmp_qt_90.RData"))
  save(trcpmp_prob, file=paste0(cpmf_path,"/trcpmp_prob.RData"))
  save(trcpmp_mean, file=paste0(cpmf_path,"/trcpmp_mean.RData"))
  save(trcpmp_time, file=paste0(cpmf_path,"/trcpmp_time.RData"))
  save(trcpmp_crps, file=paste0(cpmf_path,"/trcpmp_crps.RData"))
  save(trcpmp_psr, file=paste0(cpmf_path,"/trcpmp_psr.RData"))
  save(trcpmp_cvm, file=paste0(cpmf_path,"/trcpmp_cvm.RData"))
  
}




