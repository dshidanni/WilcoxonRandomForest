rm(list=ls())
source(file.path("R", "myRF.R"))

probs <- c(0.1, 0.5, 0.9)

data_generate <- function (n) {
  X1 <- rnorm(n)
  X2 <- rbinom(n, 1, 0.5)
  X3 <- rt(n, df = 10)
  X4 <- runif(n)
  X5 <- rbeta(n, shape1 = 3, shape2 = 2)
  X6 <- rpois(n, lambda = 1)
  Z <- rchisq(n, df = 3)
  eps <- Z - 3 # centered chi-square (right-skewed)
  mu <- 0.5 + X1 + 0.8 * X2 + 0.3 * X3 + 0.7 * sin(2 * pi * X4)-0.5*X1*X2
  sigma <- 0.4 + 0.3 * (X4 > 0.5) + 0.5 * (X1 * X3 > 0)
  Y <- mu + sigma * eps
  return(data.frame(X1,X2,X3,X4,X5,X6,Y))
}

sim_path <- "Simulations"
load(file.path(sim_path, "dat_test2.RData"))

PSR_to_CvM <- function(PSR, num) {
  U_vec <- 0.5*(PSR+1)
  U_sorted <- sort(U_vec)
  i <- seq_len(num)
  term_sum <- sum((U_sorted - (2 * i - 1) / (2 * num))^2)
  W2 <- term_sum + 1 / (12 * num)
  return(W2)
}

ntrains <- c(1000, 2000, 5000)
s2_path <- paste0(sim_path, "/s2/")
probs <- c(0.1, 0.5, 0.9)
nodesizes <- c(5, 10, 25)

mu_test <- with(dat_test2, 0.5 + X1 + 0.8 * X2 + 0.3 * X3 + 0.7 * sin(2 * pi * X4)-0.5*X1*X2)
sigma_test <- with(dat_test2, 0.4 + 0.3 * (X4 > 0.5) + 0.5 * (X1 * X3 > 0))

form_Y <- Y ~ X1 + X2 + X3 + X4 + X5 + X6

for (n in 1:3) {
  nt <- ntrains[n]
  sse_qt_10 <- sse_qt_50 <- sse_qt_90 <- sse_mean <- 
    wil_qt_10 <- wil_qt_50 <- wil_qt_90 <- wil_mean <- 
    sse_crps <- sse_psr <- wil_crps <- wil_psr <- test_Y <- 
    data.frame(matrix(0, nrow=nrow(dat_test2), ncol=250))
  sse_time <- wil_time <- 
    sse_cvm <- wil_cvm <- data.frame(matrix(0, ncol=1, nrow=250))
  
  node_size <- nodesizes[n]
  for (sim in 1:250) {
    set.seed(42+sim)
    dat_train <- data_generate(nt)
    
    set.seed(42+sim)
    Y <- mu_test + sigma_test*(rchisq(1000,df=3)-3)
    test_Y[,sim] <- Y
    
    set.seed(42+sim)
    t_Y_sse <- system.time({ m_Y_sse <- myRF(form_Y, data=dat_train, ntree=1000, mtry=2, 
                                             replace=FALSE, nodesize=node_size, 
                                             split_rule="SSE") })[["elapsed"]]
    set.seed(42+sim)
    t_Y_wil <- system.time({ m_Y_wil <- myRF(form_Y, data=dat_train, ntree=1000, mtry=2, 
                                             replace=FALSE, nodesize=node_size, 
                                             split_rule="Wilcoxon") })[["elapsed"]]
    
    sse_time[sim,1] <- t_Y_sse
    wil_time[sim,1] <- t_Y_wil
    
    sse_qt <- predict(m_Y_sse, what="quantiles", newdata=dat_test2, probs=probs)
    sse_qt_10[,sim] <- sse_qt[,1]; sse_qt_50[,sim] <- sse_qt[,2]; sse_qt_90[,sim] <- sse_qt[,3]
    sse_mean[,sim] <- predict(m_Y_sse, what="mean", newdata=dat_test2)
    sse_crps[,sim] <- predict(m_Y_sse, what="CRPS", newdata=dat_test2, y_obsvd=Y)
    sse_psr[,sim] <- predict(m_Y_sse, what="PSR", newdata=dat_test2, y_obsvd=Y)
    sse_cvm[sim,1] <- PSR_to_CvM(sse_psr[,sim], 1000)

    wil_qt <- predict(m_Y_wil, what="quantiles", newdata=dat_test2, probs=probs)
    wil_qt_10[,sim] <- wil_qt[,1]; wil_qt_50[,sim] <- wil_qt[,2]; wil_qt_90[,sim] <- wil_qt[,3]
    wil_mean[,sim] <- predict(m_Y_wil, what="mean", newdata=dat_test2)
    wil_crps[,sim] <- predict(m_Y_wil, what="CRPS", newdata=dat_test2, y_obsvd=Y)
    wil_psr[,sim] <- predict(m_Y_wil, what="PSR", newdata=dat_test2, y_obsvd=Y)
    wil_cvm[sim, 1] <- PSR_to_CvM(wil_psr[,sim], 1000)

    save(test_Y, file=paste0(s2_path,nt,"/test_Y.RData"))
    save(sse_qt_10, file=paste0(s2_path,nt,"/sse_qt_10.RData"))
    save(sse_qt_50, file=paste0(s2_path,nt,"/sse_qt_50.RData"))
    save(sse_qt_90, file=paste0(s2_path,nt,"/sse_qt_90.RData"))
    save(sse_mean, file=paste0(s2_path,nt,"/sse_mean.RData"))
    save(sse_crps, file=paste0(s2_path,nt,"/sse_crps.RData"))
    save(sse_psr, file=paste0(s2_path,nt,"/sse_psr.RData"))
    save(sse_cvm, file=paste0(s2_path,nt,"/sse_cvm.RData"))
    save(sse_time, file=paste0(s2_path,nt,"/sse_time.RData"))
    
    save(wil_qt_10, file=paste0(s2_path,nt,"/wil_qt_10.RData"))
    save(wil_qt_50, file=paste0(s2_path,nt,"/wil_qt_50.RData"))
    save(wil_qt_90, file=paste0(s2_path,nt,"/wil_qt_90.RData"))
    save(wil_mean, file=paste0(s2_path,nt,"/wil_mean.RData"))
    save(wil_time, file=paste0(s2_path,nt,"/wil_time.RData"))
    save(wil_crps, file=paste0(s2_path,nt,"/wil_crps.RData"))
    save(wil_psr, file=paste0(s2_path,nt,"/wil_psr.RData"))
    save(wil_cvm, file=paste0(s2_path,nt,"/wil_cvm.RData"))
    
  }
}




