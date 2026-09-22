rm(list=ls())

##### Scenario 1 #####
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

# generate the fixed test data
set.seed(42)
dat_test <- datgen(1000)
dat_test <- dat_test[, -c(7,8)] # remove the generated outcomes
# table(dat_test$group)

# calculate the true quantiles / probabilities / means
h_inv_func <- function(c) qnorm(pchisq(c, df = 15))


group_mus <- c("1" = -1.5, "2" = -1, "3" = -0.5, 
               "4" =  0, "5" =  0.5, "6" =  1, "7" =  1.5)

set.seed(42)
large_eps <- rnorm(1e7) # run large simulations to approx means for h(Y)
# estimate conditional means
true_mean_hY_dict <- sapply(group_mus, function(mu) {
  mean(h_func(mu + large_eps))
})
#        1         2         3         4         5         6         7 
# 8.310363 10.240903 12.465696 15.002845 17.868951 21.079136 24.647112 

library(dplyr)
mu_X = group_mus[as.character(dat_test$group)]
dat_test <- dat_test %>%
  mutate(
    truemean = mu_X,
    true10 = qnorm(0.1, mean = mu_X, sd = 1),
    true50 = qnorm(0.5, mean = mu_X, sd = 1),
    true90 = qnorm(0.9, mean = mu_X, sd = 1),
    trueprob = pnorm(0, mean = mu_X, sd = 1),
    trtruemean = true_mean_hY_dict[as.character(group)],
    trtrue10 = h_func(true10),
    trtrue50 = h_func(true50),
    trtrue90 = h_func(true90))

dir_path <- "Simulations"
save(dat_test, file=paste0(dir_path, "/dat_test.RData"))
write.csv(dat_test, file=paste0(dir_path, "/dat_test.csv"), row.names = FALSE)


##### Scenario 2 #####
load(file.path(dir_path ,"dat_test.RData"))
dat_test2 <- dat_test[,1:6]
set.seed(42)
Z <- rchisq(1000, df=3)
eps <- Z-3
mu <- with(dat_test2, 0.5 + X1 + 0.8 * X2 + 0.3 * X3 + 0.7 * sin(2 * pi * X4)-0.5*X1*X2)
sigma <- with(dat_test2, 0.4 + 0.3 * (X4 > 0.5) + 0.5 * (X1 * X3 > 0))
dat_test2$Y <- mu + sigma * eps
dat_test2$truemean <- mu

calc_true_qt <- function(tau, mu, sigma) {
  q_eps <- qchisq(tau, df = 3) - 3
  return(mu + sigma * q_eps)
}
dat_test2$true10 <- calc_true_qt(0.1,mu,sigma)
dat_test2$true50 <- calc_true_qt(0.5,mu,sigma)
dat_test2$true90 <- calc_true_qt(0.9,mu,sigma)
save(dat_test2, file=paste0(dir_path, "/dat_test2.RData"))

