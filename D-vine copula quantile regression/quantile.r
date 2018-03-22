rm(list=ls())

library(MASS)
library(VineCopula)
library(ks)


n <- 500
d <- 3
mean <- rep(0, 4)
cov_mat <- matrix(c(1, 0.4, 0.8, 0,
                    0.4, 1, 0.32, 0,
                    0.8, 0.32, 1, 0, 
                    0, 0, 0, 1),
               nrow=4, ncol=4)
sample <- data.frame(mvrnorm(n, mean, cov_mat))
colnames(sample) <-  c('Y','X1','X2','X3')

v <- kcde(sample$Y, eval.points=sample$Y)$estimate
mu1 <- kcde(sample$X1, eval.points=sample$X1)$estimate
mu2 <- kcde(sample$X2, eval.points=sample$X2)$estimate
mu3 <- kcde(sample$X3, eval.points=sample$X3)$estimate


cop_vx1 <- BiCopSelect(v, mu1, indeptest=TRUE)
cop_vx2 <- BiCopSelect(v, mu2, indeptest=TRUE) 
cop_vx3 <- BiCopSelect(v, mu3, indeptest=TRUE)

cll <- c(cop_vx1$AIC, cop_vx2$AIC, cop_vx3$AIC)
max_cll <- min(cll)
max_idx <- which.min(cll)

cop_u1u2 <- BiCopSelect(mu1, mu2, indeptest=TRUE )


