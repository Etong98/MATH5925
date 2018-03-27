rm(list=ls())

library(MASS)
library(VineCopula)
library(copula)
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
u1 <- kcde(sample$X1, eval.points=sample$X1)$estimate
u2 <- kcde(sample$X2, eval.points=sample$X2)$estimate
u3 <- kcde(sample$X3, eval.points=sample$X3)$estimate


cop.vu1 <- BiCopSelect(v, u1, indeptest=TRUE)
cop.vu2 <- BiCopSelect(v, u2, indeptest=TRUE)
cop.vu3 <- BiCopSelect(v, u3, indeptest=TRUE)

cll <- c(cop.vu1$AIC, cop.vu2$AIC, cop.vu3$AIC)
max_cll <- min(cll)
max_idx <- which.min(cll)

cop.u1u2 <- BiCopSelect(u1, u2, indeptest=TRUE )
F1 <- BiCopHfunc2(v, u2, obj=cop.vu2)
F2 <- BiCopHfunc2(u1, u2, obj=cop.u1u2)
cop.vu1_u2 <- BiCopSelect(F1, F2, indeptest=TRUE)
