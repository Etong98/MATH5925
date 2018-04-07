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
data <- data.frame(mvrnorm(n, mean, cov_mat))
colnames(data) <-  c('Y','X1','X2','X3')


  d <- dim(data)[2]-1
  n <- dim(data)[1]
  V <- kcde(data[,1], eval.points = data[, 1])$estimate
  U <- matrix(nrow=n, ncol=d)
  for (i in 1:d){
    U[,i] <- kcde(data[,i+1], eval.points=data[,i+1])$estimate
  }
  cop.VU.cll <- vector(length=d)
  I <- 1:d
  global.max.cll <- -Inf
  for (i in I){
    cop.VU.cll[i] <- BiCopSelect(V, U[,i], indeptest=TRUE)$logLik
  }
  l <- which.max(cop.VU.cll)
  global.max.cll <- max(cop.VU.cll)
  I <- I[! I %in% l]
  
  return(I)

