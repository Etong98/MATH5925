library(sn)
library(quantreg)
library(mboost)
library(np)
library(VineCopula)

rm(list = ls())

R <- 100
n <- 300

scenario <- c('C3', 't5', 'M5')

Y <- rnorm(n)
X1 <- rt(n, df=4)
X2 <- rnorm(n, mean=1, sd=4)
sim <- data.frame(cbind(Y, X1, X2))

LQR <- rq(Y~X1+X2, data=sim)
BAQR <- gamboost(Y~X1+X2, data=sim)
NPQR <- npcdistbw(Y~X1+X2, data=sim)

rvine <- RVineCopSelect(sim)