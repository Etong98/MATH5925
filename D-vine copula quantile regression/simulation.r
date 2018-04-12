rm(list = ls())

library(sn)
library(quantreg)
library(mboost)
library(np)
library(VineCopula)
library(copula)
library(ks)

R <- 11
n.train <- 300
alpha <- 0.5

MISE.LQR <- 0
MISE.BAQR <- 0

clayton.cop <- claytonCopula(param=0.86,dim=3)
mv.C3 <- mvdc(clayton.cop, margins=c('norm', 't', 'norm'), 
              paramMargins=list(list(mean=0, sd=1), list(df=1,ncp=0), list(mean=1, sd=4)))

source("D-Vine copula quantile regression\\DVine_regression.R")


for (r in 1:R){
  sim.train <- data.frame(rMvdc(n.train, mv.C3))
  colnames(sim.train) <- c('Y', 'X1', 'X2')
  DVine_regression(sim.train)
  
  # LQR <- rq(Y~X1+X2, data=sim.train)
  # BAQR <- gamboost(Y~X1+X2, data=sim.train)

  n.eval <- n.train/2
  sim.eval <- data.frame(rMvdc(n.eval, mv.C3))
  colnames(sim.eval) <- c('Y', 'X1', 'X2')
  
  # q.LQR <- predict.rq(LQR, newdata=sim.eval[2:3])
  # q.BAQR <- predict(BAQR, newdata=sim.eval[2:3])
  # 
  # MISE.LQR <- MISE.LQR+1/(n.eval*R)*sum((sim.eval$Y-q.LQR)^2)
  # MISE.BAQR <- MISE.BAQR+1/(n.eval*R)*sum((sim.eval$Y-q.BAQR)^2)
  
  
}


