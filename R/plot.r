rm(list = ls())

library(sn)
library(quantreg)
library(mboost)
library(np)
library(VineCopula)
library(copula)
library(ks)
library(graphics)

source('DVine_quantreg/simulation.r')

num_train = 300

clayton_quantile <- function(v, w, tau, theta){
  prob <- ((tau^(-theta/(1+2*theta))-1)*(v^(-theta)+w^(-theta)-1)+1)^(-1/theta)
  return(qnorm(prob))
}

clayton.cop <- claytonCopula(param=0.86,dim=3)
mv.C3 <- mvdc(clayton.cop, margins=c('norm', 't', 'norm'), 
              paramMargins=list(list(mean=0, sd=1), list(df=1,ncp=0), list(mean=1, sd=4)))

train <- data.frame(rMvdc(num_train, mv.C3))
colnames(train) <- c('Y', 'X1', 'X2')

x1 = seq(-3, 3, length.out = 50)
x2 = seq(-3, 5, length.out = 50)
u1 = pt(x1, df=1, ncp=0)
u2 = pnorm(x2, mean=1, sd=4)

plot.new()
fun1 <- function(v, w){
  tau=0.5
  theta=0.86
  return(clayton_quantile(v, w, tau, theta))
}
z = outer(u1, u2, fun1)
persp3Drgl(x1, x2, z, theta=330, phi=30, ticktype='detailed', main='median', col='red')
fun2 <- function(v, w){
  tau=0.95
  theta=0.86
  return(clayton_quantile(v, w, tau, theta))
}
z = outer(u1, u2, fun2)
persp3Drgl(x1, x2, z, theta=330, phi=30, ticktype='detailed', main='95%', add=TRUE)
