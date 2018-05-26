# # Testing
library(MASS)

rm(list=ls())
# n <- 500
# d <- 3
# mean <- rep(0, 4)
# cov_mat <- matrix(c(1, 0.4, 0.8, 0,
#                     0.4, 1, 0.32, 0,
#                     0.8, 0.32, 1, 0,
#                     0, 0, 0, 1),
#                   nrow=4, ncol=4)
# data <- data.frame(mvrnorm(n, mean, cov_mat))
# colnames(data) <-  c('Y','X1','X2','X3')

set.seed(100)

d=4
n=300

sigma = 0.1
fY = function(x){
  y = sqrt(abs(2*x[,1]-x[,2]+0.5))+(-0.5*x[,3]+1)*(0.1*x[,4]^3)
  return(y)
}
S_X = matrix(NA, 4, 4)
for (i in 1:d){
  for (j in 1:d){
    S_X[i,j] = 0.5^abs(i-j)
  }
}
mu_X = rep(0, 4)

X.train = mvrnorm(n, mu_X, S_X)
Y.train = fY(X.train)+sigma*rnorm(n)
data = data.frame(cbind(Y.train, X.train))

source('R/dvineqreg2.r', local=TRUE)
DVQR2 = DVine(data)

source('R/dvineqreg1.r', local=TRUE)
DVQR1 = Dvine1(data)
#q.DVQR = DVQR.quantile1(X.train, DVQR2, 0.5)


#DVQR1 = CDVineCopSelect(cbind(V,U[,c(4,1,2,3)]), type=2, indeptest=TRUE)
