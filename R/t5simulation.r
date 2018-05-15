rm(list = ls()) # clear workspace variable

library(sn)
library(quantreg)
library(mboost)
library(np)
library(copula)
library(MASS)

library(parallel)
library(doSNOW)
library(foreach)

num_core = detectCores()
cl = makeCluster(num_core)
registerDoSNOW(cl)

# Simulation parameter
R = 100
d = 4
D = 5

n.train = 300
alpha = 0.95 # quantile level

# True quantile function
quantile_C3 = function(x1, x2, tau){
  theta = 0.86
  #u1 = pt(x1, df=4)
  #u2 = pnorm(x2, mean=1, sd=sqrt(4))
  u1 = psn(x1, xi=-2, omega=sqrt(0.5), alpha=3); u2 = pst(x2, xi=1, omega=sqrt(2), alpha=5, nu=3)
  prob = ((tau^(-theta/(1+2*theta))-1)*(u1^(-theta)+u2^(-theta)-1)+1)^(-1/theta)
  return(qnorm(prob))
}

# quantile2 = function(X, tau, R, dist, para){
#   d = dim[X][2]
#   n = dim[X][1]
#   sigma_yx = R[1, 2:(d+1)]
#   sigma_xx = R[2:(d+1),2:(d+1)]
#   
# }

quantile_N5 = function(newdata, tau, func, sd){
  return(qnorm(tau, mean=func(newdata), sd=sd))
}

### Scenario 1: C3 ----
# delta
delta = 0.86; 
#delta = 4.67
# Margin
#dist = c('norm', 't', 'norm')
#para = list(list(mean=0, sd=1), list(df=4), list(mean=1, sd=sqrt(4)))
dist = c('st', 'sn', 'st')
para = list(list(xi=0, omega=1, alpha=2, nu=4), list(xi=-2, omega=sqrt(0.5), alpha=3),
                 list(xi=1, omega=sqrt(2), alpha=5, nu=3))

### Scenario 2: t5 ----
# # R
#R = matrix(c(1, 0.6, 0.5, 0.5, 0.4, 0.6, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 
#              0.5, 0.5, 0.5, 1, 0.5, 0.4, 0.5, 0.5, 0.5, 1), nrow=5, ncol=5)
## Margin
# dist = c('norm', 't', 'norm', 't', 'norm')
#para = list(list(mean=0, sd=1), list(df=4), list(mean=1,sd=sqrt(4)), list(df=4), list(mean=1, sd=sqrt(4)))
# dist = c('st', 'sn', 'st', 'sn', 'st')
# para = list(list(xi))

### Scenario 3: M5 ----
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


# Initalize Integrated Sqauare Error vec
# ISE.LQR = rep(0,R); ISE.BAQR = rep(0,R); ISE.NPQR = rep(0,R); ISE.DVQR = rep(0, R)

clusterEvalQ(cl, list(library(sn), library(quantreg), library(mboost), library(np), library(copula), library(MASS)))
clusterExport(cl, list('R', 'n.train', 'alpha', 'quantile_N5', 'delta', 'dist', 'para', 'S_X', 'mu_X', 'd', 'D','fY', 'sigma'))

start_time = proc.time()
ISE <- foreach(r=1:R) %dopar% {
  #copula = claytonCopula(param=delta,dim=3)
  #mv <- mvdc(copula, margins=dist,paramMargins=para)

  # Genearte training data
  #sim.train = data.frame(rMvdc(n.train, mv))
  X.train = mvrnorm(n.train, mu_X, S_X)
  Y.train = fY(X.train)+sigma*rnorm(n.train)
  sim.train = data.frame(cbind(Y.train, X.train))
  # Genearte evaulation data 
  n.eval <- n.train/2
  #sim.eval <- data.frame(rMvdc(n.eval, mv))
  X.eval = mvrnorm(n.eval, mu_X, S_X)
  Y.eval = fY(X.eval)+sigma*rnorm(n.train)
  sim.eval = data.frame(cbind(Y.eval, X.eval))
  
  # Linear quanitle regression (LQR)
  LQR <- rq(Y.train~., tau=alpha, data=sim.train)
  
  # Boosting additive (BAQR)
  it <- 100
  bc <- boost_control(mstop = it, nu=0.25, trace = TRUE, risk = "oob")
  BAQR <- gamboost(Y.train~., data=sim.train, control=bc, family=QuantReg(tau=alpha))
  
  # D-vine quantile regression (DVQR)
  source("R/dvineqreg.R")
  DVQR = DVine(sim.train)
  
  # Non parametric quantile regression
  # bw = npcdistbw(ydat=sim.train$Y, xdat=sim.train[,2:3])
  
  # Predict quanile/calculate true quantile
  q.LQR = predict.rq(LQR, newdata=sim.eval[,2:D])
  q.BAQR = predict(BAQR, newdata=sim.eval[,2:D])
  #q.NPQR = npqreg(bws=bw, newdata=sim.eval[,2:3])$quantile
  q.DVQR = DVQR.quantile(DVQR, newdata=sim.eval[,2:D], tau=alpha)
  #q = quantile(x1=sim.eval$X2, x2=sim.eval$X3, tau=alpha)
  q = quantile_N5(newdata=sim.eval[,2:D], tau=alpha, func=fY, sd=sigma)
  
  # Integrated square error
  ISE.LQR = sum((q-q.LQR)^2)/n.eval
  ISE.BAQR = sum((q-q.BAQR)^2)/n.eval
  ISE.DVQR = sum((q-q.DVQR)^2)/n.eval
  # ISE.NPQR[r] = sum((q-q.NPQR)^2/n.eval)
  output = list(ISE.DVQR, ISE.LQR, ISE.BAQR, DVQR)
  output
}
# MISE.LQR = sum(ISE.LQR)/R
# MISE.LQR
# MISE.BAQR = sum(ISE.BAQR)/R
# MISE.BAQR
# # MISE.NPQR = sum(ISE.NPQR)/R
# # MISE.NPQR
# MISE.DVQR = sum(ISE.DVQR)/R
MISE.DVQR = 0; MISE.LQR = 0; MISE.BAQR = 0
for (i in 1:R){
  MISE.DVQR =  MISE.DVQR+ISE[[i]][[1]]/R
  MISE.LQR =  MISE.LQR+ISE[[i]][[2]]/R
  MISE.BAQR = MISE.BAQR+ISE[[i]][[3]]/R
}
MISE.DVQR
MISE.LQR
MISE.BAQR

# Stop parallel
stopCluster(cl)
proc.time()-start_time

