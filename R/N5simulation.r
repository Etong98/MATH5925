# SIMULATION: SCENARIO N5

rm(list = ls()) # clear workspace variable

library(parallel)
library(doSNOW)
library(foreach)

num_core = detectCores()
cl = makeCluster(num_core)
registerDoSNOW(cl)

# Simulation parameters
R = 10
d = 4
D = d+1

n.train = 300
alpha = 0.95 # quantile level

# True quantile function
quantile_true = function(newdata, tau, func, sd){
  return(qnorm(tau, mean=func(newdata), sd=sd))
}

# Scenario paramters
sigma = 1
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

clusterEvalQ(cl, list(library(sn), library(quantreg), library(mboost), library(np), library(MASS)))
clusterExport(cl, list('R', 'n.train', 'alpha', 'quantile_true', 'S_X', 'mu_X', 'd', 'D','fY', 'sigma'))

start_time = proc.time()
ISE <- foreach(r=1:R) %dopar% {
  # Genearte training data
  X.train = mvrnorm(n.train, mu_X, S_X)
  Y.train = fY(X.train)+sigma*rnorm(n.train)
  sim.train = data.frame(cbind(Y.train, X.train))
  
  # Genearte evaulation data 
  n.eval <- n.train/2
  X.eval = mvrnorm(n.eval, mu_X, S_X)
  Y.eval = fY(X.eval)+sigma*rnorm(n.eval)
  sim.eval = data.frame(cbind(Y.eval, X.eval))
  
  # Linear quanitle regression (LQR)
  LQR <- rq(Y.train~., tau=alpha, data=sim.train)
  
  # Boosting additive (BAQR)
  it <- 100
  bc <- boost_control(mstop = it, nu=0.25, trace = TRUE, risk = "oob")
  BAQR <- gamboost(Y.train~., data=sim.train, control=bc, family=QuantReg(tau=alpha))
  
  # Semiparametric quantile regression (SPQR)
  source('R/SPqreg.R')
  SPQR = SPvine(sim.train)
  
  # D-vine quantile regression (DVQR)
  source("R/dvineqreg.R")
  DVQR = Dvine(sim.train)
  
  # Predict quanile/calculate true quantile
  q.LQR = predict.rq(LQR, newdata=sim.eval[,2:D])
  q.BAQR = predict(BAQR, newdata=sim.eval[,2:D])
  #q.DVQR = DVQR.quantile(DVQR, newdata=sim.eval[,2:D], tau=alpha)
  q.DVQR = DVQR.quantile(obj=DVQR, newdata=sim.eval[,2:D], tau=alpha)
  q.SPQR = SPQR.quantile(obj=SPQR, newdata=sim.eval[,2:D], Y=sim.train[,1], tau=alpha)
  q = quantile_true(newdata=sim.eval[,2:D], tau=alpha, func=fY, sd=sigma)
  
  # Integrated square error
  ISE.LQR = sum((q-q.LQR)^2)/n.eval
  ISE.BAQR = sum((q-q.BAQR)^2)/n.eval
  ISE.DVQR = sum((q-q.DVQR)^2)/n.eval
  ISE.SPQR = sum((q-q.SPQR)^2)/n.eval
  output = list(ISE.SPQR, ISE.DVQR, ISE.LQR, ISE.BAQR)
  output
}

# Mean Integrated Square Error
MISE.SPQR = 0; MISE.DVQR = 0; MISE.LQR = 0; MISE.BAQR = 0
for (i in 1:R){
  MISE.SPQR = MISE.SPQR+ISE[[i]][[1]]/R
  MISE.DVQR = MISE.DVQR+ISE[[i]][[2]]/R
  MISE.LQR = MISE.LQR+ISE[[i]][[3]]/R
  MISE.BAQR = MISE.BAQR+ISE[[i]][[4]]/R
}
MISE.SPQR
MISE.DVQR
MISE.LQR
MISE.BAQR

# Stop parallel
stopCluster(cl)
proc.time()-start_time

