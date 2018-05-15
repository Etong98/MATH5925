# SIMULATION: SCENARIO C3

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
R = 10
d = 2
D = d+1

n.train = 300
alpha = 0.95 # quantile level

# True quantile function
quantile_true = function(data, tau, dist, par, theta){
  u1 = do.call(paste(c('p', dist[2]), collapse=''), c(list(x=data[,1]), par[[2]]))
  u2 = do.call(paste(c('p', dist[3]), collapse=''), c(list(x=data[,2]), par[[3]]))
  prob = ((tau^(-theta/(1+2*theta))-1)*(u1^(-theta)+u2^(-theta)-1)+1)^(-1/theta)
  return(do.call(paste(c('q', dist[1]), collapse=''), c(list(p=prob), par[[1]])))
}

# Scenario parameters
delta = 0.86; 
#delta = 4.67
#dist = c('norm', 't', 'norm')
#para = list(list(mean=0, sd=1), list(df=4), list(mean=1, sd=sqrt(4)))
dist = c('st', 'sn', 'st')
para = list(list(xi=0, omega=1, alpha=2, nu=4), list(xi=-2, omega=sqrt(0.5), alpha=3),
                 list(xi=1, omega=sqrt(2), alpha=5, nu=3))


clusterEvalQ(cl, list(library(sn), library(quantreg), library(mboost), library(np), library(copula), library(MASS)))
clusterExport(cl, list('R', 'n.train', 'alpha', 'quantile_true', 'delta', 'dist', 'para'))

start_time = proc.time()
ISE <- foreach(r=1:R) %dopar% {
  # Create copula object
  copula = claytonCopula(param=delta,dim=D)
  mv <- mvdc(copula, margins=dist,paramMargins=para)

  # Genearte training data
  sim.train = data.frame(rMvdc(n.train, mv))
  # Genearte evaulation data 
  n.eval <- n.train/2
  sim.eval <- data.frame(rMvdc(n.eval, mv))
  
  # Linear quanitle regression (LQR)
  LQR <- rq(X1~., tau=alpha, data=sim.train)
  
  # Boosting additive (BAQR)
  it <- 100
  bc <- boost_control(mstop = it, nu=0.25, trace = TRUE, risk = "oob")
  BAQR <- gamboost(X1~., data=sim.train, control=bc, family=QuantReg(tau=alpha))
  
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
  q = quantile_true(data=sim.eval[,2:D], dist=dist, par=para, tau=alpha, theta=delta)

  # Integrated square error
  ISE.LQR = sum((q-q.LQR)^2)/n.eval
  ISE.BAQR = sum((q-q.BAQR)^2)/n.eval
  ISE.DVQR = sum((q-q.DVQR)^2)/n.eval
  output = list(ISE.DVQR, ISE.LQR, ISE.BAQR, DVQR)
  output
}
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

