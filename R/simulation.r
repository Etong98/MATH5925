rm(list = ls()) # clear workspace variable

library(sn)
library(quantreg)
library(mboost)
library(np)
library(copula)

library(parallel)
library(doSNOW)
library(foreach)

num_core = detectCores()
cl = makeCluster(num_core)
registerDoSNOW(cl)

# Simulation parameter
R = 10
n.train = 300
alpha = 0.5 # quantile level

# True quantile function
quantile = function(x1, x2, tau){
  theta = 0.86
  u1 = pt(x1, df=1, ncp=0)
  u2 = pnorm(x2, mean=1, sd=4)
  # u1 = psn(x1, xi=-2, omega=0.5, alpha=3); u2 = pst(x2, xi=1, omega=2, alpha=5, nu=3)
  prob = ((tau^(-theta/(1+2*theta))-1)*(u1^(-theta)+u2^(-theta)-1)+1)^(-1/theta)
  return(qnorm(prob))
}

### Scenario 1: C3 ----
# delta
delta = 0.86; 
#delta = 4.67
# Margin
dist = c('norm', 't', 'norm')
para = list(list(mean=0, sd=1), list(df=1, ncp=0), list(mean=1, sd=4))
#dist = c('st', 'sn', 'st')
#para = list(list(xi=0, omega=1, alpha=2, nu=4), list(xi=-2, omega=0.5, alpha=3),
#                 list(xi=1, omega=2, alpha=5, nu=3))

### Scenario 2: t5 ----
# # R
# R = matrix(c(1, 0.6, 0.5, 0.5, 0.4, 0.6, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 
#              0.5, 0.5, 0.5, 1, 0.5, 0.4, 0.5, 0.5, 0.5, 1), nrow=5, ncol=5)
# # Margin
# dist = c('norm', 't', 'norm')
# para = list(list(mean=0, sd=1), list(df=1, ncp=0), list(mean=1,sd=4))
# dist = c('norm', 't', 'norm')
# para = list(list(mean=0, sd=1), list(df=1, ncp=0), list(mean=1,sd=4))

### Scenario 3: M5 ----

# Initalize Integrated Sqauare Error vec
# ISE.LQR = rep(0,R); ISE.BAQR = rep(0,R); ISE.NPQR = rep(0,R); ISE.DVQR = rep(0, R)

clusterEvalQ(cl, list(library(sn), library(quantreg), library(mboost), library(np), library(copula)))
clusterExport(cl, list('R', 'n.train', 'alpha', 'quantile', 'delta', 'dist', 'para'))

start_time = proc.time()
ISE <- foreach(r=1:R) %dopar% {
  copula = claytonCopula(param=delta,dim=3)
  #copula = tCopula
  mv <- mvdc(copula, margins=dist,
                paramMargins=para)

  # Genearte training data
  sim.train <- data.frame(rMvdc(n.train, mv))
  colnames(sim.train) <- c('Y', 'X1', 'X2')
  V <- pnorm(sim.train$Y)
  
  # Linear quanitle regression (LQR)
  LQR <- rq(Y~X1+X2, tau=alpha, data=sim.train)
  
  # Boosting additive (BAQR)
  it <- 100
  bc <- boost_control(mstop = it, nu=0.25, trace = TRUE, risk = "oob")
  BAQR <- gamboost(Y~X1+X2, data=sim.train, control=bc, family=QuantReg(tau=alpha))
  
  # D-vine quantile regression (DVQR)
  source("R/dvineqreg.R")
  DVQR = DVine(sim.train)
  
  # Non parametric quantile regression
  # bw = npcdistbw(ydat=sim.train$Y, xdat=sim.train[,2:3])
  
  # Genearte evaulation data 
  n.eval <- n.train/2
  sim.eval <- data.frame(rMvdc(n.eval, mv))
  colnames(sim.eval) <- c('Y', 'X1', 'X2')
  
  # Predict quanile/calculate true quantile
  q.LQR = predict.rq(LQR, newdata=sim.eval[,2:3])
  q.BAQR = predict(BAQR, newdata=sim.eval[,2:3])
  #q.NPQR = npqreg(bws=bw, newdata=sim.eval[,2:3])$quantile
  q.DVQR = DVQR.quantile(DVQR, newdata=sim.eval[,2:3], tau=alpha)
  q = quantile(x1=sim.eval$X1, x2=sim.eval$X2, tau=alpha)
  
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

