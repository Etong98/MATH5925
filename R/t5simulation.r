# SIMULATION; SCENARIO t5

rm(list = ls()) # clear workspace variable

library(parallel)
library(doSNOW)
library(foreach)

num_core = detectCores()
cl = makeCluster(num_core)
registerDoSNOW(cl)

# Simulation parameter
R = 50
d = 4; D = d+1
n.train = 300
alpha = 0.95 # quantile level

# True quantile function
quantile_true = function(data, tau, df, corr, sd, mu){
  n = dim(data)[1]
  d = dim(data)[2]
  D = d+1
  df_c = df+d
  
  S = diag(sd)%*%corr%*%diag(sd)
  S_XX = as.matrix(S[2:D,2:D])
  S_YX = as.matrix(S[1,2:D])
  x_bar = as.matrix(sweep(data,2,mu[2:D]))
  
  mu_c = rep(NA, n); var_c = rep(NA, n)
  for (i in 1:n){
    mu_c[i] = mu[1]+t(S_YX)%*%solve(S_XX)%*%x_bar[i,]
    var_c[i] = (df+t(x_bar[i,])%*%solve(S_XX)%*%x_bar[i,])*
                (S[1,1]-t(S_YX)%*%solve(S_XX)%*%S_YX)/df_c
  }
  return(qt.scaled(tau, df=df_c, mean=mu_c, sd=sqrt(var_c)))  
}


# Scenario paramters
v = 3 # df
#R_vec = c(0.6, 0.5, 0.5, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
R_vec = c(0.27, 0.74, 0.72, 0.41, 0.28, 0.29, 0.27, 0.74, 0.42, 0.40)
#dist = c('norm', 't', 'norm', 't', 'norm')
#para = list(list(mean=0), list(df=4), list(mean=1,sd=2), list(df=4), list(mean=1,sd=2))
dist = c('st', 'sn', 'st', 'sn', 'st')
para = list(list(alpha=2,nu=4), list(xi=-2,omega=sqrt(0.5),alpha=3), list(xi=1,omega=sqrt(2),alpha=5,nu=3), 
            list(xi=-2,omega=sqrt(0.5),alpha=3), list(xi=1,omega=sqrt(2),alpha=5,nu=3))
mu = c(0, 0, 1, 0, 1)
sd = c(1, sqrt(2), 2, sqrt(2),2)

clusterEvalQ(cl, list(library(sn), library(quantreg), library(mboost), library(np), library(copula), library(MASS), library(metRology)))
clusterExport(cl, list('R', 'n.train', 'alpha', 'd', 'D', 'quantile_true', 'dist', 'para', 'mu', 'sd', 'v'))

start_time = proc.time()
ISE <- foreach(r=1:R) %dopar% {
  copula = tCopula(R_vec, dim=D, dispstr='un', df=v)
  mv <- mvdc(copula, margins=dist, paramMargins=para)

  # Genearte training data
  sim.train = data.frame(rMvdc(n.train, mv))

  # Genearte evaulation data 
  n.eval = n.train/2
  sim.eval = data.frame(rMvdc(n.eval, mv))
  
  # Linear quanitle regression (LQR)
  LQR <- rq(X1~., tau=alpha, data=sim.train)
  
  # Boosting additive (BAQR)
  it <- 100
  bc <- boost_control(mstop = it, nu=0.25, trace = TRUE, risk = "oob")
  BAQR <- gamboost(X1~., data=sim.train, control=bc, family=QuantReg(tau=alpha))
  
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
  q = quantile_true(data=sim.eval[,2:D], tau=alpha, df=v, corr=getSigma(copula), sd=sd, mu=mu)
  
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
  MISE.DVQR =  MISE.DVQR+ISE[[i]][[2]]/R
  MISE.LQR =  MISE.LQR+ISE[[i]][[3]]/R
  MISE.BAQR = MISE.BAQR+ISE[[i]][[4]]/R
}
MISE.SPQR
MISE.DVQR
MISE.LQR
MISE.BAQR

# Stop parallel
stopCluster(cl)
proc.time()-start_time

