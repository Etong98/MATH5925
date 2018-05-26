library(VineCopula)
library(CDVine)
library(logspline)
library(ks)


Dvine = function(data){
  d = ncol(data)-1
  n = nrow(data)
  
  # PIT with continous kernel smoothing estimator
  U = matrix(nrow=n, ncol=d)
  for (i in 1:d){
    U[,i] = kcde(data[,i+1], eval.points=data[,i+1])$estimate
  }
  #V = kcde(data[,1], eval.points = data[, 1])$estimate
  
  # PIT with logspline
  ls = logspline(data[,1], error.action=2)
  V = plogspline(data[,1], fit=ls)
  
  order = 0[-1]
  global_AIC = Inf
  M = 1:d
  for (step in 1:d){
    min_AIC = Inf
    for (i in M){
      if (step == 1){
        dv = BiCopSelect(V, U[,i], indeptest=TRUE)
      }else{
        dv = CDVineCopSelect(cbind(V,U[,append(order,i)]), type=2, indeptest=TRUE)
      }
      AIC = CDVineAIC(cbind(V,U[,append(order,i)]), dv$family, dv$par, dv$par2, type=2)$AIC
      if (AIC < min_AIC){
        min_AIC = AIC
        m = i
        DV = dv
      }
    }
    if (min_AIC == 0){
      return(list(DV, order, ls))
    }
    global_AIC = min_AIC
    M = M[!M %in% m]
    order = append(order,m)
  }
  return(list(DV, order, ls))
}


DVQR.quantile = function(newdata, obj, tau){
  order = obj[[2]]; fam = obj[[1]]$family; par = obj[[1]]$par; par2 = obj[[1]]$par2
  
  d = length(order); n = nrow(newdata)
  U = matrix(nrow=n, ncol=d); U2 = matrix(nrow=n, ncol=d)
  
  # Probability integral transform
  j=1
  for (i in order){
    U[,j] = kcde(newdata[,i], eval.points=newdata[,i])$estimate
    j = j+1
  }
  
  # Calculate conditioned distribution
  UU = matrix(nrow=n, ncol=d); U2 = matrix(nrow=n, ncol=d)
  for (step in 1:d){
    if (step == 1){
      UU[,1] = U[,1]
      U2[,step] = UU[,1]
    }else{
      UU_prev = UU
      h1 = U[, step]; h2 = U[,step-1]; UU[,1] = U[,step-1]
      inc = d-1; index = step
      for (k in 2:step){
        if (k==2){
          UU[, k] = BiCopHfunc1(h1, h2, family=fam[step], par=par[step], par2=par2[step])
          h1 = BiCopHfunc2(h1, h2, family=fam[step], par=par[step], par2=par2[step])
          h2 = UU_prev[,k]
        }else{
          index = index+inc
          UU[, k] = BiCopHfunc1(h1, h2, family = fam[index], par=par[index], par2=par2[index])
          h1 = BiCopHfunc2(h1, h2, family = fam[index], par=par[index], par2=par2[index])
          h2 = UU_prev[,k]
          inc = inc-1
        }
      }
      U2[,step] = h1
    }
  }
  
  # Calculate conditional quantile
  cqf = rep(tau,n)
  index = (d+1)*d/2; inc = 2
  for (step in d:1){
    if (step == d){
      cqf = BiCopHinv2(cqf, U2[,step], family=fam[index], par=par[index], par2=par2[index])
    }else{
      index = index-inc
      cqf = BiCopHinv2(cqf, U2[,step], family=fam[index], par=par[index], par2=par2[index])
      inc = inc+1
    }
  }
  quantile = qlogspline(cqf, fit=obj[[3]])
  return(quantile)
}