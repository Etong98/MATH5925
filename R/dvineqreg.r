library(VineCopula)
library(logspline)
library(ks)

# # Testing
# library(MASS)
# rm(list=ls())
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

DVine = function(data){
  d = dim(data)[2]-1; n = dim(data)[1]
  U = matrix(nrow=n, ncol=d)
  
  # PIT with continous kernel smoothing estimator
  for (i in 1:d){
    U[,i] <- kcde(data[,i+1], eval.points=data[,i+1])$estimate
  }
  V <- kcde(data[,1], eval.points = data[, 1])$estimate
  
  # PIT with logspline
  ls = logspline(data[,1], error.action=2)
  # V = plogspline(data[,1], fit=ls)
  # for (i in 1:d){
  #   ls_x = logspline(data[,i], error.action=1)
  #   if (!is.null(ls_x)){
  #     U[,i] = plogspline(data[,i], fit=ls_x)
  #   }
  # }
  
  # Allocate empty vector and matrix
  cop.VU.cll = vector(length=d)
  cop.VU = vector('list', length=d)
  cop.UU = vector('list', length=d)
  U1 = vector(length=n)
  U2 = matrix(nrow=n, ncol=d)
  order = vector(mode='integer', length=d)

  # Initialize global variable
  I <- 1:d
  global.max.cll = -Inf
  cll = 0
  for (step in 1:d){
    if (step == 1){
      loglik_max = -Inf
      for (i in 1:d){
        cop.VU_temp = BiCopSelect(V, U[,i], indeptest=TRUE)
        if (cop.VU_temp$logLik > loglik_max){
          loglik_max = cop.VU_temp$logLik
          l = i
          cop.VU[[step]] = cop.VU_temp
          U2[,1] = U[,i]
        }
      }
      U1 = V
    }else{
      U2_prev = U2; U1_prev = U1
      
      U1 = BiCopHfunc2(U1_prev, U2_prev[,step-1], obj=cop.VU[[step-1]])
      U2_temp = matrix(nrow=n, ncol=d)
      
      loglik_max = -Inf
      cop.UU_temp = vector('list', length=d)
      for (i in I){
        h1 = U[, i]; h2 = U[,l]
        for (k in 2:step){
          cop.UU_temp[[k]] = BiCopSelect(h1, h2, indeptest=TRUE)
          U2_temp[, k] = BiCopHfunc1(h1, h2, cop.UU_temp[[k]])
          h1 = BiCopHfunc2(h1, h2, obj=cop.UU_temp[[k]])
          h2 = U2_prev[, k]
        }
        cop.VU_temp = BiCopSelect(U1, h1, indeptest=TRUE)
        if (cop.VU_temp$logLik > loglik_max){
          loglik_max = cop.VU_temp$logLik
          U2 = U2_temp
          cop.UU[[step]] = cop.UU_temp
          cop.VU[[step]] = cop.VU_temp
          l = i
        }
      }
    }
    
    cll = cll+loglik_max
    if (cll <= global.max.cll){
      order = order[1:step-1]
      return(list(ord = order, copVU = cop.VU, copUU = cop.UU, logSpline=ls))
    }
    global.max.cll = cll
    I = I[! I %in% l]
    order[step] = l
  }
  return(list(ord = order, copVU = cop.VU, copUU = cop.UU, logSpline=ls))
}


DVQR.quantile = function(obj, newdata, tau){
  order = obj$ord; cop.VU = obj$copVU; cop.UU = obj$copUU
  
  d = length(order); n = dim(newdata)[1]
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
      for (k in 2:step){
        UU[, k] = BiCopHfunc1(h1, h2, obj=cop.UU[[step]][[k]])
        h1 = BiCopHfunc2(h1, h2, obj=cop.UU[[step]][[k]])
        h2 = UU_prev[, k]
      }
      U2[,step] = h1
    }
  }
  
  # Calculate conditional quantile
  cqf = rep(alpha,n)
  for (step in d:1){
    cqf = BiCopHinv2(cqf, U2[,step], obj=cop.VU[[step]])
  }
  quantile = qlogspline(cqf, fit=obj$logSpline)
  return(quantile)
}
