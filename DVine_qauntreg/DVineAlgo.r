library(MASS)
library(VineCopula)
library(copula)
library(ks)

n <- 500
d <- 3
mean <- rep(0, 4)
cov_mat <- matrix(c(1, 0.4, 0.8, 0,
                    0.4, 1, 0.32, 0,
                    0.8, 0.32, 1, 0, 
                    0, 0, 0, 1),
                  nrow=4, ncol=4)
data <- data.frame(mvrnorm(n, mean, cov_mat))
colnames(data) <-  c('Y','X1','X2','X3')

DVine <- function(data){
  d <- dim(data)[2]-1
  n <- dim(data)[1]
  V <- kcde(data[,1], eval.points = data[, 1])$estimate
  U <- matrix(nrow=n, ncol=d)
  for (i in 1:d){
    U[,i] <- kcde(data[,i+1], eval.points=data[,i+1])$estimate
  }
  cop.VU.cll <- vector(length=d)
  cop.VU <- vector('list', length=d)
  cop.UU <- vector('list', length=d)
  
  I <- 1:d
  global.max.cll <- -Inf
  order <- rep(0, d)
  k=1
  for (i in I){
    if (k == 1){
      for (i in 1:d){
        cop.VU[[i]] <- BiCopSelect(V, U[,i], indeptest=TRUE)
        cop.VU.cll[i] <- cop.VU[[i]]$logLik
      }      
    }else{
      V <- BiCopHfunc2(V, U[,l], obj=cop.VU[[l]])
      I <- I[! I %in% l]
      for (j in I){
        cop.UU[[j]] <- BiCopSelect(U[,j], U[,l], indeptest=TRUE)
        U[,j] <- BiCopHfunc2(U[,j], U[,l], obj=cop.UU[[j]])
        cop.VU[[j]] <- BiCopSelect(V, U[,j], indeptest=TRUE)
        cop.VU.cll[j] <- cop.VU[[j]]$logLik
      }  
    }
    if (max(cop.VU.cll) <= 0){
      return(order)
    }
    l <- which.max(cop.VU.cll)
    cop.VU.cll[l] <- -Inf
    order[k] <- l
    k <- k+1
  }
  return(order)
}
