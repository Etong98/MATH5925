library(VineCopula)
library(ks)

SPvine = function(data){
  n = nrow(data); d = ncol(data)-1
  
  U = matrix(nrow=n, ncol=d)
  for (i in 1:d){
    U[,i] = kcde(data[,i+1], eval.points=data[,i+1])$estimate
  }
  V = kcde(data[,1], eval.points=data[,1])$estimate
  
  RV = RVineStructureSelect(cbind(V,U), indeptest=TRUE)
  return(RV)
}

SPQR.quantile = function(newdata, obj, Y, tau){
  n.train = length(Y); d=ncol(newdata); n = nrow(newdata)
  U_mat = matrix(nrow=n, ncol=d)
  for (i in 1:d){
    U_mat[,i] = kcde(newdata[,i], eval.points=newdata[,i])$estimate
  }
  V = kcde(Y, eval.points=Y)$estimate
  
  obj_func = function(q, tau, y, pdf){
    n = length(y)
    rho = rep(NA,n)
    for (i in 1:n){
      if ((y[i]-q) < 0){rho[i] = (y[i]-q)*(tau-1)}else{rho[i] = (y[i]-q)*tau}  
    }
    return(sum(rho*pdf))
  }
  
  est_quantile = rep(NA, n)
  for (i in 1:n){
    U = matrix(U_mat[i,], nrow=n.train, ncol=d, byrow=TRUE)
    cop_pdf = RVinePDF(cbind(V,U),RVM=obj)
    est_quantile[i] = optimize(obj_func, interval=c(min(Y), max(Y)), tau=tau, y=Y, pdf=cop_pdf)$minimum  
  }
  return(est_quantile)
}