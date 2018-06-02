require(openxlsx)
require(quantreg)
require(mboost)


ASX30 = read.xlsx('data/ASX30.xlsx',detectDates = TRUE)

N = ncol(ASX30)-1
n = nrow(ASX30)

M = matrix(nrow=n,ncol=N)
for (j in 1:N){
  for (i in 1:n){
    M[i,j] = as.numeric(ASX30[i,j+1])
    if (is.na(M[i,j])){
      M[i,j] = M[i-1,j]
    }
  }
}
ASX30[,1:N+1] = M

price_company = ASX30[1:N+1]

R_company = matrix(NA, nrow=n,ncol=N)
for (j in 1:N){
  for (i in 2:n){
    R_company[i,j] = 100*(log(price_company[i,j])-log(price_company[i-1,j]))
  }
}
R_company = R_company[2:n,]

n = nrow(R_company)
n.train = n*3/4

alpha = 0.95
selection = c(8,9,10)
d = length(selection)-1
data.train = data.frame(R_company[1:n.train, selection])
data.eval = data.frame(R_company[(n.train+1):n, selection])

#
source('R/SPqreg.r')
RV = SPvine(data.train)
q.SPQR = SPQR.quantile(newdata=data.eval[,1:d+1], Y=data.train[, 1], obj=RV, tau=alpha)
#
source('R/dvineqreg.r')
DVQR = Dvine(data.train)
q.DVQR = DVQR.quantile(newdata=data.eval[,1:d+1], obj=DVQR, tau=alpha)
#
LQR = rq(X1~., tau=alpha, data=data.train)
q.LQR = predict.rq(LQR, newdata=data.eval[,1:d+1])
# 
it = 100
bc = boost_control(mstop = it, nu=0.25, trace = TRUE, risk = "oob")
BAQR = gamboost(X1~., data=data.train, control=bc, family=QuantReg(tau=alpha))
q.BAQR = predict(BAQR, newdata=data.eval[,1:d+1])

L = function(q, tau, y){
  n = length(y)
  rho = rep(NA,n)
  for (i in 1:n){
    if ((y[i]-q) < 0){rho[i] = (y[i]-q)*(tau-1)}else{rho[i] = (y[i]-q)*tau}  
  }
  return(sum(rho))
}

L.SPQR = L(q.SPQR, alpha, data.eval[,1])/(n-n.train)
L.DVQR = L(q.DVQR, alpha, data.eval[,1])/(n-n.train)
L.LQR = L(q.LQR, alpha, data.eval[,1])/(n-n.train)
L.BAQR = L(q.BAQR, alpha, data.eval[,1])/(n-n.train)

L.SPQR
L.DVQR
L.LQR
L.BAQR
