source("SP_quantreg.R")
bw = npcdistbw(ydat=data.train$Y, xdat=data.train[,2:3])
q.NPQR = npqreg(bws=bw, newdata=data.eval[,2:3])
