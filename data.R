
library(matrixcalc)
library(mvtnorm)
library(Matrix)
library(MCMCpack)
library(R2Cuba)
library(abind)
library(ggplot2)


Vcros=as.matrix(read.csv(file="Vcros.csv",sep=",",header=TRUE))
Vcr=as.matrix(read.csv(file="Vcr.csv",sep=",",header=TRUE))


cros_cov=as.matrix(read.csv(file="CROS_COV.csv",sep=",",header=TRUE,stringsAsFactors=FALSE))

indxct=c(5,6,8,16:18,20:22,50,51,58,
         78:84,113,115,121,
         158:161,171, 172, 181,185,187:192,194,195,
         196:199)
indxos = setdiff(1:205, indxct)
indx=c(indxct, indxos)
(p=length(indx))
(k1=length(indxct))
cros_cov = cros_cov[indx,]
lamcthat=as.numeric(cros_cov[1:k1,7])
nucthat=log(lamcthat)

psi_os=as.numeric(cros_cov[,8])
lam_os=-log(1-psi_os)

ph=as.numeric(cros_cov[,3])
ap=as.numeric(cros_cov[,4])
lp=as.numeric(cros_cov[,5])
hi=as.numeric(cros_cov[,6])

X1=cbind(lam_os)
X2=cbind(lam_os, hi)
X3=cbind(lam_os, ap)
X4=cbind(lam_os, lp)
X5=cbind(lam_os, ph)
X6=cbind(lam_os, hi, ap)
X7=cbind(lam_os, hi, lp)
X8=cbind(lam_os, hi, ph)
X9=cbind(lam_os, ap, ph)
X10=cbind(lam_os, lp, ph)
X11=cbind(lam_os, hi, ap, ph)
X12=cbind(lam_os, hi, lp, ph)

XX=list(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12)

Xall=cbind(lam_os, ap, lp, hi, ph)

