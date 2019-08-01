

ftauetaus.pred=function(taue, taus, U011)
{
  U11=(1/taus)*U011
  mat1=chol2inv(chol((1/taue)*diag(k1)+U11))
  
  Sig.beta=(  t(XT)%*%mat1%*%XT  )
  Sig.beta.inv= chol2inv(chol(Sig.beta))
  Sig.beta.inv=as.Symmetric.matrix(Sig.beta.inv)
  detdet = det(Sig.beta.inv)
  beta0hat=Sig.beta.inv%*%t(XT)%*%mat1%*%nucthat
  hhh=t(nucthat)%*%mat1%*%nucthat-t(beta0hat)%*%Sig.beta%*%beta0hat
  
  fac1=dgamma(taue, shape=a.e, rate=b.e)*dgamma(taus, shape=a.s, rate=b.s)
  fac2=(detdet^0.5)*exp(-0.5*hhh)
  num=fac1*fac2
  den=(det((1/taue)*diag(k1)+U11)^0.5)*((2*pi)^(0.5*(k1-r)))
  return(num/den)
}


ftauetaus.log.present=function(taue, taus, XT, mat1, V.curl)
{
  
  Sig.gamma=as.matrix(rbind(cbind(taue*t(XT)%*%XT,taue*t(XT)%*%t(K), matrix(0, r, k1-r)), 
                            cbind(rbind(taue*K%*%XT, matrix(0, k1-r, r)), taue*diag(k1)+taus*V.curl)))
  Sig.gamma.inv=chol2inv(chol(Sig.gamma))
  Sig.gamma.inv=as.Symmetric.matrix(Sig.gamma.inv)
  
  mat2=diag(k1)- taue*t(rbind(t(XT), mat1))%*%Sig.gamma.inv%*%rbind(t(XT), mat1)
  hhh=taue*t(nucthat)%*%mat2%*%nucthat
  
  t1=dgamma(taue, shape=a.e, rate=b.e, log=T) + dgamma(taus, shape=a.s, rate=b.s, log=T)
  t2=(0.5*k1)*log(taue) + (0.5*(k1-G))*log(taus) 
  
  detdet=det(Sig.gamma)
  if(detdet==0) t3= 0.5*log(10^(-5)) + 0.5*hhh 
  if(detdet!=0) t3= 0.5*log(detdet) + 0.5*hhh 
  
  return(t1+t2-t3)
}

ftauetaus.log.absent=function(taue, taus, XT, L, V.curl)
{
  
  Sig.gamma=as.matrix(rbind(cbind(taue*t(XT)%*%XT,matrix(0, r, k1-r)), 
                            cbind(matrix(0, k1-r, r), taue*diag(k1-r)+taus*V.curl)))
  Sig.gamma.inv=chol2inv(chol(Sig.gamma))
  Sig.gamma.inv=as.Symmetric.matrix(Sig.gamma.inv)
  
  mat2=diag(k1)- taue*t(rbind(t(XT), L))%*%Sig.gamma.inv%*%rbind(t(XT), L)
  hhh=taue*t(nucthat)%*%mat2%*%nucthat
  
  t1=dgamma(taue, shape=a.e, rate=b.e, log=T) + dgamma(taus, shape=a.s, rate=b.s, log=T)
  t2=(0.5*k1)*log(taue) + (0.5*(k1-r-G))*log(taus) 
  
  detdet=det(Sig.gamma)
  if(detdet==0) t3= 0.5*log(10^(-5)) + 0.5*hhh 
  if(detdet!=0) t3= 0.5*log(detdet) + 0.5*hhh 
  
  return(t1+t2-t3)
}


fmarg.unconstrained=function(param, XT, K, L, V.curl)
{
  taue = param[1]; taus = param[2]
  A11 = taue*t(XT)%*%XT ; A12 = cbind(taue*t(XT)%*%t(K), matrix(0, r, k1-r))
  A21 = t(A12); A22 = taue*diag(k1)+taus*V.curl
  Sig.gamma=as.matrix(rbind(cbind(A11, A12), 
                            cbind(A21, A22)))
  Sig.gamma.inv=solve(Sig.gamma) #chol2inv(chol(Sig.gamma))
  
  Sig.gamma.inv=as.Symmetric.matrix(Sig.gamma.inv)
  mat2=diag(k1)- taue*t(rbind(t(XT), K, L))%*%Sig.gamma.inv%*%rbind(t(XT), K, L)
  hhh=taue*t(nucthat)%*%mat2%*%nucthat
  
  t1=dgamma(taue, shape=a.e, rate=b.e) * dgamma(taus, shape=a.s, rate=b.s)
  t2=(taue^(0.5*k1))* (taus^(0.5*(k1-G)))
  t3 = (det(Sig.gamma)^(-0.5))*exp(-0.5*hhh)
  return(t1 * t2 * t3)
}

fmarg.constrained=function(param, XT, L, V22.curl)
{
  taue = param[1]; taus = param[2]
  Sig.gamma=as.matrix(rbind(cbind(taue*t(XT)%*%XT,matrix(0, r, k1-r)), 
                            cbind(matrix(0, k1-r, r), taue*diag(k1-r)+taus*V22.curl)))
  Sig.gamma.inv=chol2inv(chol(Sig.gamma))
  Sig.gamma.inv=as.Symmetric.matrix(Sig.gamma.inv)
  mat2=diag(k1)- taue*t(rbind(t(XT), L))%*%Sig.gamma.inv%*%rbind(t(XT), L)
  hhh=taue*t(nucthat)%*%mat2%*%nucthat
  
  t1=dgamma(taue, shape=a.e, rate=b.e) * dgamma(taus, shape=a.s, rate=b.s)
  t2=(taue^(0.5*k1))* (taus^(0.5*(k1-r-G)))
  t3 = (det(Sig.gamma)^(-0.5))*exp(-0.5*hhh)
  return(t1 * t2 * t3)
}

qstar = function(Q)
{
  eigQ=eigen(Q, symmetric = T)
  Z=eigQ$vectors
  D.diag=D.Q.star=eigQ$values
  D.Q.star[D.Q.star<0]=10^(-8)
  Q.star=Z%*%diag(D.Q.star)%*%t(Z)
  Q.star = as.Symmetric.matrix(Q.star)
  return(Q.star)
}


mcse = function(x, FUNC,  batchsize=floor(sqrt(length(x)))) {
  n = length(x)
  nb = n - batchsize + 1
  bval = rep(NA,nb)
  for (j in 1:nb) {
    ind = j:(batchsize+j-1)
    bval[j] = FUNC(x[ind])
  }
  ##  var.bval = var(bval)*(nb-1) * n * batchsize / ( (n-batchsize) * nb )    #  OLBM method
  var.bval = var(bval)*(nb-1) * batchsize / nb
  
  list(se=sqrt(var.bval / n), batchsize=batchsize)
}

medQuantile = function(z) 
{  
  quantile(z, probs=0.5)  
} # end of fuction 'lowerQuantile'

lowerQuantile = function(z) 
{  
  quantile(z, probs=0.025)  
} # end of fuction 'lowerQuantile'

upperQuantile = function(z) 
{ 
  quantile(z, probs=0.975)  
} # end of fuction 'upperQuantile'

Mode <- function(x){ x = x[!is.na(x)]; ux <- unique(x);  ux[which.max(tabulate(match(x, ux)))]}
width = function(x, pl, pu, na.rm = T, type = 1){quantile(x, probs = pu, na.rm = T, type = 1) - quantile(x, probs = pl, na.rm = T, type = 1)}
as.Symmetric.matrix = function(A)
{ 
  if(!is.symmetric.matrix(A)) {A[upper.tri(A, diag = F)]=t(A)[upper.tri(t(A), diag=F)]}
  return(A)
}
trim.fn=function(x, frac, tail)
{
  if(tail == 'both')
  {
    qs=quantile(x, probs=c(frac, 1-frac), na.rm=T)
    x1=x[x>qs[1] & x<qs[2]]
    x2=c(x1, rep(NA, length(x)-sum(x>qs[1] & x<qs[2])))
  }
  if(tail == 'left')
  {
    qs=quantile(x, probs=frac, na.rm=T)
    x1=x[x>qs]
    x2=c(x1, rep(NA, length(x)-sum(x>qs)))
  }
  if(tail == 'right')
  {
    qs=quantile(x, probs=1-frac, na.rm=T)
    x1=x[x<qs]
    x2=c(x1, rep(NA, length(x)-sum(x<qs)))
  }
  return(x2)
}
