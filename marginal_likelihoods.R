
rm(list=ls())

set.seed(2016)
options(digits=15)

source('utility_functions.R')
source('data.R')
# source('marginal_likelihoods.R')

V=Vcr

G = k1- rankMatrix(V)[1] # number of islands
diag.V = diag(V)

ts = format(Sys.time(), "%d%m%y_%H%M%S")
folderName = paste0("marginal_", ts, sep = "")
dir.create(path = folderName)

start.time=Sys.time()

M=12

rr=rep(0, M)
for(i in 1:M) rr[i]=dim(XX[[i]])[2]


#==========================================
# Marginal densities
#==========================================
marg_num1.constrained=marg_num1.unconstrained = rep(0,12)
a.e=0.01; b.e=0.01
a.s=0.01; b.s=0.01

for( i in 1:M){
  # i=3
  r=dim(XX[[i]])[2]
  
  X=as.matrix(XX[[i]])
  XT=as.matrix(X[1:k1,])
  XS=as.matrix(X[(k1+1):p,])
  
  B=as.Symmetric.matrix(XT %*% solve(t(XT)%*% XT) %*% t(XT)) 
  eigB=eigen(B)
  
  if(i == 1) 
  {
    K=as.matrix(t(eigB$vectors))
    K=as.matrix(t(K[eigB$values >0.9, ]))  
    KT=K[,1:k1 ]
  }
  if(i != 1) 
  {
    K=as.matrix(t(eigB$vectors))
    K=as.matrix(K[eigB$values >0.9, ])  
    KT=K[,1:k1 ]
  }
  
  Bc=diag(k1)-B
  eigBc=eigen(Bc)
  L=t(eigBc$vectors)
  L=L[eigBc$values >0.9, ] 
  LT=L[,1:k1 ]
  
  V.curl= rbind(K, L) %*% V %*% t(rbind(K, L))
  
  V22.curl=L%*%V%*%t(L)
  
  numint.unconstrained=cuhre(ndim=2, ncomp=1, integrand=fmarg.unconstrained, XT = XT, K = K, L = L, V.curl = V.curl,
                             lower=rep(10^(-10),2), upper=rep(20,2),
                             rel.tol=10^(-10), abs.tol = 0,
                             flags=list(verbose=0, final=0, pseudo.random=0, mersenne.seed=NULL),
                             min.eval=0) #, max.eval=30000)
  
  marg_num1.unconstrained[i] = numint.unconstrained$value
  print(paste0('Model ',i, ': Unonstrained marginal = ',marg_num1.unconstrained[i]))
  
  numint.constrained=cuhre(ndim=2, ncomp=1, integrand=fmarg.constrained, XT = XT, L = L, V22.curl = V22.curl,
                           lower=rep(10^(-10),2), upper=rep(20,2),
                           rel.tol=10^(-10), abs.tol = 0,
                           flags=list(verbose=0, final=0, pseudo.random=0, mersenne.seed=NULL),
                           min.eval=0, max.eval=20000)
  
  marg_num1.constrained[i]=numint.constrained$value
  print(paste0('Model ',i, ': Constrained marginal = ',marg_num1.constrained[i]))
  
} # end of for(i in 1:M)

BF_num1.unconstrained=marg_num1.unconstrained/max(marg_num1.unconstrained)
BF_num1.constrained=marg_num1.constrained/max(marg_num1.constrained)
bf.out=cbind(marg_num1.unconstrained, BF_num1.unconstrained, marg_num1.constrained, BF_num1.constrained)
dimnames(bf.out)=list(as.character(c(1:12)),
                      c('marg_cuhre.unconstrained',' BF_cuhre.unconstrained', 'marg_cuhre.constrained',' BF_cuhre.constrained'))

print(bf.out)
fname6=paste(folderName, "/marginal_densities.csv", sep = "")
write.csv(bf.out, file=fname6, quote = F,row.names=T)
