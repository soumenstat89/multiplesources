
rm(list=ls())

set.seed(2016)
options(digits=15)

source('utility_functions.R')
source('data.R')


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
marg_num1=rep(0,12)
a.e=0.01; b.e=0.01
a.s=0.01; b.s=0.01

for( i in 1:M)
{
# i=7
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

V.curl=L%*%V%*%t(L)

topaste=paste0('Model ',i)
print(topaste)

numint=cuhre(ndim=2, ncomp=1, integrand=fmarg, XT = XT, L = L, V.curl = V.curl,
             lower=rep(10^(-10),2), upper=rep(20,2),
             rel.tol=10^(-10), abs.tol = 0,
             flags=list(verbose=0, final=0, pseudo.random=0, mersenne.seed=NULL),
             min.eval=0, max.eval=20000)

marg_num1[i]=numint$value
print(marg_num1[i])
} # end of for(i in 1:M)

BF_num1=marg_num1/max(marg_num1)
bf.out=cbind(marg_num1, BF_num1)
dimnames(bf.out)=list(as.character(c(1:12)), c('marg_cuhre',' BF_cuhre'))
print(bf.out)
fname6=paste(folderName, "/marginal_densities.csv", sep = "")
write.csv(bf.out, file=fname6, quote = F,row.names=T)
