
rm(list=ls())

options(digits=15)
ndraws = 50000
burnin = 20000

source('utility_functions.R')
source('data.R')

V=Vcros

A_HB=seAhb=rep(0,12)

ts = format(Sys.time(), "%d%m%y_%H%M%S")
folderName = paste0("constrained_predictive_hb_", ts, sep = "")
dir.create(path = folderName)

start.time=Sys.time()

post.out=NULL

M=12

lamnew.est=matrix(NA, nrow=p, ncol=M)
lamnew.se=matrix(NA, nrow=p, ncol=M)
muhb.est=matrix(NA, nrow=p, ncol=M)
tot.abundance.hb=matrix(NA, nrow=M, ncol=2)
 i=6
  
  print(paste0('Model ',i))
  
  r=dim(XX[[i]])[2]
  
  X=as.matrix(XX[[i]])
  XT=as.matrix(X[1:k1,])
  XS=as.matrix(X[(k1+1):p,])
   
  B=as.Symmetric.matrix(X %*% solve(t(X)%*% X) %*% t(X)) 
  
  Bc=diag(p)-B
  eigBc=eigen(Bc)
  L=t(eigBc$vectors)
  L=L[eigBc$values >0.9, ] 

  U0=as.Symmetric.matrix(t(L)%*%ginv(L%*% V %*% t(L))%*%L)
  U011=U0[1:k1, 1:k1]
  
  taue=1
  taus=1
  beta=rep(1, r)
  
  #---------------------------------------------------
  # .... assign values for adaptive Metropolis-Hastings sampling
  batch = 0
  batchsize = 5000
  minepsLogSigmaForRW = 0.01
  eps = minepsLogSigmaForRW
  
  sigmaOfProposal.taue = 0.1
  sigmaOfProposal.taus = 0.1
  
  a.s=0.01; b.s=0.01
  a.e=0.01; b.e=0.01
  
  target.curr=ftauetaus.pred(taue, taus, U011)
  #--------------------------
  # Compute MH draws
  #--------------------------
  
  cat('Begin MCMC', '\n', '\n')
  
  continueMCMC = TRUE
  draw = 0
  naccept.taue=0
  naccept.taus=0
  
  fname1=paste0(folderName, '/markovchain_M', i, '.txt', sep="")
  cat(c('abundance', 'taue', 'taus', paste('beta.', 0:(r-1), sep='')), sep=',', file=fname1, append=TRUE)
  
  fname11=paste0(folderName, '/markovchain_lambda_M', i, '.txt', sep="")
  cat(c(paste('lambda', 1:p, sep='')), sep=',', file=fname11, append=TRUE)
 
  while(continueMCMC) 
  {
    
    draw = draw + 1
    drawinterval = batchsize
    
    if (draw == round(draw/drawinterval)*drawinterval)  cat('..... drawing sample #', draw, '\n')
   
    # draw nu
    
    U=(1/taus)*U0
    U11=U[1:k1, 1:k1]
    U12=U[1:k1, (k1+1):p]
    U21=t(U12) 
    U22=U[(k1+1):p, (k1+1):p]
    mat3=solve((1/taue)*diag(k1)+U11)
    
    mu=X%*%beta+rbind(U11, U21)%*%mat3%*%(nucthat-XT%*%beta)
    P=as.Symmetric.matrix(U-rbind(U11, U21)%*%mat3%*%cbind(U11, U12))
    
    nu=as.vector(rmvnorm(1, mu, P, method=c( "svd")))
    lam=exp(nu)
    abundance=sum(lam)
    
    # draw beta
    
    Sig0=(  t(XT)%*%mat3%*%XT )
    Sig.beta= as.Symmetric.matrix(solve(Sig0))
    
    beta0hat=Sig.beta%*%t(XT)%*%mat3%*%nucthat
    beta=as.vector(rmvnorm(1, beta0hat, Sig.beta))
    
    # draw taue
    
    taue.cand=rlnorm(1, log(taue), sigmaOfProposal.taue)
    
    target.cand=ftauetaus.pred(taue.cand, taus, U011)
    num=target.cand * dlnorm(taue,      log(taue.cand), sigmaOfProposal.taue) 
    den=target.curr * dlnorm(taue.cand, log(taue),      sigmaOfProposal.taue) 
    
    if(den==0)
    {
      taue=taue.cand
      target.curr=target.cand
      naccept.taue=naccept.taue+1
    }
    if(den>0)
    {
      R=num/den
      if (runif(1,0,1) <= R) 
      {
        taue=taue.cand
        target.curr=target.cand
        naccept.taue=naccept.taue+1
      }  
    }
    
    if (floor(draw/batchsize)*batchsize == draw) 
    {
      SigmaDiff = ifelse(naccept.taue > 0.44*batchsize, exp(2*eps), exp(-2*eps))
      sigmaOfProposal.taue = sigmaOfProposal.taue * SigmaDiff
      cat("naccept.taue = ", naccept.taue, "proposal sd of taue = ", sigmaOfProposal.taue, '\n')
      naccept.taue = 0 # reset counter for next batch
    }
    
    # draw taus
    
    taus.cand=rlnorm(1, log(taus), sigmaOfProposal.taus)
    
    target.cand=ftauetaus.pred(taue, taus.cand, U011)
    num=target.cand * dlnorm(taus,      log(taus.cand),  sigmaOfProposal.taus) 
    den=target.curr * dlnorm(taus.cand, log(taus),       sigmaOfProposal.taus) 
    
    if(den==0)
    {
      taus=taus.cand
      target.curr=target.cand
      naccept.taus=naccept.taus+1
    }
    if(den>0)
    {
      R=num/den
      if (runif(1,0,1) <= R) 
      {
        taus=taus.cand
        target.curr=target.cand
        naccept.taus=naccept.taus+1
      }  
    }
    
    if (floor(draw/batchsize)*batchsize == draw) 
    {
      SigmaDiff = ifelse(naccept.taus > 0.44*batchsize, exp(2*eps), exp(-2*eps))
      sigmaOfProposal.taus = sigmaOfProposal.taus * SigmaDiff
      cat("naccept.taus = ", naccept.taus, "proposal sd of taus = ", sigmaOfProposal.taus, '\n')
      naccept.taus = 0 # reset counter for next batch
    }
    
    cat('\n', file=fname1, append=TRUE)
    cat(c(abundance, taue, taus, beta), sep=',', file=fname1, append=TRUE)
    
    cat('\n', file=fname11, append=TRUE)
    cat(lam, sep=',', file=fname11, append=TRUE)
    
    
    if (draw == ndraws) 
    {    
      cat('Completed ', ndraws, ' draws of MCMC algorithm', '\n')
      end.time=Sys.time()
      time_taken=end.time-start.time
      print(time_taken)
      # numOfDraws = as.integer(readline(prompt='Enter additional number of MCMC draws -> '))
      numOfDraws =0      
      if (numOfDraws == 0) continueMCMC = FALSE
      if (numOfDraws > 0) ndraws = ndraws + numOfDraws
    }
    
  }  # end of while loop
  
  cat('MCMC is completed!', '\n', '\n')
  
  
  #-------------------------
  # Obtaining Estimates
  #-------------------------
  
  post1 = as.matrix(read.csv(fname1, sep=",", header = T))
  ndraws = dim(post1)[1]
  post1=post1[(burnin+1):ndraws,]
  post = apply(post1, 2, trim.fn, frac = 0.1, tail = 'both')
  
  prob.quantiles = c(0.025, 0.5, 0.975)  
  prob.names = paste(as.character(100*prob.quantiles), '%', sep='')
  post.hb = cbind(apply(post,2,mean,na.rm=T), apply(post,2,sd,na.rm=T), apply(post,2,Mode), 
                  t(apply(post, 2, quantile, probs=prob.quantiles,na.rm=T)), 
                  apply(post, 2, width, pl = prob.quantiles[1], pu = prob.quantiles[2],na.rm=T) )
  
  taue_hb=post.hb[2,1]
  taus_hb=post.hb[3,1]
  beta_hb=post.hb[-c(1:3),1]
 
  U=(1/taus_hb)*U0
  U11=U[1:k1, 1:k1]
  U12=U[1:k1, (k1+1):p]
  U21=t(U12) 
  U22=U[(k1+1):p, (k1+1):p]
  mat4=solve((1/taue_hb)*diag(k1)+U11)
  
  mu_hb=X%*%beta_hb+rbind(U11, U21)%*%mat4%*%(nucthat-XT%*%beta_hb)
  P_hb=as.Symmetric.matrix(U-rbind(U11, U21)%*%mat4%*%cbind(U11, U12))
  muhb.est[,i]=mu_hb
 
  post.lam1 = read.csv(fname11, sep=",", header = T)
  post.lam1=post.lam1[(burnin+1):ndraws,]
  post.lam = apply(post.lam1, 2, trim.fn, frac = 0.1, tail = 'both')
  lamnew.est[,i]=apply(post.lam, 2, mean)
  lamnew.se[,i]=apply(post.lam, 2, sd)
    A_HB[i]=post.hb[1,1]
    seAhb[i]=post.hb[1,2]
  
  cat('HB abundance estimate -> \n')
  tot.abundance.hb[i,]=c(A_HB[i], seAhb[i])
  print(tot.abundance.hb[i,])
  cat('\n')
  
  post.hb=round(post.hb,3)
  dimnames(post.hb) = list(c('abundance', 'taue', 'taus', paste('beta_', 0:(r-1), sep='')),
                           c('Mean', 'SD', 'Mode', prob.names, 'CI Width'))

  post.out=rbind(post.out, c('', '', '', '','', '',''), 
                 c(paste0('Model', i, sep=""), '', '', '','', '',''),  
                 c('Mean', 'SD', 'Mode', prob.names, 'CI Width'), 
                 cbind(post.hb))
 
  end.time=Sys.time()
  (time.taken=end.time-start.time); print(time.taken)

write.csv(post.out, file=paste0(folderName, '/EstimatesOfParameters.csv', sep=''), quote = F, row.names=T)
bf.out=cbind(A_HB, seAhb) 
dimnames(bf.out)=list(as.character(c(1:12)), c('Abundance estimate.hb', 'SE.hb'))

print(bf.out)
fname6=paste0(folderName, "/abundance_estimates.csv", sep = "")
write.csv(bf.out, file=fname6, quote = F,row.names=T) 

#================================================================

A=cbind(indx,c(lamcthat, rep('', p-k1)), lam_os, exp(muhb.est))
A=rbind(A,  
        c('', sum(lamcthat), sum(lam_os), tot.abundance.hb[,1]),
        c('', '', '', apply(exp(muhb.est), 2, sum)))
dimnames(A)[2]=list(c('CT sites', 'lamcthat', 'lamoshat', paste('mu.hb_M', 1:M, sep='')))
fname4=paste0(folderName, '/ComparingSitewiseLamEst.From.mu_hb.csv', sep="")
write.csv(A, file=fname4, quote = F,row.names=F)

A=cbind(indx,c(lamcthat, rep('', p-k1)), lam_os, lamnew.est)
A=rbind(A, 
        c('', sum(lamcthat), sum(lam_os), tot.abundance.hb[,1]),
        c('', '', '', apply(lamnew.est, 2, sum)))
dimnames(A)[2]=list(c('CT sites', 'lamcthat', 'lamoshat', paste('M', 1:M, sep='')))
fname4=paste0(folderName, '/ComparingSitewiseLamEst.csv', sep="")
write.csv(A, file=fname4, quote = F,row.names=F)

