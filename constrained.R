
rm(list=ls())

set.seed(2016)
options(digits=15)
ndraws = 50000
burnin = 25000

source('utility_functions.R')
source('data.R')


V=Vcr

G = k1- rankMatrix(V)[1] # number of islands
diag.V = diag(V)

#------------------------------------------------
# Names of the paramters for different models
#------------------------------------------------
i=6
print(paste0('Model ',i))

ts = format(Sys.time(), "%d%m%y_%H%M%S")
folderName = paste0("constrained_M", i, "_", ndraws,"run_", ts, sep = "")
dir.create(path = folderName)

start.time=Sys.time()

post.out=cor.mat=cor.mat.full=NULL

M=12

rr=rep(0, M)
for(jj in 1:M) rr[jj]=dim(XX[[jj]])[2]


  r=dim(XX[[i]])[2]
  
  X=as.matrix(XX[[i]])
  XT=as.matrix(X[1:k1,])
  XS=as.matrix(X[(k1+1):p,])
 
  B=XT %*% solve(t(XT)%*% XT) %*% t(XT) # capital B
  eigB=eigen(B)
  
  if(i == 1) 
  {
    K=as.matrix(t(eigB$vectors))
    K=as.matrix(t(K[eigB$values >0.9, ]))  # eigen vectors corrsp nonzero eigen values (i.e 1) of B
    KT=K[,1:k1 ]
  }
  if(i != 1) 
  {
    K=as.matrix(t(eigB$vectors))
    K=as.matrix(K[eigB$values >0.9, ])  # eigen vectors corrsp nonzero eigen values (i.e 1) of B
    KT=K[,1:k1 ]
  }
 
  Bc=diag(k1)-B
  eigBc=eigen(Bc)
  L=t(eigBc$vectors)
  L=L[eigBc$values >0.9, ] # eigen vectors corrsp nonzero eigen values (i.e 1) of Bc
  LT=L[,1:k1 ]
 
  V.curl=L%*%V%*%t(L)
  
  
  batch = 0
  batchsize = 5000
  minepsLogSigmaForRW = 0.01
  eps = minepsLogSigmaForRW
  
  sigmaOfProposal.taue = 0.2
  sigmaOfProposal.taus = 0.1
  fname1 = paste(folderName, '/sigma_proposals.txt', sep = "")
  cat(c(paste("sigmaOfProposal.taue = ", sigmaOfProposal.taue, sep = ''), '\n'), sep = ',', file = fname1, append = TRUE)
  cat(c(paste("sigmaOfProposal.taus = ", sigmaOfProposal.taus, sep = ''), '\n'), sep = '', file = fname1, append = TRUE)
  
  a.e=0.01; b.e=0.01
  a.s=0.01; b.s=0.01
  
  taue = 1 
  taus = 1
  beta = rep(1, r)
  gamma=c(beta,rnorm(k1-r,mean=0,sd=0.1))
  
  logtaue = log(taue)
  logtaus = log(taus)
  log.target.curr=ftauetaus.log.absent(exp(logtaue), exp(logtaus), XT, L, V.curl)
  #--------------------------
  # Compute MH draws
  #--------------------------
  
  cat('Begin MCMC', '\n', '\n')
  
  continueMCMC = TRUE
  draw = 0
  naccept.tauetaus=0
  
  fname1=paste0(folderName, '/markovchain_M', i, '.txt', sep="")
  cat(c('taue', 'taus', paste('beta.', 1:r, sep='')),  paste('theta2.', 1:(k1-r), sep=''),
      sep=',', file=fname1, append=TRUE)
  
  while(continueMCMC) 
  {
    
    draw = draw + 1
    drawinterval = batchsize
    
    if (draw == round(draw/drawinterval)*drawinterval)  cat('..... drawing sample #', draw, '\n')
    
    # update the increment/decrement for adaptive Metropolis-Hastings samplers
    if (floor(draw/batchsize)*batchsize == draw) 
    {
      batch = batch + 1
      if (1/sqrt(batch) < minepsLogSigmaForRW)  eps = 1/sqrt(batch)
    }
    
    #############################################################
    
    # draw gamma
    
    Sig.gamma=as.matrix(rbind(cbind(taue*t(XT)%*%XT, matrix(0, r, k1-r)), 
                              cbind( matrix(0, k1-r, r),  taue*diag(k1-r)+taus*V.curl)))
    Sig.gamma.inv=solve(Sig.gamma)
    Sig.gamma.inv = as.Symmetric.matrix(Sig.gamma.inv)
    gamma0hat=taue*Sig.gamma.inv%*%rbind(t(XT), L)%*%nucthat
    gamma=as.vector(rmvnorm(1, gamma0hat, Sig.gamma.inv))

    #############################################################
    # draw taue and taus
    
      logtaue.cand=rnorm(1, logtaue, sigmaOfProposal.taue)
      logtaus.cand=rnorm(1, logtaus, sigmaOfProposal.taus)
      
      log.target.cand=ftauetaus.log.absent(exp(logtaue.cand), exp(logtaus.cand), XT, L, V.curl)
      
      lognum=log.target.cand 
      
      logden=log.target.curr 
      
      if(!is.nan(log.target.cand)) c.es=0
    
    if(lognum!=-Inf | logden!=-Inf)
    {
      logR =lognum-logden
      if (runif(1,0,1) <= exp(logR)) 
      {
        logtaue=logtaue.cand
        logtaus=logtaus.cand
        log.target.curr=log.target.cand
        naccept.tauetaus=naccept.tauetaus+1
      }
    }
    
    if (floor(draw/batchsize)*batchsize == draw) 
    {
      SigmaDiff = ifelse(naccept.tauetaus > 0.44*batchsize, exp(2*eps), exp(-2*eps))
      sigmaOfProposal.taue = sigmaOfProposal.taue * SigmaDiff
      sigmaOfProposal.taus = sigmaOfProposal.taus * SigmaDiff
      cat("naccept.tauetaus = ", naccept.tauetaus, '\n',
          "proposal sd of taue = ", sigmaOfProposal.taue,
          "proposal sd of taus = ", sigmaOfProposal.taus,'\n')
      naccept.tauetaus = 0 # reset counter for next batch
    }
   
    #############################################################
    
      cat('\n', file=fname1, append=TRUE)
    cat(c(exp(logtaue), exp(logtaus), gamma), sep=',', file=fname1, append=TRUE)
  
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

  post = as.matrix(read.csv(fname1, sep=",", header = T))
  ndraws = dim(post)[1]
  post.mcmc=as.mcmc(post)
  post=post[(burnin+1):ndraws,]
  
  
  # Compute posterior means and quantiles
 
  prob.quantiles = c(0.025,0.25, 0.5,0.75, 0.975)  # for credible limits

  prob.names = paste(as.character(100*prob.quantiles), '%', sep='')
  
  post.hb = cbind(apply(post,2,mean, na.rm=T), apply(post,2,sd, na.rm=T),
                  apply(post,2,Mode),
                  t(apply(post, 2, quantile, probs=prob.quantiles, na.rm=T)),
                  apply(post,2,width, na.rm = T, pl = prob.quantiles[1], pu = prob.quantiles[5]) )
  
  post.hb=round(post.hb,3)
  dimnames(post.hb) = list(c('taue', 'taus', paste('beta_', 1:r, sep=''),
                             paste('theta2_', 1:(k1-r), sep='')), 
                           c('Mean', 'SD', 'Mode', prob.names, 'CI width'))
  print(post.hb)
  post.out=rbind(post.out, rep('', 9), 
                 c(paste0('Model', i, sep=""), rep('', 8)),  
                 c('Mean', 'SD', 'Mode', prob.names, 'CI width'), 
                 post.hb)
  
  cor.comp=round(cor(post), 5)
  dimnames(cor.comp) = list(c('taue', 'taus', paste('beta.', 1:rr[i], sep=''), 
                              paste('theta2.', 1:(k1-rr[i]), sep='')),
                            c('taue', 'taus', paste('beta.', 1:rr[i], sep=''), 
                              paste('theta2.', 1:(k1-rr[i]), sep='')))

    cor.mat.full=round(cor(post), 5)
 
  fname6=paste0(folderName, "/traceplots_taue_M", i, ".jpeg", sep = "")
  jpeg(fname6, width = 1000, height = 1000, units = "px", pointsize = 12, quality = 100)
  traceplot(post.mcmc[,1], xlab = "Iterations", ylab = "taue", main=paste('Traceplot of taue_M', i, '.txt', sep=""))
  dev.off()
  
  fname6=paste0(folderName, "/traceplots_taus_M", i, ".jpeg", sep = "")
  jpeg(fname6, width = 1000, height = 1000, units = "px", pointsize = 12, quality = 100)
  traceplot(post.mcmc[,2], xlab = "Iterations", ylab = "taus", main=paste('Traceplot of taus_M', i, '.txt', sep=""))
  dev.off()
  
  end.time=Sys.time()
  (time.taken=end.time-start.time); print(time.taken)


write.csv(post.out, file=paste0(folderName, '/EstimatesOfParameters.csv', sep=''), quote = F, row.names=T)
write.csv(cor.mat.full, file=paste0(folderName, '/Correlation.matrices.full.csv', sep=''), quote = F, row.names=T)

