## This first function is from Alex Luedtke's paper: 
## Luedtke, Alex, Marco Carone, and Mark J. van der Laan. 
## "An omnibus non‐parametric test of equality in distribution for unknown functions." 
## Journal of the Royal Statistical Society: Series B (Statistical Methodology) 81.1 (2019): 75-99.

mmd.test = function(R,S,D.R,D.S,sig.meth='eig',num.reps=1e4,return.cutoff=FALSE){
  n = length(R)
  
  D.R.mat1 = matrix(rep(D.R,n),nrow=n)
  D.R.mat2 = matrix(rep(D.R,each=n),nrow=n)
  
  D.S.mat1 = matrix(rep(D.S,n),nrow=n)
  D.S.mat2 = matrix(rep(D.S,each=n),nrow=n)
  
  R.mat1 = matrix(rep(R,n),nrow=n)
  R.mat2 = matrix(rep(R,each=n),nrow=n)
  
  S.mat1 = matrix(rep(S,n),nrow=n)
  S.mat2 = matrix(rep(S,each=n),nrow=n)
  
  EE = ((2*(R.mat1-R.mat2)*(D.R.mat2-D.R.mat1) + 1 - (4*(R.mat1-R.mat2)^2-2)*D.R.mat1*D.R.mat2)*exp(-(R.mat1-R.mat2)^2)
        - ((2*(S.mat1-R.mat2)*(D.R.mat2-D.S.mat1) + 1 - (4*(S.mat1-R.mat2)^2-2)*D.S.mat1*D.R.mat2)*exp(-(S.mat1-R.mat2)^2))
        - ((2*(R.mat1-S.mat2)*(D.S.mat2-D.R.mat1) + 1 - (4*(R.mat1-S.mat2)^2-2)*D.R.mat1*D.S.mat2)*exp(-(R.mat1-S.mat2)^2))
        + (2*(S.mat1-S.mat2)*(D.S.mat2-D.S.mat1) + 1 - (4*(S.mat1-S.mat2)^2-2)*D.S.mat1*D.S.mat2)*exp(-(S.mat1-S.mat2)^2))
  
  # EE = exp(-(R.mat1-R.mat2)^2) - 2*exp(-(S.mat1-R.mat2)^2) + exp(-(S.mat1-S.mat2)^2)
  
  if(sig.meth=='eig'){
    line.means = rowMeans(EE)
    EE.ctrd = EE - matrix(rep(line.means,n),nrow=n) - matrix(rep(line.means,each=n),nrow=n) + matrix(rep(mean(line.means),n^2),nrow=n)
    num.eigs = min(200,n)
    # tmp = eigen(EE.ctrd)$values/n
    require('rARPACK')
    tmp = eigs_sym(EE.ctrd,num.eigs,which='LA')$values/n
    num.pos.eigs = num.eigs # sum(tmp>0)
    draws=c(matrix(rnorm(num.reps*num.pos.eigs)^2-1,nrow=num.reps,ncol=num.pos.eigs)%*%cbind(tmp[1:num.pos.eigs]))
  }
  
  # U-statistic
  diag(EE) = 0
  est = (rbind(rep(1/(n-1),n)) %*% EE %*% cbind(rep(1/n,n)))[1,1]
  # V-statistic
  # est = (rbind(rep(1/n,n)) %*% EE %*% cbind(rep(1/n,n)))[1,1]
  
  if(sig.meth=='eig'){
    pval = mean(draws>n*est)
  } else if(sig.meth=='var'){
    pval = pchisq(est/(2*var(D.R)/n)+1,df=1,lower.tail=FALSE)
  }
  
  return(if(!return.cutoff){
    c(est,pval)
  }else{
    c(est,pval,
      if(sig.meth=='eig'){
        quantile(draws,0.95)
      }else{
        2*var(D.R)*(qchisq(0.95,df=1)-1)})})
}


# Tests if S is used in the regression function, i.e. if
# E[M|A,W,S] = E[M|A,W] almost surely, or equivalently (under positivity) if
# E[Y|S=1,A,W] = E[Y|S=0,A,W] almost surely
# Returns an estimate of Psi, and a p-value
# If est.g is false, then just uses predefined function g0

est.psi.probM = function(S,W,A,Z,M, astar, sig.meth='var',est.g=TRUE){
  require('SuperLearner')
  #SL.library=c('SL.rpart','SL.glm.interaction','SL.glm','SL.earth','SL.nnet')
  SL.library=c('SL.glmnet', 'SL.glm')
  n=length(S)
  
  # estimate Z regressions for astar and S = 0 and S = 1
  X = data.frame(W=W,A=A,S=S)
  newX = data.frame(W=rbind(W,W), A = astar, S=rep(c(0,1),each=n))
  
  Z.est = SuperLearner(Y=Z,X=X,newX = newX, SL.library=SL.library)
  Z.est.S0 = Z.est$SL.predict[,1][1:n]
  Z.est.S1 = Z.est$SL.predict[,1][(n+1):(2*n)]
  
  X = data.frame(W=W,Z=Z,S=S)
  newX = data.frame(W=rbind(W,W,W,W),Z=rep(c(0,1,0,1),each=n), S = rep(c(0,0,1,1), each = n))
  # Estimate outcome regressions
  M.est = SuperLearner(Y=M, X=X, newX=newX, SL.library=SL.library)
  M.est.Z0S0 = M.est$SL.predict[,1][1:n]
  M.est.Z1S0 = M.est$SL.predict[,1][1:n+n]
  M.est.Z0S1 = M.est$SL.predict[,1][1:n+2*n]
  M.est.Z1S1 = M.est$SL.predict[,1][1:n+3*n]
  
  M.est.0 = M.est.Z0S0*(1 - Z.est.S0) + M.est.Z1S0*Z.est.S0
  M.est.1 = M.est.Z0S1*(1 - Z.est.S1) + M.est.Z1S1*Z.est.S1
  
  if(est.g){
    gg = SuperLearner(Y=S,X=W,SL.library=SL.library,family='binomial')$SL.predict[,1]
  } else {
    gg = g0(W)
  }
  
  # Plug-in estimate of blip
  R = M.est.1 - M.est.0
  SS = rep(0,n)
  
  D.R = S/gg * (M-M.est.1) - (1-S)/(1-gg) * (M-M.est.0)
  D.S = rep(0,n)
  
  return(mmd.test(R,SS,D.R,D.S,sig.meth=sig.meth))
}
