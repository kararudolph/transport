#' Omnibus test
#'
#' @importFrom SuperLearner SuperLearner All
#'
#' @author Alex Luedtke
#'
#' @references
#'
#' Luedtke et al. 2017
#'
#' @export
est_psi_prob_binom =
  function(W, A, Y, W.train = NULL, A.train = NULL, Y.train = NULL,
           sig.meth = 'var', est.g = TRUE,
           g0 = NULL,
           SL.library = c('SL.glm', 'SL.step', 'SL.earth', 'SL.gam', 'SL.nnet', 'SL.mean','SL.glm.interaction')) {
  #SL.library=c('SL.rpart','SL.glm.interaction','SL.glm','SL.earth','SL.nnet')

  #SL.library <- c("SL.glm", "SL.step", "SL.earth", "SL.gam", "SL.nnet", "SL.mean", "SL.bayesglm")
  n=length(A)

  if(is.null(W.train) | is.null(A.train) | is.null(Y.train)){
    W.train = W
    A.train = A
    Y.train = Y
  }

  # Estimate outcome regressions
  Qbar.est = SuperLearner(Y=Y.train,X=data.frame(W=W.train,A=A.train),newX=data.frame(W=rbind(W,W),A=rep(c(0,1),each=n)),SL.library=SL.library, family='binomial')
  Qbar.est.0 = Qbar.est$SL.predict[,1][1:n]
  Qbar.est.1 = Qbar.est$SL.predict[,1][(n+1):(2*n)]

  if(est.g){
    gg = SuperLearner(Y=A,X=W,SL.library=SL.library,family='binomial')$SL.predict[,1]
  } else {
    gg = g0(W)
  }

  # Plug-in estimate of blip
  R = Qbar.est.1 - Qbar.est.0
  S = rep(0,n)

  D.R = A/gg * (Y-Qbar.est.1) - (1-A)/(1-gg) * (Y-Qbar.est.0)
  D.S = rep(0,n)

  return(mmd_test(R,S,D.R,D.S,sig.meth=sig.meth))
}
