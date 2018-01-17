#' Tests if A is used in the regression function, i.e. if
#' E[Y|A,W] = E[Y|W] almost surely, or equivalently (under positivity) if
#' E[Y|A=1,W] = E[Y|A=0,W] almost surely
#' Returns an estimate of Psi, and a p-value
#' If est.g is false, then just uses predefined function g0
#' NOTE: Theory all uses bounded Y, so if Y is not bounded then at least try to make sure most of its
#'	mass falls in the interval [-1,1] --- can do this by scaling by a constant
#' @importFrom SuperLearner SuperLearner All
#'
#' @author Alex Luedtke
#'
#' @references
#'
#' Luedtke et al. 2017
#'
#' @export
est_psi_prob =
  function(W, A, Y, W.train = NULL, A.train = NULL, Y.train = NULL, sig.meth = 'var',
           est.g = TRUE, g0 = NULL,
           SL.library=c('SL.rpart','SL.glm','SL.earth','SL.nnet','SL.gam')) {
  #SL.library=c('SL.rpart','SL.glm.interaction','SL.glm','SL.earth','SL.nnet')

  n=length(A)

  if(is.null(W.train) | is.null(A.train) | is.null(Y.train)){
    W.train = W
    A.train = A
    Y.train = Y
  }

  # Estimate outcome regressions
  Qbar.est = SuperLearner(Y=Y.train,X=data.frame(W=W.train,A=A.train),
                          newX=data.frame(W=rbind(W,W),A=rep(c(0,1),each=n)),
                          SL.library=SL.library)
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
