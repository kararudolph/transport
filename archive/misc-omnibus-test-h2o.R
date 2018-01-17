est.psi.prob.binom.h2o = function(W,A,Y,W.train=NULL,A.train=NULL,Y.train=NULL,sig.meth='var',est.g=TRUE, data,
                                  learner = c("h2o.glm.wrapper", "h2o.randomForest.wrapper", "h2o.gbm.wrapper"),
                                  # learner = c("h2o.glm.wrapper"),
                                  metalearner = "h2o.glm.wrapper"){
  require('h2oEnsemble')
  require('h2o')
  #why do I have to run these outside of the function?
  h2o.gbmint2=function(..., ntrees=100, interaction.depth=2, seed=1)h2o.gbm.wrapper(...,ntrees=ntrees, interaction.depth=interaction.depth, seed=seed)
  h2o.lasso = function(..., alpha=1)h2o.glm.wrapper(...,alpha=alpha)
  h2o.reg = function(..., alpha=.5)h2o.glm.wrapper(...,alpha=alpha)
  h2o.ridge = function(..., alpha=0)h2o.glm.wrapper(...,alpha=alpha)
  learner=c("h2o.glm.wrapper", "h2o.gbm.wrapper", "h2o.gbmint2", "h2o.lasso",  "h2o.reg", "h2o.ridge")
  # learner=c("h2o.glm.wrapper", "h2o.gbm.wrapper", "h2o.gbmint2", "h2o.lasso")
  #learner=c("h2o.gbmint2", "h2o.lasso",  "h2o.reg", "h2o.ridge")
  metalearner="h2o.glm.wrapper"


  tmpdat<-data[,c(A,Y)]
  x<-c(names(W), A)
  dat<-data.frame(cbind(tmpdat, W))
  dat.h2o<-as.h2o(dat)
  dat.h2o[,A] <-as.factor(dat.h2o[,A])
  dat.h2o[,Y] <-as.factor(dat.h2o[,Y])

  dat.h2o.a0<-dat.h2o.a1<-dat.h2o
  dat.h2o.a0[,A]<-0
  dat.h2o.a1[,A]<-1
  n=nrow(dat)

  if(is.null(W.train) | is.null(A.train) | is.null(Y.train)){
    W.train = W
    A.train = A
    Y.train = Y
  }

  # Estimate outcome regressions
  Qbar.est = h2o.ensemble(x=x, y=Y, training_frame=dat.h2o, family="binomial", learner=learner, metalearner=metalearner, cvControl=list(V=5))

  Qbar.est.0 = as.data.frame(predict(Qbar.est, dat.h2o.a0)$pred)[,3]
  Qbar.est.1 = as.data.frame(predict(Qbar.est, dat.h2o.a1)$pred)[,3]

  if(est.g){
    x<-c(names(W))
    gfit<-h2o.ensemble(x=x, y=A, training_frame=dat.h2o, family="binomial", learner=learner, metalearner=metalearner, cvControl=list(V=5))
    gg = as.data.frame(predict(gfit, dat.h2o)$pred)[,3]
  } else {
    gg = g0(W)
  }

  # Plug-in estimate of blip
  R = Qbar.est.1 - Qbar.est.0
  S = rep(0,n)

  D.R = dat[,A]/gg * (dat[,Y]-Qbar.est.1) - (1-dat[,A])/(1-gg) * (dat[,Y]-Qbar.est.0)
  D.S = rep(0,n)

  return(mmd.test(R,S,D.R,D.S,sig.meth=sig.meth))
}