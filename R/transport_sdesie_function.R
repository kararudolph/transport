#this function uses TMLE to estimate the transported stochastic direct and indirect effect 
#for observed data O=(S,W,A,Z,M,Y), where
#S=1 is the source population and S=0 is the target population
#W are covariates
#A is an instrumental variable
#Z is an intermediate variable that is affected by instrument A
#M is a mediator variable and a function of W, Z under Model 1
# and a function of W, Z, A under Model 2
#Y is an outcome variable and a function of W, Z, M under Model 1
# and a function of W, Z, M, A under Model 2 
# Model 1 includes the exclusion restriction assumptions (there is no direct effect of A on M or of A on Y)
# Model 2 allows the exclusion restriction not to hold

get_gstarM  = function(data, forms, Wnames, pooled, gstar_S) {  
  dataZ1 = dataZ0 = dataA1 = dataA0 = data
  # Define the parameter as to the mechanism used for M and Z
  
  dataZ1$Z = 1
  dataZ0$Z = 0
  dataA1$A = 1
  dataA0$A = 0
  
  if (!pooled) {
    Mstarfit = glm(formula=forms$MformStratS, data=data[data$S==gstar_S, ], family = "binomial")
    Zstarfit = glm(formula=forms$ZformStratS, data=data[data$S==gstar_S, ], family = "binomial")
    
    predMz1 = predict(Mstarfit, newdata = dataZ1, type = 'response')
    predMz0 = predict(Mstarfit, newdata = dataZ0, type = 'response')
    
    predZa0 = predict(Zstarfit, newdata = dataA0, type = 'response')
    predZa1 = predict(Zstarfit, newdata = dataA1, type = 'response')
    
  } else {
    Mstarfit = glm(formula=forms$Mstarform, data=data, family = "binomial")
    Zstarfit = glm(formula=forms$Zstarform, data=data, family = "binomial")
    
    dataZ1$S = gstar_S
    dataZ0$S = gstar_S
    dataA1$S = gstar_S
    dataA0$S = gstar_S
    
    predMz1 = predict(Mstarfit, newdata = dataZ1, type = 'response')
    predMz0 = predict(Mstarfit, newdata = dataZ0, type = 'response')
    
    predZa0 = predict(Zstarfit, newdata = dataA0, type = 'response')
    predZa1 = predict(Zstarfit, newdata = dataA1, type = 'response')
    
  }
  
  gstarM_astar0 = predMz1*predZa0 + predMz0*(1 - predZa0)
  gstarM_astar1 = predMz1*predZa1 + predMz0*(1 - predZa1)
  
  return(list(gstarM_astar1 = gstarM_astar1, gstarM_astar0 = gstarM_astar0))
}


transportmedtmle<-function(data, forms, Wnames, iv, gma1, gma0){
  
  #get inital fit Q_Y
  yfit<-glm(formula=forms$Yform, family="binomial", data=data[data$S==1,])
  qyinit<-cbind(predict(yfit, newdata=data, type="response"), 
                predict(yfit, newdata=data.frame(cbind(data[,c(Wnames, "A", "Z")], M=0)), type="response"),
                predict(yfit, newdata=data.frame(cbind(data[,c(Wnames, "A", "Z")], M=1)), type="response"))
  
  #estimate weights for targeting
  afit<-glm(formula=forms$Aform, family="binomial", data=data)
  psa1s1<-predict(afit, newdata=data.frame(cbind(data[,Wnames], S=1)), type="response")
  psa1<-predict(afit, newdata=data, type="response")
  
  mz<-predict(glm(formula=forms$MformStratS, family="binomial", data=data[data$S==1,]), newdata=data, type="response")
  psm<-(mz*data$M) + ((1-mz)*(1-data$M))
  
  ps1w<-predict(glm(formula=forms$Sform, family="binomial", data=data), newdata=data, type="response")
  ps0<-mean(1-data$S) 
  
  zfit<-glm(formula=forms$Zstarform, family="binomial", data=data)
  pzs0<-predict(zfit, newdata=data.frame(cbind(data[,Wnames, drop = FALSE], A=data$A, S=0)), type="response")
  pz<-predict(zfit, newdata=data, type="response")
  pza1s0<-predict(zfit, newdata=data.frame(cbind(data[,Wnames, drop = FALSE], A=1, S=0)), type="response")
  pza0s0<-predict(zfit, newdata=data.frame(cbind(data[,Wnames, drop=FALSE], A=0, S=0)), type="response")
  pza1s1<-predict(zfit, newdata=data.frame(cbind(data[,Wnames, drop=FALSE], A=1, S=1)), type="response")
  pza0s1<-predict(zfit, newdata=data.frame(cbind(data[,Wnames, drop=FALSE], A=0, S=1)), type="response")
  pzs1<-pza1s1*psa1s1 + pza0s1*(1-psa1s1)
  
  if(iv==TRUE){
    ha1gma1<-(((data$M*gma1 + (1-data$M)*(1-gma1))* ( (data$Z*pza1s0) + ((1-data$Z)*(1-pza1s0)) ) *(1-ps1w)) / (psm* ( (data$Z*pzs1) + ((1-data$Z)*(1-pzs1)) )*ps1w*ps0))*I(data$S==1)*data$weights
    ha1gma0<-(((data$M*gma0 + (1-data$M)*(1-gma0))* ( (data$Z*pza1s0) + ((1-data$Z)*(1-pza1s0)) ) *(1-ps1w)) / (psm* ( (data$Z*pzs1) + ((1-data$Z)*(1-pzs1)) )*ps1w*ps0))*I(data$S==1 )*data$weights
    ha0gma0<-(((data$M*gma0 + (1-data$M)*(1-gma0))* ( (data$Z*pza0s0) + ((1-data$Z)*(1-pza0s0)) ) *(1-ps1w)) / (psm* ( (data$Z*pzs1) + ((1-data$Z)*(1-pzs1)) )*ps1w*ps0))*I(data$S==1 )*data$weights
  }
  else{
    ha1gma1<-(((data$M*gma1 + (1-data$M)*(1-gma1))* ( (data$Z*pza1s0) + ((1-data$Z)*(1-pza1s0)) ) *(1-ps1w)) / (psm* ( (data$Z*pzs1) + ((1-data$Z)*(1-pzs1)) ) * psa1s1*ps1w*ps0))*I(data$S==1 & data$A==1)*data$weights
    ha1gma0<-(((data$M*gma0 + (1-data$M)*(1-gma0))* ( (data$Z*pza1s0) + ((1-data$Z)*(1-pza1s0)) ) *(1-ps1w)) / (psm* ( (data$Z*pzs1) + ((1-data$Z)*(1-pzs1)) ) * psa1s1*ps1w*ps0))*I(data$S==1 & data$A==1)*data$weights
    ha0gma0<-(((data$M*gma0 + (1-data$M)*(1-gma0))* ( (data$Z*pza0s0) + ((1-data$Z)*(1-pza0s0)) ) *(1-ps1w)) / (psm* ( (data$Z*pzs1) + ((1-data$Z)*(1-pzs1)) ) * (1-psa1s1)*ps1w*ps0))*I(data$S==1 & data$A==0)*data$weights
  }
  
  #target Q_Y
  #for E(Y_{1,gma0})
  epsilonma1g0<-coef(glm(Y ~  1 , weights=ha1gma0, offset=(qlogis(qyinit[,1])), family="quasibinomial", data=data))
  qystara1gma0<-cbind(plogis(qlogis(qyinit[,2]) + epsilonma1g0), plogis(qlogis(qyinit[,3]) + epsilonma1g0))
  
  #for E(Y_{1,gma1})
  epsilonma1g1<-coef(glm(Y ~  1 , weights=ha1gma1, offset=(qlogis(qyinit[,1])), family="quasibinomial", data=data))
  qystara1gma1<-cbind(plogis(qlogis(qyinit[,2]) + epsilonma1g1), plogis(qlogis(qyinit[,3]) + epsilonma1g1))
  
  #for E(Y_{0,gma0})
  epsilonma0g0<-coef(glm(Y ~  1 , weights=ha0gma0, offset=(qlogis(qyinit[,1])), family="quasibinomial", data=data))
  qystara0gma0<-cbind(plogis(qlogis(qyinit[,2]) + epsilonma0g0), plogis(qlogis(qyinit[,3]) + epsilonma0g0))
  
  #estimate Q_M
  data$Qma1g0<-qystara1gma0[,1]*(1-gma0) + qystara1gma0[,2]*gma0
  data$Qma1g1<-qystara1gma1[,1]*(1-gma1) + qystara1gma1[,2]*gma1
  data$Qma0g0<-qystara0gma0[,1]*(1-gma0) + qystara0gma0[,2]*gma0
  
  #estimate Q_Z
  Qzfita1g0<-glm(formula=paste("Qma1g0", forms$QZform, sep="~"), data=data[data$A==1,], family="quasibinomial")
  Qzfita1g1<-glm(formula=paste("Qma1g1", forms$QZform, sep="~"), data=data[data$A==1,], family="quasibinomial")
  Qzfita0g0<-glm(formula=paste("Qma0g0", forms$QZform, sep="~"), data=data[data$A==0,], family="quasibinomial")
  
  Qza1g0<-predict(Qzfita1g0, type="response", newdata=data)
  Qza1g1<-predict(Qzfita1g1, type="response", newdata=data)
  Qza0g0<-predict(Qzfita0g0, type="response", newdata=data)
  
  #update Q_Z 
  #Note: only need to do the update step if A is nonrandom
  ha1<-(I(data$S==0 & data$A==1)/(psa1*ps0))*data$weights
  ha0<-(I(data$S==0 & data$A==0)/((1-psa1)*ps0))*data$weights
  
  epsilonza1g0<-coef(glm(Qma1g0~ 1 , data=data, weights=ha1, offset=qlogis(Qza1g0), family="quasibinomial"))
  epsilonza1g1<-coef(glm(Qma1g1~ 1 , data=data, weights=ha1, offset=qlogis(Qza1g1), family="quasibinomial"))
  epsilonza0g0<-coef(glm(Qma0g0~ 1 , data=data, weights=ha0, offset=qlogis(Qza0g0), family="quasibinomial"))
  
  Qzupa1g0<-plogis(qlogis(Qza1g0) + epsilonza1g0)
  Qzupa1g1<-plogis(qlogis(Qza1g1) + epsilonza1g1)
  Qzupa0g0<-plogis(qlogis(Qza0g0) + epsilonza0g0)
  
  #estimate psi
  tmlea1m0<-sum(Qzupa1g0*data$weights*I(data$S==0))/sum(data$weights[data$S==0])
  tmlea1m1<-sum(Qzupa1g1*data$weights*I(data$S==0))/sum(data$weights[data$S==0])
  tmlea0m0<-sum(Qzupa0g0*data$weights*I(data$S==0))/sum(data$weights[data$S==0])
  
  #estimands
  nde<-tmlea1m0-tmlea0m0
  
  nie<-tmlea1m1-tmlea1m0
  
  scaling <- sum(data$S==0)/sum(data$weights[data$S==0])
  
  #EIC
  D_Y_11 <- scaling*ha1gma1*(data$Y - (qystara1gma1[,1]*(1-data$M) + qystara1gma1[,2]*data$M))
  D_Y_10 <- scaling*ha1gma0*(data$Y - (qystara1gma0[,1]*(1-data$M) + qystara1gma0[,2]*data$M))
  D_Y_00 <- scaling*ha0gma0*(data$Y - (qystara0gma0[,1]*(1-data$M) + qystara0gma0[,2]*data$M))
  
  D_Z_11 <- scaling*ha1*(data$Qma1g1 - Qzupa1g1)
  D_Z_10 <- scaling*ha1*(data$Qma1g0 - Qzupa1g0)
  D_Z_00 <- scaling*ha0*(data$Qma0g0 - Qzupa0g0)
  
  D_W_11 <- (Qzupa1g1*data$weights - tmlea1m1*(sum(data$weights[data$S==0])/nrow(data[data$S==0,]))) * I(data$S==0)/ps0
  D_W_10 <- (Qzupa1g0*data$weights - tmlea1m0*(sum(data$weights[data$S==0])/nrow(data[data$S==0,]))) * I(data$S==0)/ps0
  D_W_00 <- (Qzupa0g0*data$weights - tmlea0m0*(sum(data$weights[data$S==0])/nrow(data[data$S==0,]))) * I(data$S==0)/ps0
  
  nde_eic<- (D_Y_10 + D_Z_10 + D_W_10) - (D_Y_00 + D_Z_00 + D_W_00)  
  nie_eic<- (D_Y_11 + D_Z_11 + D_W_11) - (D_Y_10 + D_Z_10 + D_W_10) 
  
  #sample variance
  ndevar<-var(nde_eic)/nrow(data)
  nievar<-var(nie_eic)/nrow(data)
  
  return(list("nde"=nde, "nie"=nie, "ndevar"=ndevar, "nievar"=nievar))
}
