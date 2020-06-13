library(devtools)
library(SuperLearner)

source(transport_sdesie_highdimM.R)
source(myslfunction.R)

n<-10000 

w0<-rbinom(n, 1, .5)
w1<-rbinom(n,1, .4 + (.2*w0))

probsel<-plogis(-1+ log(4)*w1 + log(4)*w0)
psel<-rbinom(n, 1, probsel)

probsite<-plogis(log(1.2)*w1 + log(1.2)*w0 + log(1.2)*w0*w1)
site<-rbinom(n, 1, probsite)

#instrument
a<-rbinom(n, 1, .5)
  
#exposure
z<-rbinom(n,1,plogis(-log(2) + log(4)*a - log(2)*w1 + log(1.4)*site + log(1.43)*site*a))

#mediator
m<-rbinom(n, 1, plogis(-log(2) + log(4)*z - log(1.4)*w1 + log(1.4)*site))

#outcomes
y<-rbinom(n,1, plogis(-log(5) + log(8)*z + log(4)*m - log(1.2)*w1 + log(1.2)*w1*z) )

dat<-data.frame(W0=w0, W1=w1, A=a, Z=z, M=m, Y=y, S=site, psel=psel, probsel=probsel)
obsdat<-dat[dat$psel==1,]

#obsdat$weights <- mean(1 - obsdat$S) / (mean((1 - obsdat$S) / obsdat$psel) * obsdat$psel)
 
candidates<-c('SL.glm.interaction', 'SL.glm', 'SL.mean') 
mediation(data=obsdat, candidates=candidates, nfolds=5, family.outcome="binomial") 
