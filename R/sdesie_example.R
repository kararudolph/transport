############################################################################
# This script estimates SDE and SIE for each site and also estimates 
#   transported SDE and SIE for target site (S=0)
# Uses functions from present paper Rudolph et al. Transporting 
#   Stochastic Direct and Indirect Effects to New Populations. 2019.
# Uses functions from Rudolph et al. Robust and Flexible Estimation of 
#   Stochastic Mediation Effects: A Proposed Method and Example in a 
#   Randomized Trial Setting. Epidemiologic Methods 2017.
############################################################################

##### Load required packages and functions  #####
library(devtools)
source_url("https://raw.githubusercontent.com/kararudolph/SDE-SIE/master/ivmedtmle.R")
source("transportmedtmle.R") 

##### Create data #####
set.seed(2350) 
n <-5000 

w0 <-rbinom(n, 1, .5)
w1 <-rbinom(n,1, .4 + (.2*w0))

probsel <- plogis(-1+ log(4)*w1 + log(4)*w0)
psel <- rbinom(n, 1, probsel)
svywt <- mean(probsel)/probsel

probsite <- plogis(-1+ log(8)*w1 + log(8)*w0 + log(8)*w0*w1)
site <- rbinom(n, 1, probsite)

#instrument
a <- rbinom(n, 1, .5)
  
#exposure
z0 <- rbinom(n,1,plogis(       - log(2)*w1 + log(2)*site)) 
z1 <- rbinom(n,1,plogis(log(4) - log(2)*w1+ log(8)*site))
z <- ifelse(a==1, z1, z0)

#mediator
m0<-rbinom(n, 1, plogis(-log(3)          - log(1.4)*w1 + log(1.4)*site)) 
m1<-rbinom(n, 1, plogis(-log(3) + log(10)- log(1.4)*w1 + log(1.4)*site + log(4)*site))
m<-ifelse(z==1, m1, m0)

#outcomes
y<-rbinom(n,1, plogis(log(1.2)  + (log(3)*z)  + log(3)*m - log(1.2)*w1 + log(1.2)*w1*z) )

dat<-data.frame(w0=w0, w1=w1, a=a, z=z, m=m, y=y, s=site, psel=psel, svywt=svywt, radid_person=seq(1,n,1))
obsdat<-dat[dat$psel==1,]

############################################################################
# Estimate SDE and SIE for each site
# Uses functions from Rudolph et al. Robust and Flexible Estimation of 
#   Stochastic Mediation Effects: A Proposed Method and Example in a 
#   Randomized Trial Setting. Epidemiologic Methods 2017.
############################################################################

site0dat<-obsdat[obsdat$s==0,]
site1dat<-obsdat[obsdat$s==1,] 
zmodel<-"z ~ a + w1 "
mmodel<-"m ~ z + w1 "
ymodel<-"y ~ m + z*w1"
qmodel<-"w1"

# For S=0
zfit<-glm(formula=zmodel, family="binomial", data=site0dat)
mfit<-glm(formula=mmodel, family="binomial", data=site0dat)

za0<-predict(zfit, newdata=data.frame(w1=site0dat$w1, a=0), type="response")
za1<-predict(zfit, newdata=data.frame(w1=site0dat$w1, a=1), type="response")

mz1<-predict(mfit, newdata=data.frame(w1=site0dat$w1, z=1), type="response")
mz0<-predict(mfit, newdata=data.frame(w1=site0dat$w1, z=0), type="response")

gm<-(mz1*za0) + (mz0*(1-za0))
gma1<-(mz1*za1) + (mz0*(1-za1))

ressite0<-ivmedtmle(a=site0dat$a, z=site0dat$z, m=site0dat$m, y=site0dat$y, 
                    w=data.frame(w1=site0dat$w1), svywt=site0dat$svywt, 
                    zmodel=zmodel, mmodel=mmodel, ymodel=ymodel, 
                    qmodel=qmodel, gm=gm, gma1=gma1)

# For S=1
zfit<-glm(formula=zmodel, family="binomial", data=site1dat)
mfit<-glm(formula=mmodel, family="binomial", data=site1dat)

za0<-predict(zfit, newdata=data.frame(w1=site1dat$w1, a=0), type="response")
za1<-predict(zfit, newdata=data.frame(w1=site1dat$w1, a=1), type="response")

mz1<-predict(mfit, newdata=data.frame(w1=site1dat$w1, z=1), type="response")
mz0<-predict(mfit, newdata=data.frame(w1=site1dat$w1, z=0), type="response")

gm<-(mz1*za0) + (mz0*(1-za0))
gma1<-(mz1*za1) + (mz0*(1-za1))

ressite1<-ivmedtmle(a=site1dat$a, z=site1dat$z, m=site1dat$m, y=site1dat$y, 
                    w=data.frame(w1=site1dat$w1), svywt=site1dat$svywt, 
                    zmodel=zmodel, mmodel=mmodel, ymodel=ymodel, 
                    qmodel=qmodel, gm=gm, gma1=gma1)

############################################################################
# Estimate transported SDE and SIE for S=0
# Uses functions from Rudolph et al. Transporting Stochastic Direct and 
#   Indirect Effects to New Populations. 2019.
############################################################################

colnames(obsdat)<-c("W0", "W1", "A", "Z", "M", "Y", "S", "psel", "weights", "radid_person")

forms=list(Sform = formula("S ~  W1*W0"),
           Aform = formula("A ~ 1"),
           ZformStratS = formula("Z ~ A + W1"),
           Zstarform = formula("Z ~ A*S + W1 "),
           QZform = "W1*S",
           MformStratS = formula("M ~  Z+ W1"),
           Mstarform = formula("M ~  Z*S + W1"),
           Yform = formula("Y ~ M + Z*W1"))

Wnames=c("W1")

gstarMpooled<-get_gstarM(obsdat, forms, Wnames, pooled=TRUE, gstar_S=0)
gstarMnp<-get_gstarM(obsdat, forms, Wnames, pooled=FALSE, gstar_S=0)

# Model 1
mod1_pooled <- transportmedtmle(obsdat, forms, Wnames, iv=TRUE, gma1=gstarMpooled[[1]], gma0=gstarMpooled[[2]])
mod1_notpooled <- transportmedtmle(obsdat, forms, Wnames, iv=TRUE, gma1=gstarMnp[[1]], gma0=gstarMnp[[2]])

# Model 2
mod2_pooled <- transportmedtmle(obsdat, forms, Wnames, iv=FALSE, gma1=gstarMpooled[[1]], gma0=gstarMpooled[[2]])
mod2_notpooled <- transportmedtmle(obsdat, forms, Wnames, iv=FALSE, gma1=gstarMnp[[1]], gma0=gstarMnp[[2]])
