## a variable needs to be named a and have values 0/1
## site variable needs to be named 'site' and needs to have value 0 for the site where the outcome data is not used and value 1 for the site where the outcome data is used
## z variable needs to be named z and have values 0/1
## y variable needs to be named y and have values 0/1
## w variables in a dataframe named w and with names w1:wx

#-----------------------------------------------------------------------------------------
#' fit ittate_TMLE with SuperLearner
#' @importFrom stats coef glm plogis predict var
#' @importFrom SuperLearner SuperLearner All
#' @export
old_ittatetmle_multi_z_sl <-
  function(a, z, y, site, w, aamodel,
           asitemodel,s_awz_model, s_aw_model, aoutmodel, aq2model) {
    datw  <- w
    n.dat <- nrow(datw)

    #calculate components of clever covariate
    cpa <- predict(glm(formula=aamodel, family="binomial", data=data.frame(cbind(datw, a=a))), newdata=datw, type="response")
    cps <- predict(glm(formula=asitemodel, data=data.frame(cbind(site=site, datw)), family="binomial"), type="response")

    # ratio component 1
    # library(SuperLearner)
    # SL.library <- c("SL.glmnet", "SL.glm", "SL.gam", "SL.mean")
    # awz <- data.frame(cbind(a, w, z))
    # s_awz <- SuperLearner::SuperLearner(Y = site, X = awz, family = 'binomial', SL.library = SL.library)
    awz <- data.frame(cbind(site, a, w, z))
    s_awz <- glm(formula = s_awz_model, family = 'binomial', data = awz)

    a1wz   <- a0wz <- awz
    a0wz$a <- 0
    a1wz$a <- 1

    # s_pred_a1_wz <- predict(s_awz, newdata = a1wz, type="response")[[1]]
    s_pred_a1_wz <- predict(s_awz, newdata = a1wz, type="response")
    s_pred_1_a1_wz <- s_pred_a1_wz
    s_pred_0_a1_wz <- 1 - s_pred_a1_wz

    # s_pred_a0_wz <- predict(s_awz, newdata = a0wz, type="response")[[1]]
    s_pred_a0_wz <- predict(s_awz, newdata = a0wz, type="response")
    s_pred_1_a0_wz <- s_pred_a0_wz
    s_pred_0_a0_wz <- 1 - s_pred_a0_wz

    # ratio component 2
    # aw <- data.frame(cbind(a, w))
    # s_aw <- SuperLearner::SuperLearner(Y = site, X = aw, family = 'binomial', SL.library = SL.library)
    aw <- data.frame(cbind(site, a, w))
    s_aw <- glm(formula = s_aw_model, family = 'binomial', data = aw)

    a1w   <- a0w <- aw
    a0w$a <- 0
    a1w$a <- 1

    # s_pred_a1_w <- predict(s_aw, newdata = a1w, type="response")[[1]]
    s_pred_a1_w <- predict(s_aw, newdata = a1w, type="response")
    s_pred_1_a1_w <- s_pred_a1_w
    s_pred_0_a1_w <- 1 - s_pred_a1_w

    # s_pred_a0_w <- predict(s_aw, newdata = a0w, type="response")[[1]]
    s_pred_a0_w <- predict(s_aw, newdata = a0w, type="response")
    s_pred_1_a0_w <- s_pred_a0_w
    s_pred_0_a0_w <- 1 - s_pred_a0_w

    # replaces dga1s0/dga1s1
    gz_ratio_a1 <- s_pred_0_a1_wz / s_pred_1_a1_wz * s_pred_1_a1_w / s_pred_0_a1_w
    # replaces dga0s0/dga0s1
    gz_ratio_a0 <- s_pred_0_a0_wz / s_pred_1_a0_wz * s_pred_1_a0_w / s_pred_0_a0_w
    #calculate clever covariate
    g0w    <- (1-cpa)*(gz_ratio_a0)*(cps/(1-cps))
    g1w    <- cpa*(gz_ratio_a1)*(cps/(1-cps))
    h0w    <- ((1-a)*I(site==1))/g0w
    h1w    <- (a*I(site==1))/g1w

    ymodel <- glm(formula=aoutmodel,
                  family="binomial",
                  data=data.frame(cbind(datw, a=a, zmat, site=site, y=y)),
                  subset=site==1)

    #initial prediction
    q      <- cbind(predict(ymodel, type="link", newdata=data.frame(cbind(datw, a=a,zmat))),
                    predict(ymodel, type="link", newdata=data.frame(cbind(datw, a=0,zmat))),
                    predict(ymodel, type="link", newdata=data.frame(cbind(datw, a=1,zmat)))
                    )

    epsilon <- coef(glm(y ~ -1 + offset(q[,1]) + h0w + h1w , family="binomial", subset=site==1 ))

    #update initial prediction
    q1          <- q + c((epsilon[1]*h0w + epsilon[2]*h1w), epsilon[1]/g0w, epsilon[2]/g1w)

    predmodela0 <- suppressWarnings(glm(formula=paste("plogis(q1)", aq2model, sep="~"), data=data.frame(cbind(w, a=a,site=site, q1=q1[,2])), subset=site==0 & a==0 , family="binomial"))
    predmodela1 <- suppressWarnings(glm(formula=paste("plogis(q1)", aq2model, sep="~"), data=data.frame(cbind(w,a=a, site=site, q1=q1[,3])), subset=site==0 & a==1 , family="binomial"))
    predmodelaa <- suppressWarnings(glm(formula=paste("plogis(q1) ~", aq2model, "+a", sep=""), data=data.frame(cbind(w, site=site, q1=q1[,1], a=a)), subset=site==0, family="binomial"))

    #get initial prediction for second regression model
    q2pred   <- cbind(predict(predmodelaa, type="link", newdata=data.frame(cbind(datw, a=a))), predict(predmodela0, type="link", newdata=datw), predict(predmodela1, type="link", newdata=datw))

    cz       <- cbind(ifelse(a==0,I(site==0)/(1-cpa), I(site==0)/cpa),
                      I(site==0)/(1-cpa),
                      I(site==0)/cpa
                      )

    epsilon2 <- suppressWarnings(
        coef(glm(plogis(q1[,1]) ~ -1 + offset(q2pred[,1]) + cz[,2] + cz[,3], family="binomial", subset= site==0))
        )
    for(k in 1:2){
        epsilon2[k] <- ifelse(is.na(epsilon2[k]), 0, epsilon2[k])
    }

    q2      <- q2pred + c((epsilon2[1]*cz[,2] + epsilon2[2]*cz[,3]),
                         epsilon2[1]/(1-cpa),
                         epsilon2[2]/cpa)

    tmleest <- mean(plogis(q2[,3][site==0])) - mean(plogis(q2[,2][site==0]))

    ps0     <- mean(I(site==0))

    eic     <- (((h1w/ps0) - (h0w/ps0))*(y - plogis(q[,1]))) +
        (((a*cz[,3]/ps0) - ((1-a)*cz[,2]/ps0))* (plogis(q[,1]) - plogis(q2pred[,1]))) +
        ((I(site==0)/ps0)*((plogis(q2pred[,3]) - plogis(q2pred[,2])) - tmleest))

    return(list("est"=tmleest, "var"=var(eic)/n.dat, "eic"=eic))

}

## a variable needs to be named a and have values 0/1
## site variable needs to be named 'site' and needs to have value 0 for the site where the outcome data is not used and value 1 for the site where the outcome data is used
## z variable needs to be named z and have values 0/1
## y variable needs to be named y and have values 0/1
## w variables in a dataframe named w and with names w1:wx

old_ittatetmle_multi_z <- function(a, z, y, site, w, aamodel, asitemodel,s_awz_model, s_aw_model, aoutmodel, aq2model){
    datw  <- w
    n.dat <- nrow(datw)

    #calculate components of clever covariate
    cpa <- predict(glm(formula=aamodel, family="binomial", data=data.frame(cbind(datw, a=a))), newdata=datw, type="response")
    cps <- predict(glm(formula=asitemodel, data=data.frame(cbind(site=site, datw)), family="binomial"), type="response")

    # ratio component 1
    # library(SuperLearner)
    # SL.library <- c("SL.glmnet", "SL.glm", "SL.gam", "SL.mean")
    # awz <- data.frame(cbind(a, w, z))
    # s_awz <- SuperLearner::SuperLearner(Y = site, X = awz, family = 'binomial', SL.library = SL.library)
    awz <- data.frame(cbind(site, a, w, z))
    s_awz <- glm(formula = s_awz_model, family = 'binomial', data = awz)

    a1wz   <- a0wz <- awz
    a0wz$a <- 0
    a1wz$a <- 1

    # s_pred_a1_wz <- predict(s_awz, newdata = a1wz, type="response")[[1]]
    s_pred_a1_wz <- predict(s_awz, newdata = a1wz, type="response")
    s_pred_1_a1_wz <- s_pred_a1_wz
    s_pred_0_a1_wz <- 1 - s_pred_a1_wz

    # s_pred_a0_wz <- predict(s_awz, newdata = a0wz, type="response")[[1]]
    s_pred_a0_wz <- predict(s_awz, newdata = a0wz, type="response")
    s_pred_1_a0_wz <- s_pred_a0_wz
    s_pred_0_a0_wz <- 1 - s_pred_a0_wz

    # ratio component 2
    # aw <- data.frame(cbind(a, w))
    # s_aw <- SuperLearner::SuperLearner(Y = site, X = aw, family = 'binomial', SL.library = SL.library)
    aw <- data.frame(cbind(site, a, w))
    s_aw <- glm(formula = s_aw_model, family = 'binomial', data = aw)

    a1w   <- a0w <- aw
    a0w$a <- 0
    a1w$a <- 1

    # s_pred_a1_w <- predict(s_aw, newdata = a1w, type="response")[[1]]
    s_pred_a1_w <- predict(s_aw, newdata = a1w, type="response")
    s_pred_1_a1_w <- s_pred_a1_w
    s_pred_0_a1_w <- 1 - s_pred_a1_w

    # s_pred_a0_w <- predict(s_aw, newdata = a0w, type="response")[[1]]
    s_pred_a0_w <- predict(s_aw, newdata = a0w, type="response")
    s_pred_1_a0_w <- s_pred_a0_w
    s_pred_0_a0_w <- 1 - s_pred_a0_w

    # replaces dga1s0/dga1s1
    gz_ratio_a1 <- s_pred_0_a1_wz / s_pred_1_a1_wz * s_pred_1_a1_w / s_pred_0_a1_w
    # replaces dga0s0/dga0s1
    gz_ratio_a0 <- s_pred_0_a0_wz / s_pred_1_a0_wz * s_pred_1_a0_w / s_pred_0_a0_w
    #calculate clever covariate
    g0w    <- (1-cpa)*(gz_ratio_a0)*(cps/(1-cps))
    g1w    <- cpa*(gz_ratio_a1)*(cps/(1-cps))
    h0w    <- ((1-a)*I(site==1))/g0w
    h1w    <- (a*I(site==1))/g1w

    ymodel <- glm(formula=aoutmodel,
                  family="binomial",
                  data=data.frame(cbind(datw, a=a, zmat, site=site, y=y)),
                  subset=site==1)

    #initial prediciton
    q      <- cbind(predict(ymodel, type="link", newdata=data.frame(cbind(datw, a=a,zmat))),
                    predict(ymodel, type="link", newdata=data.frame(cbind(datw, a=0,zmat))),
                    predict(ymodel, type="link", newdata=data.frame(cbind(datw, a=1,zmat)))
                    )

    epsilon <- coef(glm(y ~ -1 + offset(q[,1]) + h0w + h1w , family="binomial", subset=site==1 ))

    #update initial prediction
    q1          <- q + c((epsilon[1]*h0w + epsilon[2]*h1w), epsilon[1]/g0w, epsilon[2]/g1w)

    predmodela0 <- suppressWarnings(glm(formula=paste("plogis(q1)", aq2model, sep="~"), data=data.frame(cbind(w, a=a,site=site, q1=q1[,2])), subset=site==0 & a==0 , family="binomial"))
    predmodela1 <- suppressWarnings(glm(formula=paste("plogis(q1)", aq2model, sep="~"), data=data.frame(cbind(w,a=a, site=site, q1=q1[,3])), subset=site==0 & a==1 , family="binomial"))
    predmodelaa <- suppressWarnings(glm(formula=paste("plogis(q1) ~", aq2model, "+a", sep=""), data=data.frame(cbind(w, site=site, q1=q1[,1], a=a)), subset=site==0, family="binomial"))

    #get initial prediction for second regression model
    q2pred   <- cbind(predict(predmodelaa, type="link", newdata=data.frame(cbind(datw, a=a))), predict(predmodela0, type="link", newdata=datw), predict(predmodela1, type="link", newdata=datw))

    cz       <- cbind(ifelse(a==0,I(site==0)/(1-cpa), I(site==0)/cpa),
                      I(site==0)/(1-cpa),
                      I(site==0)/cpa
                      )

    epsilon2 <- suppressWarnings(
        coef(glm(plogis(q1[,1]) ~ -1 + offset(q2pred[,1]) + cz[,2] + cz[,3], family="binomial", subset= site==0))
        )
    for(k in 1:2){
        epsilon2[k] <- ifelse(is.na(epsilon2[k]), 0, epsilon2[k])
    }

    q2      <- q2pred + c((epsilon2[1]*cz[,2] + epsilon2[2]*cz[,3]),
                         epsilon2[1]/(1-cpa),
                         epsilon2[2]/cpa)

    tmleest <- mean(plogis(q2[,3][site==0])) - mean(plogis(q2[,2][site==0]))

    ps0     <- mean(I(site==0))

    eic     <- (((h1w/ps0) - (h0w/ps0))*(y - plogis(q[,1]))) +
        (((a*cz[,3]/ps0) - ((1-a)*cz[,2]/ps0))* (plogis(q[,1]) - plogis(q2pred[,1]))) +
        ((I(site==0)/ps0)*((plogis(q2pred[,3]) - plogis(q2pred[,2])) - tmleest))

    return(list("est"=tmleest, "var"=var(eic)/n.dat, "eic"=eic))

}

