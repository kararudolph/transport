#' Transported intent-to-treat average treatment effect
#'
#' More details to be added.
#'
#' @param a needs to be named a and have values 0/1
#' @param z matrix, each column be an individual z covariate and can be either
#'   binary or continuous (or a mix of both types) this code calculates the
#'   ratio dga1s0/dga1s1 with a trick in Zheng(2012)
#' @param y needs to be named y and have values 0/1
#' @param site The site, value 0 for the site where the outcome data is not used
#'   and value 1 for the site where the outcome data is used.
#' @param w in a dataframe named w and with names w1:wx
#' @param weights vector with length = nrow(data). p(Delta_2 | A, W, S)
#' @param aamodel The aamodel
#' @param asitemodel The asitemodel
#' @param s_awz_model formula for fitting S ~ A + W + Z
#' @param s_aw_model formula for fitting S ~ A + W
#' @param aoutmodel The aoutmodel
#' @param aq2model The aq 2 model
#
#' @return A list with three elements:
#'
#' \begin{enumerate}
#' \item est: the parameter estimate.
#' \item var: the variance estimate of the EIC.
#' \item eic: the efficient influnce curve unit-level values.
#' \end{enumerate}
#'
#' @examples TBD.
#'
#' @references TBD.
#'
#' @importFrom stats coef glm plogis predict var
#' @export
transport_ittate =
  function(a, z, y, site, w, weights = NULL,
           superlearner_lib = NULL,
           aamodel,
           lasitemodel,
           s_awz_model,
           s_aw_model,
           aoutmodel,
           aq2model) {
    datw  <- w
    n.dat <- nrow(datw)

    # calculate components of clever covariate
    cpa_glm = glm(formula = aamodel, family = "binomial", data = data.frame(cbind(datw, a = a)))
    cpa <- predict(cpa_glm, newdata = datw, type="response")

    cps_glm = glm(formula = asitemodel, data = data.frame(cbind(site = site, datw)), family = "binomial")
    cps <- predict(cps_glm, type="response")

    #
    # ratio component 1
    # ---------------------------------------------------------------------------------------
    awz <- data.frame(cbind(site, a, w, z))
    s_awz <- glm(formula = s_awz_model, family = 'binomial', data = awz)

    a1wz   <- a0wz <- awz
    # data: intervention A = 0
    a0wz$a <- 0
    # data: intervention A = 1
    a1wz$a <- 1

    # S^hat @ A = 1
    s_pred_a1_wz <- predict(s_awz, newdata = a1wz, type="response")
    s_pred_1_a1_wz <- s_pred_a1_wz
    s_pred_0_a1_wz <- 1 - s_pred_a1_wz

    # S^hat @ A = 0
    s_pred_a0_wz <- predict(s_awz, newdata = a0wz, type="response")
    s_pred_1_a0_wz <- s_pred_a0_wz
    s_pred_0_a0_wz <- 1 - s_pred_a0_wz

    #
    # ratio component 2
    # ---------------------------------------------------------------------------------------
    aw <- data.frame(cbind(site, a, w))
    s_aw <- glm(formula = s_aw_model, family = 'binomial', data = aw)

    a1w   <- a0w <- aw
    a0w$a <- 0
    a1w$a <- 1

    s_pred_a1_w <- predict(s_aw, newdata = a1w, type="response")
    s_pred_1_a1_w <- s_pred_a1_w
    s_pred_0_a1_w <- 1 - s_pred_a1_w

    s_pred_a0_w <- predict(s_aw, newdata = a0w, type="response")
    s_pred_1_a0_w <- s_pred_a0_w
    s_pred_0_a0_w <- 1 - s_pred_a0_w

    # replaces dga1s0/dga1s1
    gz_ratio_a1 <- s_pred_0_a1_wz / s_pred_1_a1_wz * s_pred_1_a1_w / s_pred_0_a1_w
    # replaces dga0s0/dga0s1
    gz_ratio_a0 <- s_pred_0_a0_wz / s_pred_1_a0_wz * s_pred_1_a0_w / s_pred_0_a0_w

    #
    #calculate clever covariate
    # ---------------------------------------------------------------------------------------
    g0w    <- (1-cpa)*(gz_ratio_a0)*(cps/(1-cps))
    g1w    <- cpa*(gz_ratio_a1)*(cps/(1-cps))
    h0w    <- ((1-a)*I(site==1))/g0w
    h1w    <- (a*I(site==1))/g1w

    ymodel <- glm(formula = aoutmodel,
                  family = "binomial",
                  data = data.frame(cbind(datw, a = a, z, site = site, y = y)),
                  subset = site == 1)

    #
    #initial prediciton
    # ---------------------------------------------------------------------------------------
    q      <- cbind(predict(ymodel, type = "link",
                            newdata = data.frame(cbind(datw, a = a, z))),
                    predict(ymodel, type = "link",
                            newdata = data.frame(cbind(datw, a = 0, z))),
                    predict(ymodel, type = "link",
                            newdata = data.frame(cbind(datw, a = 1, z)))
                    )

    epsilon <- coef(glm(y ~ -1 + offset(q[, 1]) + h0w + h1w , family = "binomial",
                        subset = site == 1 ))

    #
    #update initial prediction
    # ---------------------------------------------------------------------------------------
    q1          <- q + c((epsilon[1]*h0w + epsilon[2]*h1w), epsilon[1]/g0w, epsilon[2]/g1w)

    predmodela0 <- suppressWarnings(glm(formula=paste("plogis(q1)", aq2model, sep="~"), data=data.frame(cbind(w, a=a,site=site, q1=q1[,2])), subset=site==0 & a==0 , family="binomial"))
    predmodela1 <- suppressWarnings(glm(formula=paste("plogis(q1)", aq2model, sep="~"), data=data.frame(cbind(w,a=a, site=site, q1=q1[,3])), subset=site==0 & a==1 , family="binomial"))
    predmodelaa <- suppressWarnings(glm(formula=paste("plogis(q1) ~", aq2model, "+a", sep=""), data=data.frame(cbind(w, site=site, q1=q1[,1], a=a)), subset=site==0, family="binomial"))

    #
    #get initial prediction for second regression model
    # ---------------------------------------------------------------------------------------
    q2pred   <- cbind(predict(predmodelaa, type="link", newdata=data.frame(cbind(datw, a=a))),
     predict(predmodela0, type="link", newdata=datw),
      predict(predmodela1, type="link", newdata=datw))

    # clever covariates about site
    cz       <- cbind(ifelse(a==0,I(site==0)/(1-cpa), I(site==0)/cpa), # A = a
                      I(site==0)/(1-cpa), # a = 0
                      I(site==0)/cpa # a = 1
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

    # ITTATE @ site = 1
    tmleest <- mean(plogis(q2[,3][site==0])) - mean(plogis(q2[,2][site==0]))

    ps0     <- mean(I(site==0))

    # inverse weight by the censoring probability p(Delta_2 | A,W,S)
    if (!is.null(weights)) {
        # if(nrow(wmat) != length(weights)) warning("length of weights not match nrow(wmat)!")
        ps0 <- ps0 * weights
    }

    # y - E[Y|A = a, W]
    eic     <- (((h1w/ps0) - (h0w/ps0))*(y - plogis(q[,1]))) +
    # y - E[Y|A = a, W]
        (((a*cz[,3]/ps0) - ((1-a)*cz[,2]/ps0))* (plogis(q[,1]) - plogis(q2pred[,1]))) +
        ((I(site==0)/ps0)*((plogis(q2pred[,3]) - plogis(q2pred[,2])) - tmleest))

    results = list("est" = tmleest,
                   "var" = var(eic) / n.dat,
                   "eic" = eic)
    return(results)

}

