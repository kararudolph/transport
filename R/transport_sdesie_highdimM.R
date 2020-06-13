#this function estimates the transported stochastic direct and indirect effects
#for observed data O=(S,W,A,Z,M,Y), where
#S=1 is the source population and S=0 is the target population
#W are covariates
#A is a treatment or instrument and is 0/1
#Z is an intermediate variable. for this implementation, it is 0/1, but the function may be adapted to allow for alternative distributions
#M are mediating variables that may follow any distribution; they are considered jointly
#Y is an outcome variable and may follow any distribution

#superlearner is used for model fitting
#source myslfunction.R

estimatheta <- function(data, contrast, weights,
                        candidates, 
                        fitg, fite, fitm, fitq, fitr, fitb, 
                        valSets) {
    ## contrasts
    aprime <- contrast[1]
    astar <- contrast[2]

    ## extract data
    A <- data[, "A"]
    M <- data[, substr(names(data), 1, 1) == "M"]
    Z <- data[, substr(names(data), 1, 1) == "Z"]
    Y <- data[, "Y"]
    W <- data[, substr(names(data), 1, 1) == "W"]
	S <- data[, "S"]
    weights <- weights
    ## observations and cross-validation
    n <- length(A)
  
    ## compute u pseudo outcome and fit u function
    gone <- stats::predict(fitg,
                                newdata = data.frame(W, S=0))$pred[, 1]
    eone <- stats::predict(fite,
                                newdata = data.frame(W, M, S=0))$pred[, 1]
    qoneprime <- stats::predict(fitq,
                                newdata = data.frame(W, A = aprime, S=0))$pred[, 1]
    roneprime <- stats::predict(fitr,
                                newdata = data.frame(W, M, A = aprime, S=0))$pred[, 1]
    mprime <- stats::predict(fitm,
                             newdata = data.frame(W, M, Z, A = aprime))$pred[, 1]
    mprimeone <- stats::predict(fitm,
                             newdata = data.frame(W, M, Z=1, A = aprime))$pred[, 1]
    mprimezero <- stats::predict(fitm,
                             newdata = data.frame(W, M, Z=0, A = aprime))$pred[, 1]

    boneprime <- stats::predict(fitb,
                             newdata = data.frame(W, M, Z, A = aprime))$pred[, 1]
    boneprimez1 <- stats::predict(fitb,
                             newdata = data.frame(W, M, Z=1, A = aprime))$pred[, 1]
    boneprimez0 <- stats::predict(fitb,
                             newdata = data.frame(W, M, Z=0, A = aprime))$pred[, 1]

  ## compute u pseudo outcome and fit u function

    gprime <- gone * aprime + (1 - gone) * (1 - aprime)
    gstar <- gone * astar + (1 - gone) * (1 - astar)
    eprime <- eone * aprime + (1 - eone) * (1 - aprime)
    estar <- eone * astar + (1 - eone) * (1 - astar)

    qprime <- Z * qoneprime + (1 - Z) * (1 - qoneprime)
    rprime <- Z * roneprime + (1 - Z) * (1 - roneprime)
  
    hstar  <- (1 - boneprime) / boneprime * gprime / gstar * qprime / rprime * estar / eprime

    upseudo <- mprime * gprime / gstar * qprime / rprime * estar / eprime
 ## in case values are all very close to each other
    if (sd(upseudo) < .Machine$double.eps) {
        candidateu <- c("SL.mean")
    } else {
        candidateu<-candidates
    }

    ## predict on estimated u
    fitu <- mySL(
        Y = upseudo,
        X = data.frame(W, Z, A, S),
        family = stats::gaussian(),
        SL.library = candidateu,
        validRows = valSets
    )

    uprimeone <-
        stats::predict(fitu, newdata = data.frame(W, Z = 1, A = aprime, S=0))$pred[, 1]
    uprimezero <-
        stats::predict(fitu, newdata = data.frame(W, Z = 0, A = aprime, S=0))$pred[, 1]
    
    ## compute v pseudo outcome and fit v function
    ## build v nuisance function
    vpseudo <- mprimeone * qoneprime + mprimezero * (1 - qoneprime)

    if (stats::sd(vpseudo) < .Machine$double.eps) {
        candidatev <- c("SL.mean")
    } else {
        candidatev <- candidates
    }

    ## fit nuisance function v
    fitv <- mySL(
        Y = vpseudo,
        X = data.frame(W, A, S),
        family = stats::gaussian(),
        SL.library = candidatev,
        validRows = valSets
    )
    vstar <- stats::predict(fitv, newdata = data.frame(W, A = astar, S=0))$pred[, 1]

    ## compute one step
    t<-1-mean(S)

    ipwprime1 <- as.numeric(A == aprime & S == 1) / (gprime*t)
    ipwprime0 <- as.numeric(A == aprime & S == 0) / (gprime*t)
    ipwstar0  <- as.numeric(A == astar & S == 0) / (gstar*t)

    ## compute EIF components
    uprimediff <- (uprimeone - uprimezero)
    eify <- ipwprime1 * hstar / mean(ipwprime1 * hstar) * (Y - mprime)
    eifu <- ipwprime0 / mean(ipwprime0) * uprimediff *
        (Z - qoneprime)
    eifv <- ipwstar0 / mean(ipwstar0) * (vpseudo - vstar)
    eif <- weights * (eify + eifu + eifv + (1-S)/t *vstar) 
    os <- mean(eif)

    eifos <- eif - (1-S)/t * os

    ## now, compute the TMLE
    stopcrit <- FALSE
    iter <- 1

    ## iterative TMLE
    while (!stopcrit) {
       
        hstar <- (1 - boneprime) / boneprime * gprime / gstar * qprime / rprime * estar / eprime
        hstarone <- (1 - boneprimez1) / boneprimez1 * gprime / gstar *
            qoneprime / roneprime * estar / eprime
        hstarzero <- (1 - boneprimez0) / boneprimez0* gprime / gstar *
            (1 - qoneprime) / (1 - roneprime) * estar / eprime

         ## first fluctuation/tilting
        suppressWarnings(
            tiltm <- stats::glm(
                as.formula("Y ~ -1 + offset(mprime_logit) + hstar"),
                data = data.frame(list(Y = Y, A = A, S = S,
                                       mprime_logit = stats::qlogis(mprime),
                                       hstar = hstar)),
                subset = A == aprime & S == 1,
                weights = weights / (gprime * t),
                family = stats::binomial()
            )
        )

        ## second fluctuation/tilting
        suppressWarnings(
            tiltq <- stats::glm(
                stats::as.formula("Z ~ -1 + offset(qoneprime_logit) + uprimediff"),
                data = data.frame(list(Z = Z, A = A,
                                       qoneprime_logit = stats::qlogis(qoneprime),
                                       uprimediff = uprimediff)),
                subset = A == aprime & S == 0,
                weights = weights / (gprime * t),
                family = stats::binomial()
            )
        )
        ## extract epsilon and set to zero if failed fluctuation/tilting
        coefq <- stats::coef(tiltq)
        coefm <- stats::coef(tiltm)
        if (is.na(coefq)) coefq <- 0
        if (is.na(coefm)) coefm <- 0

        mprime <- stats::plogis(stats::qlogis(mprime) + coefm * hstar)
        mprimeone <- stats::plogis(stats::qlogis(mprimeone) + coefm * hstarone)
        mprimezero <- stats::plogis(stats::qlogis(mprimezero) + coefm *
                                    hstarzero)

        qoneprime <- stats::plogis(stats::qlogis(qoneprime) + coefq *
                                   uprimediff)
        qprime <- Z * qoneprime + (1 - Z) * (1 - qoneprime)

        hstar <- (1 - boneprime) / boneprime * gprime / gstar * qprime / rprime * estar / eprime
        ipwprimeone  <- as.numeric(A == aprime & S == 1) / (t * gprime)
        ipwprimezero <- as.numeric(A == aprime & S == 0) / (t * gprime)
        ipwstar      <- as.numeric(A == astar & S == 0) / (t * gstar)
        eify <- ipwprimeone * hstar / mean(ipwprimeone * hstar) * (Y - mprime)
        eifu <- ipwprimezero / mean(ipwprimezero) * uprimediff *
            (Z - qoneprime)

        eifyu <- weights * (eify + eifu)

        stopcrit <- abs(mean(eifyu)) < sd(eifyu) / (sqrt(n) * max(10, log(n)))

        ## iterate iterator
        iter <- iter + 1
    }  

    vpseudo <- mprimeone * qoneprime + mprimezero * (1 - qoneprime)

    ## fit nuisance function v
    fitv <- mySL(
        Y = vpseudo,
        X = data.frame(W, A, S),
        family = stats::gaussian(),
        SL.library = candidates,
        validRows = valSets
    )
    vstar <- stats::predict(fitv, newdata = data.frame(W, A = astar, S=0))$pred[, 1]


    vstar[vstar < 1e-3] <- 1e-3
    vstar[vstar > 1 - 1e-3] <- 1 - 1e-3
    suppressWarnings(
        tiltv <- stats::glm(
                            as.formula("vpseudo ~ offset(vstar_logit)"),
                            data = data.frame(list(vpseudo = vpseudo, A = A, S=S, vstar_logit = stats::qlogis(vstar))),
                            subset = (A == astar & S==0),
                            weights = weights / (gprime * t),
                            family = stats::binomial()
                        )
    )
    vstar <- stats::plogis(stats::qlogis(vstar) + stats::coef(tiltv))

    qprime <- Z * qoneprime + (1 - Z) * (1 - qoneprime)
    hstar <- ((1-boneprime) * gprime * qprime * estar) / (boneprime * gstar * rprime * eprime)
    upseudo <- mprime * gprime / gstar * qprime / rprime * estar / eprime
    vpseudo <- mprimeone * qoneprime + mprimezero * (1 - qoneprime)

    ipwprime1 <- as.numeric(A == aprime & S == 1) / (gprime*t)
    ipwprime0 <- as.numeric(A == aprime & S == 0) / (gprime*t)
    ipwstar0  <- as.numeric(A == astar & S == 0) / (gstar*t)

    eify <- ipwprime1 * hstar / mean(ipwprime1 * hstar) * (Y - mprime)
    eifu <- ipwprime0 / mean(ipwprime0) * uprimediff *
        (Z - qoneprime)
    eifv <- ipwstar0  * (vpseudo - vstar)

    #eif <- (sum(I(S==0))/sum(weights[S==0]))*weights * (eify + eifu + eifv + vstar)
#should it be
     eif <- weights * (eify + eifu + eifv + (1 - S) / t * vstar)
    
    #tmle <- sum(vstar*weights*I(S==0))/sum(weights[S==0])
    tmle <- mean(weights * (1 - S) * vstar) / t
    eiftmle <- eif - weights * (1 - S) / t * tmle

    return(list(
        estimates = c(os = os, tmle = tmle),
        eifs = cbind(os = eifos, tmle = eiftmle)
    ))
}

mediation <- function(data, candidates, nfolds, family.outcome) {

    n <- nrow(data)
    
    A <- data[, "A"]
    M <- data[, substr(names(data), 1, 1) == "M"]
    Z <- data[, substr(names(data), 1, 1) == "Z"]
    Y <- data[, "Y"]
    W <- data[, substr(names(data), 1, 1) == "W"]
    S <- data[, "S"]
    weights <- mean(1 - data$S) / (mean((1 - data$S) / data$probsel) * data$probsel)

    ## preliminaries
    
    ## observations and cross-validation
    valSets <- split(sample(seq_len(n)), rep(seq_len(nfolds), length = n))
    #for S=1
    ns1 <- length(A[S==1])
    valSetss1 <- split(sample(seq_len(ns1)), rep(seq_len(nfolds), length = ns1))

    ## fit nuisance functions
    fitg <- mySL(
        Y = A,
        X = data.frame(W, S),
        family = stats::binomial(),
        SL.library = candidates,
        validRows = valSets
    )
    fite <- mySL(
        Y = A,
        X = data.frame(W, M, S),
        family = stats::binomial(),
        SL.library = candidates,
        validRows = valSets
    )

    fitm <- mySL(
        Y = Y[S==1],
        X = subset(data.frame(W, M, Z, A), S==1),
        family = family.outcome,
        SL.library = candidates,
        validRows = valSetss1
    )

    fitq <- mySL(
        Y = Z,
        X = data.frame(W, A, S),
        family = stats::binomial(),
        SL.library = candidates,
        validRows = valSets
    )
    fitr <- mySL(
        Y = Z,
        X = data.frame(W, M, A, S),
        family = stats::binomial(),
        SL.library = candidates,
        validRows = valSets
    )
    fitb <- mySL(
        Y = S,
        X = data.frame(W, M, Z, A),
        family = stats::binomial(),
        SL.library = candidates,
        validRows = valSets
    )

      ## compute under different contrasts
    theta11 <- estimatheta(data,
                           contrast = c(1, 1),
                           weights,
                           candidates, fitg, fite, fitm, fitq, fitr, fitb, valSets
                           )
    theta10 <- estimatheta(data,
                           contrast = c(1, 0), 
                           weights,
                           candidates, fitg, fite, fitm, fitq, fitr, fitb, valSets
                           )
    theta00 <- estimatheta(data,
                           contrast = c(0, 0), 
                           weights,
                           candidates, fitg, fite, fitm, fitq, fitr, fitb, valSets
                           )

    ## point estimates
    indirect <- theta11$estimates - theta10$estimates
    direct <- theta10$estimates - theta00$estimates

    ## standard error estimates
    seindirect <- sqrt(diag(stats::var(theta11$eifs - theta10$eifs)) / n)
    sedirect <- sqrt(diag(stats::var(theta10$eifs - theta00$eifs)) / n)

    ## output
    return(list(indirectest=indirect, directest=direct, indirectse=seindirect, directse=sedirect))
 
}
