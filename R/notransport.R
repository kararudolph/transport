#' No-transport estimator
#'
#' This estimator is also used in \code{\link{transport_cace}}.
#'
#' <Parameter info>
#'
#' <Return info>
#'
#' @examples
#'
#' # TBD.
#'
#' @references TBD.
#'
#' @importFrom stats coef glm plogis predict var
#' @export
notransport <- function(a, z, w, site, ntamodel, ntzmodel) {
  datw <- w
  n.dat <- nrow(datw)
  ps0 <- mean(I(site == 0))

  # Calculate components of clever covariate.
  glm_cpa = glm(formula = ntamodel,
                data = data.frame(cbind(a = a, site = site, datw)),
        subset = site == 0, family = "binomial")
  cpa <- predict(glm_cpa, newdata = datw, type = "response")

  g0w <- 1 - cpa
  g1w <- cpa

  # Clever covariates.
  h0w <- I(site == 0) * (1 - a) / (g0w * ps0)
  h1w <- I(site == 0) * a / (g1w * ps0)

  zmodel <- glm(formula = ntzmodel, family = "binomial",
                data = data.frame(cbind(a = a, z = z , site = site , datw)),
    subset = site == 0)

  data_new0 <- data_new1 <- datw
  data_new0$a <- 0
  data_new1$a <- 1

  q <- cbind(predict(zmodel, type = "link",
                     newdata = data.frame(cbind(a = a, z = z, datw))),
      predict(zmodel, type = "link", newdata = data_new0),
      predict(zmodel, type = "link", newdata = data_new1))

  epsilon <- coef(glm(z ~ -1 + offset(q[, 1]) + h0w + h1w, family = "binomial",
      subset = site == 0))

  q1 <- q + c((epsilon[1] * h0w + epsilon[2] * h1w),
              I(site == 0) * epsilon[1] / (g0w * ps0),
              I(site == 0) * epsilon[2] / (g1w * ps0))

  tmleest <- mean(plogis(q1[, 3][site == 0])) - mean(plogis(q1[, 2][site == 0]))

  eic <- (((a * h1w) - ((1 - a) * h0w)) * (z - plogis(q[, 1]))) +
    ((I(site == 0) / ps0) * ((plogis(q1[, 3]) - plogis(q1[, 2])) - tmleest))

  results = list("est" = tmleest,
                 "var" = var(eic) / n.dat,
                 "eic" = eic)

  return(results)
}
