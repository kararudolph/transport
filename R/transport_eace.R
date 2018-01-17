#' Average causal effect of the exposure on the outcome, ignoring the instrument
#'
#' @examples
#'
#' # TBD.
#'
#' @references TBD.
#'
#' @importFrom stats coef glm plogis predict var
#' @export
transport_eace <-
  function(a,
           z,
           y,
           site ,
           w,
           nsitemodel ,
           nzmodel,
           noutmodel) {
  datw <- w
  n.dat <- nrow(w)

  # Calculate components of clever covariate
  glm_cps = glm(formula = nsitemodel, data = data.frame(cbind(site = site, datw)),
                family = "binomial")
  cps <- predict(glm_cps, type = "response")

  glm_cpz = glm(formula = nzmodel, data = data.frame(cbind(a = a, z = z, datw)),
      family = "binomial")
  cpz <- predict(glm_cpz, type = "response")

  # Calculate clever covariate.
  g0w <- ((1 - cpz) * cps) / (1 - cps)
  g1w <- (cpz * cps) / (1 - cps)
  h0w <- ((1 - z) * I(site == 1)) / g0w
  h1w <- (z * I(site == 1)) / g1w

  ymodel <- glm(formula = noutmodel, family = "binomial",
      data = data.frame(cbind(datw, a = a, z = z, site = site , y = y)),
      subset = site == 1)

  data_new0 <- data_new1 <- datw
  data_new0$z <- 0

  data_new1$z <- 1

  # Initial prediction.
  q <- cbind(predict(ymodel, type = "link",
                     newdata = data.frame(cbind(datw, a = a, z = z))),
      predict(ymodel, type = "link", newdata = data_new0),
      predict(ymodel, type = "link", newdata = data_new1))

  epsilon <- coef(glm(y ~ -1 + offset(q[, 1]) + h0w + h1w, family = "binomial",
                      subset = site == 1))
  # Update initial prediction.
  q1 <- q + c((epsilon[1] * h0w + epsilon[2] * h1w),
              epsilon[1] / g0w, epsilon[2] / g1w)

  # Get efficient influence curve values for everyone
  tmleest <- mean(plogis(q1[, 3][site == 0])) - mean(plogis(q1[, 2][site == 0]))
  ps0 <- mean(I(site == 0))
  eic <- (((z * h1w / ps0) - ((1 - z) * h0w / ps0)) * (y - plogis(q[, 1]))) +
    (I(site == 0) / ps0 * plogis(q1[, 3])) - (I(site == 0) / ps0 * plogis(q1[, 2]))
    - (tmleest / ps0)

  results = list("est" = tmleest,
                 "var" = var(eic) / n.dat,
                 "eic" = eic[site == 0])
  return(results)
}