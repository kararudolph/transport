#' Transported complier average causal effect (CACE) (aka treatment effect on treated)
#'
#' Details to be added.
#'
#' <Add parameter info>
#'
#' <Add return info>
#'
#' @examples
#'
#' # TBD.
#'
#' @references
#'
#' Imbens and Rubin, 1997
#'
#' @importFrom stats coef glm plogis predict var cov
#' @export
transport_cace <-
  function(ca,
           cz,
           cy,
           csite,
           cw,
           czmodel,
           csitemodel,
           coutmodel,
           cq2model) {

  datw <- cw
  n.dat <- nrow(datw)
  ps0 <- mean(I(csite == 0))

  camodel <- "a ~ 1"

  notransportate <- notransport(a = ca, z = cz, site = csite, w = cw,
                                ntamodel = camodel , ntzmodel = czmodel)

  ittate <- transport_ittate(a = ca, z = cz, y = cy, site = csite, w = cw,
                             aamodel = camodel, asitemodel = csitemodel,
                             #azmodel = czmodel,
                             s_awz_model = czmodel,
                             aoutmodel = coutmodel,
                             aq2model = cq2model)

  cate <- ittate$est / notransportate$est
  varcate <- (ittate$est^2 / notransportate$est^2) *
    (((ittate$var * n.dat) / ittate$est^2) -
       ((2 * cov(cbind(ittate$eic , notransportate$eic))[1 , 2]) /
          (ittate$est * notransportate$est)) + ((notransportate$var * n.dat) /
                                                  notransportate$est^2))
  eic <- (ittate$eic / notransportate$est) -
    (ittate$est / (notransportate$est^2)) * notransportate$eic

  results = list("est" = cate,
                 "var" = var(eic) / n.dat,
                 "eic" = eic)
  return(results)
}