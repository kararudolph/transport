#' Test of maximum mean discrepancy
#'
#' @examples
#'
#' # TBD.
#'
#' @references
#'
#' Luedtke et al. 2017
#'
#' @export
#'
#' @author Alex Luedtke
#'
#' @importFrom rARPACK eigs_sym
#' @importFrom stats rnorm pchisq quantile qchisq
mmd_test = function(R, S, D.R, D.S, sig.meth = 'eig', num.reps = 1e4,
                    return.cutoff = FALSE) {
  n = length(R)

  D.R.mat1 = matrix(rep(D.R,n),nrow=n)
  D.R.mat2 = matrix(rep(D.R,each=n),nrow=n)

  D.S.mat1 = matrix(rep(D.S,n),nrow=n)
  D.S.mat2 = matrix(rep(D.S,each=n),nrow=n)

  R.mat1 = matrix(rep(R,n),nrow=n)
  R.mat2 = matrix(rep(R,each=n),nrow=n)

  S.mat1 = matrix(rep(S,n),nrow=n)
  S.mat2 = matrix(rep(S,each=n),nrow=n)

  EE = ((2*(R.mat1-R.mat2)*(D.R.mat2-D.R.mat1) + 1 - (4*(R.mat1-R.mat2)^2-2)*D.R.mat1*D.R.mat2)*exp(-(R.mat1-R.mat2)^2)
        - ((2*(S.mat1-R.mat2)*(D.R.mat2-D.S.mat1) + 1 - (4*(S.mat1-R.mat2)^2-2)*D.S.mat1*D.R.mat2)*exp(-(S.mat1-R.mat2)^2))
        - ((2*(R.mat1-S.mat2)*(D.S.mat2-D.R.mat1) + 1 - (4*(R.mat1-S.mat2)^2-2)*D.R.mat1*D.S.mat2)*exp(-(R.mat1-S.mat2)^2))
        + (2*(S.mat1-S.mat2)*(D.S.mat2-D.S.mat1) + 1 - (4*(S.mat1-S.mat2)^2-2)*D.S.mat1*D.S.mat2)*exp(-(S.mat1-S.mat2)^2))

  # EE = exp(-(R.mat1-R.mat2)^2) - 2*exp(-(S.mat1-R.mat2)^2) + exp(-(S.mat1-S.mat2)^2)

  if(sig.meth=='eig'){
    line.means = rowMeans(EE)
    EE.ctrd = EE - matrix(rep(line.means,n),nrow=n) - matrix(rep(line.means,each=n),nrow=n) + matrix(rep(mean(line.means),n^2),nrow=n)
    num.eigs = min(200,n)
    # tmp = eigen(EE.ctrd)$values/n
    tmp = eigs_sym(EE.ctrd,num.eigs,which='LA')$values/n
    num.pos.eigs = num.eigs # sum(tmp>0)
    draws=c(matrix(rnorm(num.reps*num.pos.eigs)^2-1,nrow=num.reps,ncol=num.pos.eigs)%*%cbind(tmp[1:num.pos.eigs]))
  }

  # U-statistic
  diag(EE) = 0
  est = (rbind(rep(1/(n-1),n)) %*% EE %*% cbind(rep(1/n,n)))[1,1]
  # V-statistic
  # est = (rbind(rep(1/n,n)) %*% EE %*% cbind(rep(1/n,n)))[1,1]

  if(sig.meth=='eig'){
    pval = mean(draws>n*est)
  } else if(sig.meth=='var'){
    pval = pchisq(est/(2*var(D.R)/n)+1,df=1,lower.tail=FALSE)
  }

  return(if(!return.cutoff){
    c(est,pval)
  }else{
    c(est,pval,
      if(sig.meth=='eig'){
        quantile(draws,0.95)
      }else{
        2*var(D.R)*(qchisq(0.95,df=1)-1)})})
}
