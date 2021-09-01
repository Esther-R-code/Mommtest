#' this function calculate the PSEs under risk difference scale and 95% confidence interval
#' @param pa rho under one intervention
#' @param pb rho under another intervention
#' @param V the variance-covariance matrix of theta
#' @export rd_ord
rd_ord=function(pa,pb,V){
  path=pa[1]-pb[1]
  der=matrix(pa[-1]-pb[-1], nrow=1)
  std.error=sqrt(der%*%V%*%t(der))
  pv=2*pnorm(abs(path/std.error),lower.tail=F)
  return(c(path, path-1.96*std.error, path+1.96*std.error, pv))
}
