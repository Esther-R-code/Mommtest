#' this function calculate the PSEs under risk ratio scale and 95% confidence interval
#' @param pa rho under one intervention
#' @param pb rho under another intervention
#' @param V the variance-covariance matrix of theta
#' @export rr_ord
rr_ord=function(pa,pb,V){
  path=pa[1]/pb[1]
  der=matrix((pa[-1]/pa[1])-(pb[-1]/pb[1]),nrow=1)#(prho/ptheta)/rho
  std.error=sqrt(der%*%V%*%t(der))
  pv=2*pnorm(abs(log(path)/std.error),lower.tail=F)
  return(c(path, exp(log(path)-1.96*std.error), exp(log(path)+1.96*std.error), pv))
}
