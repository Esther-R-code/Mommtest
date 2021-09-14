#' total_a_r the risk of total effet when the value of outcome is a
#' @title total_a_r
#' @param TEa total effect under outcome=a
#' @param TEa1 total effect under outcome=a-1
#' @param confounders values of Confounders
#' @param intval the values of exposure used for intervention. Default is c(0,1).
#' @export total_a_r

total_a_r=function(TEa,TEa1,intval=c(0,1),confounders){
  wi0<-intval[1]
  wi1<-intval[2]
  total.a.r<-c(NA,NA)
  if(sum(TEa1)==-Inf){
  total.a.r[2]<-pnorm(sum(TEa*c(1,wi1,confounders)))
  total.a.r[1]<-pnorm(sum(TEa*c(1,wi0,confounders)))
  }else if(sum(TEa,na.rm=T)==-Inf){
    total.a.r[2]<-1-pnorm(sum(TEa1*c(1,wi1,confounders)))
    total.a.r[1]<-1-pnorm(sum(TEa1*c(1,wi0,confounders)))
  }
  else{gtheta1<-pnorm(sum(TEa1*c(1,wi1,confounders)))
  gtheta0<-pnorm(sum(TEa1*c(1,wi0,confounders)))
  total.a.r[2]<-pnorm(sum(TEa*c(1,wi1,confounders)))-gtheta1
  total.a.r[1]<-pnorm(sum(TEa*c(1,wi0,confounders)))-gtheta0} #pnorm: cdf of normal dist.
  return(total.a.r)
}
