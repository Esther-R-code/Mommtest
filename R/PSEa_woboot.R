#' @title PSE without bootstrap
#' @param data.boot the data set
#' @param Para the parameter that will trated as nominal variables. Seeting as "~W+Q+S". The detail please see the R function 'clm'
#' @param n_para the number of pararmeters that used in the model of outcome
#' @param n_cate the number of categories of the outcome
#' @param ITE the total effect
#' @param theta.hat the estimate of theta (all of the parameters of the models)
#' @param sig.hat the eastimate of sigma (all of the variance of the parameters)
#' @param Vcov.matrix the variance-covariance matrix of theta
#' @param intv number of interventions. Default is 4.
#' @param intval the values of exposure used for intervention. Default is c(0,1)
#' @param confounders the values of confounders, def:c()
#' @param a the value of outcome
#' @export PSEa_woboot

PSEa_woboot=function(data_boot, Para, n_para, n_cate, ITE, theta.hat, sig.hat, Vcov.matrix, intv, intval, confounders, a){
  n_cate1=n_cate-1
  TEa<-NULL
  TEa1<-NULL
  for(con.count in 1:(n_para-2))
  {TEa<-c(TEa,ITE[(con.count-1)*(n_cate)+(a+1)])
  TEa1<-c(TEa1,ITE[(con.count-1)*(n_cate)+(a)])}

  rTEa<-total_a_r(TEa,TEa1,intval,confounders)

  rd.TEa<-rTEa[2]-rTEa[1]
  rr.TEa<-rTEa[2]/rTEa[1]
  bdnp<- c("lower(a)","upper(a)","pv(a)")
  ####################################################

  w0000<- rep(intval[1],4)
  w1000<- c(intval[2],rep(intval[1],3))
  w1100<- c(rep(intval[2],2),rep(intval[1],2))
  w1110<- c(rep(intval[2],3),intval[1])
  w1111<- rep(intval[2],4)
  p0000<- rho(theta.hat,sig.hat,Name=c("W","Q","S"),n_cate1,a,w0000,confounders)
  p1000<- rho(theta.hat,sig.hat,Name=c("W","Q","S"),n_cate1,a,w1000,confounders)
  p1100<- rho(theta.hat,sig.hat,Name=c("W","Q","S"),n_cate1,a,w1100,confounders)
  p1110<- rho(theta.hat,sig.hat,Name=c("W","Q","S"),n_cate1,a,w1110,confounders)
  p1111<- rho(theta.hat,sig.hat,Name=c("W","Q","S"),n_cate1,a,w1111,confounders)

  rho_values<- c(p0000[1],p1000[1],p1100[1],p1110[1],p1111[1],rd.TEa,rr.TEa)
  names(rho_values)<- c("p0000","p1000","p1100","p1110","p1111","total RD","total RR")
  #null hypothesis is 0 no effect in path
  RD1<- rd_ord(p1000,p0000,Vcov.matrix)
  RD2<- rd_ord(p1100,p1000,Vcov.matrix)
  RD3<- rd_ord(p1110,p1100,Vcov.matrix)
  RD4<- rd_ord(p1111,p1110,Vcov.matrix)
  #null hypothesis is log(1) no difference in log(relative risk). This is done in the log scale since it is more normally distributed
  RR1<- rr_ord(p1000,p0000,Vcov.matrix)
  RR2<- rr_ord(p1100,p1000,Vcov.matrix)
  RR3<- rr_ord(p1110,p1100,Vcov.matrix)
  RR4<- rr_ord(p1111,p1110,Vcov.matrix)

  pse_values<- c(RD1,RR1,RD2,RR2,RD3,RR3,RD4,RR4)
  names(pse_values)<- c("RD W>Y",bdnp,"RR W>Y",bdnp,"RD W>S>Y",bdnp,"RR W>S>Y",bdnp,"RD W>Q>Y",bdnp,"RR W>Q>Y",bdnp,"RD W>Q>S>Y",bdnp,"RR W>Q>S>Y",bdnp)
  ######################################################################
  return(c(rho_values,pse_values))
}
