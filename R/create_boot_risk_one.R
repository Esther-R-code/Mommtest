#' @title the risk create by bootstrapping
#' @param data_boot the data set
#' @param confounders the values of confounders, def:c()
#' @param Para the parameter that will trated as nominal variables. Seeting as "~W+Q+S". The detail please see the R function 'clm'
#' @param n_cate the number of categories of the outcome
#' @param a the value of outcome
#' @param intv number of interventions. Default is 4.
#' @param intval the values of exposure used for intervention. Default is c(0,1).
#' @export create_boot_risk_one
create_boot_risk_one=function(data_boot, confounders=c(), Para, n_cate, a, intv=4, intval=c(0,1)){
    n_cate1<- n_cate-1
    out.reg.bt<- clm(Y ~ 1, nominal=Para, data=data_boot, link="probit")
    s.reg.bt<- lm(S~.-Y, data=data_boot)
    q.reg.bt<- lm(Q~.-Y-S, data=data_boot)
    beta.hat.bt<- out.reg.bt$coefficients
    alpha.hat.bt<- s.reg.bt$coefficients
    delta.hat.bt<- q.reg.bt$coefficients
    ss.hat.bt<- sqrt(mean((s.reg.bt$residuals)^2))
    sq.hat.bt<- sqrt(mean((q.reg.bt$residuals)^2))
    theta.hat.bt<- list(beta.hat.bt, alpha.hat.bt, delta.hat.bt)
    sig.hat.bt<- c(sq.hat.bt, ss.hat.bt)
    if(intv==3){
      w000<- rep(intval[1],3)
      w100<- c(intval[2],intval[1],intval[1])
      w110<- c(intval[2],intval[2],intval[1])
      w111<- rep(intval[2],3)
      p000<- rho(theta.hat.bt,sig.hat.bt,Name=c("W","Q","S"),n_cate1,a,w000,confounders)
      p100<- rho(theta.hat.bt,sig.hat.bt,Name=c("W","Q","S"),n_cate1,a,w100,confounders) #first part of difference
      p110<- rho(theta.hat.bt,sig.hat.bt,Name=c("W","Q","S"),n_cate1,a,w110,confounders) #first part of difference
      p111<- rho(theta.hat.bt,sig.hat.bt,Name=c("W","Q","S"),n_cate1,a,w111,confounders) #first part of difference
      if(p000[1]==0|p100[1]==0|p110[1]==0|is.nan(p000[1])|is.nan(p100[1])|is.nan(p110[1])|is.nan(p111[1])|p000[1]<0|p100[1]<0|p110[1]<0|p111[1]<0){
       fail<- fail+1
     }else{
        #Stoppoint<- Stoppoint+1
        RD1<- p100[1]-p000[1]
        RD2<- p110[1]-p100[1]
        RD3<- p111[1]-p110[1]
        RR1<- p100[1]/p000[1]
        RR2<- p110[1]/p100[1]
        RR3<- p111[1]/p110[1]
        risk_boot<- c(RD1,RR1,RD2,RR2,RD3,RR3,NA,NA)}
    }else if(intv==4){
      w0000<- rep(intval[1],4)
      w1000<- c(intval[2],rep(intval[1],3))
      w1100<- c(rep(intval[2],2),rep(intval[1],2))
      w1110<- c(rep(intval[2],3),intval[1])
      w1111<- rep(intval[2],4)
      p0000<- rho(theta.hat.bt,sig.hat.bt,Name=c("W","Q","S"),n_cate1,a,w0000,confounders)
      p1000<- rho(theta.hat.bt,sig.hat.bt,Name=c("W","Q","S"),n_cate1,a,w1000,confounders) #first part of difference
      p1100<- rho(theta.hat.bt,sig.hat.bt,Name=c("W","Q","S"),n_cate1,a,w1100,confounders) #first part of difference
      p1110<- rho(theta.hat.bt,sig.hat.bt,Name=c("W","Q","S"),n_cate1,a,w1110,confounders) #first part of difference
      p1111<- rho(theta.hat.bt,sig.hat.bt,Name=c("W","Q","S"),n_cate1,a,w1111,confounders) #first part of difference
        RD1<- p1000[1]-p0000[1]
        RD2<- p1100[1]-p1000[1]
        RD3<- p1110[1]-p1100[1]
        RD4<- p1111[1]-p1110[1]
        RR1<- p1000[1]/p0000[1]
        RR2<- p1100[1]/p1000[1]
        RR3<- p1110[1]/p1100[1]
        RR4<- p1111[1]/p1110[1]
        risk_boot<- c(RD1,RR1,RD2,RR2,RD3,RR3,RD4,RR4)#}
    }
  return(risk_boot)
}
