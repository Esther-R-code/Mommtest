#' @title PSE with bootstrap
#' @param data.boot the data set
#' @param Para the parameter that will trated as nominal variables. Seeting as "~W+Q+S". The detail please see the R function 'clm'
#' @param n_para the number of pararmeters that used in the model of outcome
#' @param n.cata the number of categories of the outcome
#' @param ITE the total effect
#' @param theta.hat the estimate of theta (all of the parameters of the models)
#' @param sig.hat the eastimate of sigma (all of the variance of the parameters)
#' @param Vcov.matrix the variance-covariance matrix of theta
#' @param intv number of interventions. Default is 4.
#' @param intval the values of exposure used for intervention. Default is c(0,1)
#' @param confounders the values of confounders, def:c()
#' @param a the value of outcome
#' @param nb: number of bootstrapping
#' @param n.core: number of cores that will be used
#'
#' @export PSEa_wiboot_c
#'
PSEa_wiboot_c=function(data_boot, Para, n.para, n.cate, ITE, theta.hat, sig.hat, Vcov.matrix, intv, intval, confounders, a, nb, n.core){
  n.cate1=n.cate-1
  TEa<-NULL
  TEa1<-NULL
  for(con.count in 1:(n.para-2))
  {TEa<-c(TEa,ITE[(con.count-1)*(n.cate)+(a+1)])
  TEa1<-c(TEa1,ITE[(con.count-1)*(n.cate)+(a)])}

  rTEa<-total_a_r(TEa,TEa1,intval,confounders)

  rd.TEa<-rTEa[2]-rTEa[1]
  rr.TEa<-rTEa[2]/rTEa[1]

  ###########################################
  repeat{
    RNGkind("L'Ecuyer-CMRG")
    create_boot_risk<-mclapply(data_boot,"create_boot_risk_one", confounders=confounders, Para=Para, n.cate=n.cate, a=a, intv=intv, intval=intval,mc.cores=n.core)
    #create_boot_risk<-lapply(data_boot,"create_boot_risk_one", confounders=confounders, Para=Para, n.cate=n.cate, a=a, intv=intv, intval=intval)

    nb_r<-nb*1.2
    risk.boot<-matrix(unlist(create_boot_risk),nrow=nb_r,byrow=T)
    #print("dim of risk.boot=")
    #print(dim(risk.boot))
    if(any(risk.boot=="NaN")){
      ind<-which(risk.boot[,1]=="NaN"|risk.boot[,2]=="NaN"|risk.boot[,3]=="NaN"|risk.boot[,4]=="NaN"|risk.boot[,5]=="NaN"|	risk.boot[,6]=="NaN"|risk.boot[,7]=="NaN"|risk.boot[,8]=="NaN")
      if(length(ind)>0){
        risk.boot<-risk.boot[-ind,]}
    }
    if(any(risk.boot<0)){
      ind<-which(risk.boot[,2]<0|risk.boot[,4]<0|risk.boot[,6]<0|risk.boot[,8]<0)
      if(length(ind)>0){
        risk.boot<-risk.boot[-ind,]}
    }
    #print("dim of risk.boot=")
    #print(dim(risk.boot))
    if(dim(risk.boot)[1]>nb){break}
  }
  risk.boot<-risk.boot[1:nb,]
  total<- dim(risk.boot)[1]
  b.RD1<- c(quantile(risk.boot[,1], probs=c(0.025, 0.975),na.rm=T))
  b.RD2<- c(quantile(risk.boot[,3], probs=c(0.025, 0.975),na.rm=T))
  b.RD3<- c(quantile(risk.boot[,5], probs=c(0.025, 0.975),na.rm=T))
  b.RR1<- c(quantile(risk.boot[,2], probs=c(0.025, 0.975),na.rm=T))
  b.RR2<- c(quantile(risk.boot[,4], probs=c(0.025, 0.975),na.rm=T))
  b.RR3<- c(quantile(risk.boot[,6], probs=c(0.025, 0.975),na.rm=T))
  #null hypothesis is 0 no effect in path
  valRD1<- sum(risk.boot[,1]>0)/total
  pvalRD1<- min(valRD1,(1-valRD1))*2
  valRD2<- sum(risk.boot[,3]>0)/total
  pvalRD2<- min(valRD2,(1-valRD2))*2
  valRD3<- sum(risk.boot[,5]>0)/total
  pvalRD3<- min(valRD3,(1-valRD3))*2
  #null hypothesis is log(1) no difference in log(relative risk). This is done in the log scale since it is more normally distributed
  valRR1<- sum(log(risk.boot[,2])>0)/total
  pvalRR1<- min(valRR1,(1-valRR1))*2
  if(length(which(log(risk.boot[,2])=="NaN"))!=0){print("risk.boot[,2]");print(which(log(risk.boot[,2])=="NaN"))}
  valRR2<- sum(log(risk.boot[,4])>0)/total
  pvalRR2<- min(valRR2,(1-valRR2))*2
  if(length(which(log(risk.boot[,4])=="NaN"))!=0){print("risk.boot[,4]");print(which(log(risk.boot[,4])=="NaN"))}
  valRR3<- sum(log(risk.boot[,6])>0)/total
  pvalRR3<- min(valRR3,(1-valRR3))*2
  if(length(which(log(risk.boot[,6])=="NaN"))!=0){print("risk.boot[,6]");print(which(log(risk.boot[,6])=="NaN"))}
  bdnp<- c("lower(a)","upper(a)","pv(a)","lower(b)","upper(b)","pv(b)")
  #}
  ####################################################

  #if(intv==4){
  w0000 <- rep(intval[1],4)
  w1000 <- c(intval[2],rep(intval[1],3))
  w1100 <- c(rep(intval[2],2),rep(intval[1],2))
  w1110 <- c(rep(intval[2],3),intval[1])
  w1111 <- rep(intval[2],4)
  p0000<- rho(theta.hat,sig.hat,n.para,Name=c("W","Q","S"),n.cate1,a,w0000,confounders)
  p1000<- rho(theta.hat,sig.hat,n.para,Name=c("W","Q","S"),n.cate1,a,w1000,confounders)
  p1100<- rho(theta.hat,sig.hat,n.para,Name=c("W","Q","S"),n.cate1,a,w1100,confounders)
  p1110<- rho(theta.hat,sig.hat,n.para,Name=c("W","Q","S"),n.cate1,a,w1110,confounders)
  p1111<- rho(theta.hat,sig.hat,n.para,Name=c("W","Q","S"),n.cate1,a,w1111,confounders)

  rho_values<- c(p0000[1],p1000[1],p1100[1],p1110[1],p1111[1],rd.TEa,rr.TEa)
  names(rho_values)<- c("p0000","p1000","p1100","p1110","p1111","total RD","total RR")
  RD1<- rd_ord(p1000,p0000,Vcov.matrix)
  RD2<- rd_ord(p1100,p1000,Vcov.matrix)
  RD3<- rd_ord(p1110,p1100,Vcov.matrix)
  RD4<- rd_ord(p1111,p1110,Vcov.matrix)
  RR1<- rr_ord(p1000,p0000,Vcov.matrix)
  RR2<- rr_ord(p1100,p1000,Vcov.matrix)
  RR3<- rr_ord(p1110,p1100,Vcov.matrix)
  RR4<- rr_ord(p1111,p1110,Vcov.matrix)

  #  if(nb>0){
  b.RD4<- c(quantile(risk.boot[,7], probs=c(0.025, 0.975),na.rm=T))
  b.RR4<- c(quantile(risk.boot[,8], probs=c(0.025, 0.975),na.rm=T))
  valRD4<- sum(risk.boot[,7]>0)/total
  pvalRD4<- min(valRD4,(1-valRD4))*2
  valRR4<- sum(log(risk.boot[,8])>0)/total
  pvalRR4<- min(valRR4,(1-valRR4))*2
  if(length(which(log(risk.boot[,8])=="NaN"))!=0){print("risk.boot[,8]");print(which(log(risk.boot[,8])=="NaN"))}
  pse_values<- c(RD1,b.RD1,pvalRD1,RR1,b.RR1,pvalRR1,RD2,b.RD2,pvalRD2,RR2,b.RR2,pvalRR2,
                 RD3,b.RD3,pvalRD3,RR3,b.RR3,pvalRR3,RD4,b.RD4,pvalRD4,RR4,b.RR4,pvalRR4)
  #    }else{
  #    pse_values<- c(RD1,RR1,RD2,RR2,RD3,RR3,RD4,RR4)
  #  }
  names(pse_values)<- c("RD W>Y",bdnp,"RR W>Y",bdnp,"RD W>S>Y",bdnp,"RR W>S>Y",bdnp,"RD W>Q>Y",bdnp,"RR W>Q>Y",bdnp,"RD W>Q>S>Y",bdnp,"RR W>Q>S>Y",bdnp)
  #}
  ######################################################################
  return(c(rho_values,pse_values))
}
