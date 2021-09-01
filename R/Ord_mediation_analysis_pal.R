#' Apply mediation analysis for non-rare orinal/multinomial outcome with two gaussian mediators
#' @title Ord_mediation_analysis
#' multinomial outcome mediation analysis with two ordered mediators
#' @param Indata: Input data with outcome, exposure, mediators, and covariates as columns. The number of rows equals to the number of people in the data you used.
#' @param n_cate: the number of categories of the outcome
#' @param confounders: the values of confounder
#' @param intv: number of interventions. Default is 4 and only 4.
#' @param intval: the values of exposure used for intervention. Default is c(0,1).
#' @param nb: Number of bootstrapping. Default is 0 (i.e. no bootstrapping applied).
#' @param n_core: number of cores that will be used
#' @export Ord_mediation_analysis_pal
#' @example
#' library(ordinal)
#' library(parallel)
#' simluate a data set
#' sample_size<- 500
#' exp_cont<- c(0,1)
#' med_cont<- c(0.23,0.23)
#' n_cate<- 3
#' deltaq<- c(0.1,0.5) #Q c(d0,dw)
#' alphas<- c(0.2,0.5,0.5) #S c(al0,alw,alq)
#' betay1<- c(-1.8, 1.0, 1.2, 1.3) #c(b01,bw1,bq1,bs1)
#' betay2<- c(0.8, 0.8, 1.1, 1.2) #c(b02,bw2,bq2,bs2)
#' betay<- rbind(betay1,betay2)
#' coa<- c(0,0.4)
#' j<-56 #sample(1,1:1e+09)
#' sim_data<-lapply(j,"sim_medi_cate_data",sample_size=sample_size,exp_cont=exp_cont,med_cont=med_cont,n_cate=n_cate,deltaq=deltaq,alphas=alphas,betay=betay,coa=coa)
#' seeting of mutiple mediation analysis
#' confounders<- median(sim_data[[1]][,5])
#' intv<- 4
#' intval<- c(0, 1)
#' nb<- 0
#' n_core<- 8
#' result<- Ord_mediation_analysis_pal(Indata=sim_data[[1]], n_cate, confounders, intv, intval, nb, n_core)

Ord_mediation_analysis_pal=function(Indata, n_cate, confounders, intv=4, intval=c(0,1), nb=0, n_core){
  #dt=as.data.frame(Indata)
  colnames(Indata)[1:4]<-c("Y","W","Q","S")
  #n_para<-dim(Indata)[2]
  #n_cate<-length(names(summary(Indata$Y)))
  n_cate1<-n_cate-1
  n.coa<-length(confounders)
  n_para<-4+n.coa

  Para<-"~W+Q+S"
  ParaTE<-"~W"
  if(n_para>4){
    n.con<-dim(Indata)[2]-4
    for(con.count in 1:n.con){
      colnames(Indata)[(4+con.count)]<-paste0("Con",con.count)
      Para<-paste0(Para,"+Con",con.count)
      ParaTE<-paste0(ParaTE,"+Con",con.count)
    }}
  out.reg<-clm(Y ~ 1, nominal=Para, data=Indata, link="probit")
  s.reg<-lm(S~.-Y, data=Indata)
  q.reg<-lm(Q~.-S-Y, data=Indata)
  beta.hat<-out.reg$coefficients
  alpha.hat<-s.reg$coefficients
  delta.hat<-q.reg$coefficients
  #sigbeta.hat<-summary(at.out.reg)$coefficients[,2]
  sigs.hat<-sqrt(mean((s.reg$residuals)^2))#MLE of sigmaq  sum[(s-s^)^2]/n
  sigq.hat<-sqrt(mean((q.reg$residuals)^2))#MLE of sigmas  sum[(q-q^)^2]/n
  theta.hat<- list(beta.hat, alpha.hat, delta.hat)
  sig.hat<- c(sigq.hat, sigs.hat)
  Vcov.matrix<-create_vcovmatrix(out.reg,s.reg,q.reg)

  TE<-clm(Y ~ 1-Q-S, nominal=ParaTE, data=Indata, link="probit")$coefficients#total effect
  ITE<-NULL
  for(para.count in 1:(n_para-2)){
    startat<-(para.count-1)*(n_cate1)
    ITE<-c(ITE,-Inf,TE[(startat+1):(startat+(n_cate1))])}

  PSE<-NULL
  coln_PSE<-NULL

  #create_boot_data#####
  if(nb>0){
    nb_r<- nb*1.2
    RNGkind("L'Ecuyer-CMRG")
    data.boot<-mclapply(1:nb_r, "create_data_boot_one", Indata=Indata, j=j, mc.cores=n_core)
  }

  ####fixed###########################################
  if(nb>0){
    for(cate.count in 1:(n_cate)){
      # for each Y=a=cate.count calculate the PSEs
      PSE_a<-PSEa_wiboot_c(data.boot, Para=Para, n_para=n_para, n_cate=n_cate, ITE=ITE, theta.hat=theta.hat, sig.hat=sig.hat, Vcov.matrix=Vcov.matrix, intv=intv, intval=intval, confounders=confounders, a=cate.count, nb=nb, n_core=n_core)
      PSE<-cbind(PSE,PSE_a)#;print(PSE)
      coln_PSE<-c(coln_PSE,paste0("Outcome=",cate.count))#;print(coln_PSE)
    }
  }else{
   for(cate.count in 1:(n_cate)){
   	PSE_a<-PSEa_woboot(Indata, Para=Para, n_para=n_para, n_cate=n_cate, ITE=ITE, theta.hat=theta.hat, sig.hat=sig.hat, Vcov.matrix=Vcov.matrix, intv=intv, intval=intval, confounders=confounders, a=cate.count, n_core=n_core)
  	PSE<-cbind(PSE,PSE_a)#;print(PSE)
  	coln_PSE<-c(coln_PSE,paste0("Outcome=",cate.count))#;print(coln_PSE)
  	}
  }
  ####fixed###################################
  colnames(PSE)<-coln_PSE
  return(PSE)
}
