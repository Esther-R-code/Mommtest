#' Apply mediation analysis for non-rare orinal/multinomial outcome with two gaussian mediators
#' for plot of bias
#' @title Ord_mediation_analysis for sensitivity analysis
#' @param Indata: Input data with outcome, exposure, mediators, and covariates as columns. The number of rows equals to the number of people in the data you used.
#' @param n_cate: the number of categories of the outcome
#' @param confounders: the values of confounder
#' @param intv: Number of intervention, only 4.
#' @param intval: the value of exposure used in Intervention , Default is c(0,1).
#' @param nb: Number of bootstrapping. Default is 0 (no bootstrapping applied).
#' @param n_core: number of cores that will be used in parallel computing. default=1
#' @param w_c<- cate of Y that user have instercting in. If all in let wc="all",wc will set to be seq(1,n.cate,by=1).
#' @export Ord_mediation_analysis_pal_sc
#'
#' @return please see the help docunment of function "Ord_mediation_analysis_pal"
#'
Ord_mediation_analysis_pal_sc=function(Indata, n_cate, confounders, intv=4, intval=c(0, 1), nb=0, w_c, n.core=1){
  colnames(Indata)[1:4]<-c("Y","W","Q","S")
  if(n_cate!=length(names(summary(Indata$Y)))){stop("please make sure 'n_cate' is equal to the number of categories of the outcome!")}
  n_cate1<-n_cate-1
  n_para<-4+length(confounders)
  if(n_para!=dim(Indata)[2]){stop("please make sure the length of 'confounders' equal to the number of covariates of the model of outcome!")}

  Para<-"~W+Q+S"
  ParaTE<-"~W"
  if(n.para>4){
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
  Vcov.matrix<-create_vcovmatrix(out.reg,s.reg,q.reg)

  #TE<-clm(Y ~ 1-Q-S, nominal=ParaTE, data=Indata, link="probit")$coefficients#total effect
  #ITE<-NULL
  #for(para.count in 1:(n.para-2)){
  #  startat<-(para.count-1)*(n_cate1)
  #  ITE<-c(ITE,-Inf,TE[(startat+1):(startat+(n_cate1))])}
  theta.hat<- list(beta.hat, alpha.hat, delta.hat)
  sig.hat<- c(sigq.hat, sigs.hat)

  PSE<-NULL
  coln_PSE<-NULL

  #create_boot_data#####
  if(nb>0){
    nb_r<- nb*1.2
    RNGkind("L'Ecuyer-CMRG")
    databoot<-mclapply(1:nb_r, "create_data_boot_one", Indata=Indata, mc.cores=n_core, mc.cleanup = TRUE)
  }

  if(w_c[1]=="all"){w_c<- seq(1, n_cate, by=1)}
  if(nb>0){
    for(cate.count in w_c){
      # for each Y=a=cate.count calculate the PSEs except the cate user want to skip
      PSE_a<-PSEa_wiboot_c(databoot, Para=Para, n.para=n.para, n_cate=n_cate, theta.hat=theta.hat, sig.hat=sig.hat, Vcov.matrix=Vcov.matrix, intv=intv, intval=intval, confounders=confounders, a=cate.count, nb=nb, n.core=n.core)
      PSE<-cbind(PSE,PSE_a)#;print(PSE)
      coln_PSE<-c(coln_PSE,paste0("Outcome=",cate.count))#;print(coln_PSE)
    }
  }else{
   for(cate.count in w_c){
   	PSE_a<-PSEa_woboot(Para=Para, n.para=n.para, n_cate=n_cate, theta.hat=theta.hat, sig.hat=sig.hat, Vcov.matrix=Vcov.matrix, intv=intv, intval=intval, confounders=confounders, a=cate.count)
  	PSE<-cbind(PSE,PSE_a)#;print(PSE)
  	coln_PSE<-c(coln_PSE,paste0("Outcome=",cate.count))#;print(coln_PSE)
  	}
  }

  colnames(PSE)<-coln_PSE
  return(PSE)
}
