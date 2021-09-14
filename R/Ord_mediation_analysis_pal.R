#' Apply mediation analysis for non-rare orinal/multinomial outcome with two gaussian mediators.
#' Before using this function, please make sure the R packages "ordinal" and "parallel" have been install and libraries.
#' @title Ord_mediation_analysis
#' @param Indata: Input data with outcome, exposure, mediators, and covariates as columns. The number of rows equals to the number of people in the data you used.
#' @param n_cate: the number of categories of the outcome
#' @param confounders: the values of confounder
#' @param intv: number of interventions. Default is 4 and only 4.
#' @param intval: the values of exposure used for intervention. Default is c(0,1).
#' @param nb: Number of bootstrapping. Default is 0 (i.e. no bootstrapping applied).
#' @param n_core: number of cores that will be used in parallel computing. default=1
#' @export Ord_mediation_analysis_pal
#' @examples
#' library(ordinal)
#' library(parallel)
#' #simluate a data set
#' sample_size<- 500
#' exp_cont<- c(0,1)
#' med_cont<- c(0.23,0.23)
#' n_cate<- 3
#' deltaq<- c(0.1,0.5) #Q c(d0,dw)
#' alphas<- c(0.2,0.5,0.5) #S c(al0,alw,alq)
#' betay1<- c(-1.8, 1.0, 1.2, 1.3, 0.05) #c(b01,bw1,bq1,bs1,bc1)
#' betay2<- c(0.8, 0.8, 1.1, 1.2, 0.03) #c(b02,bw2,bq2,bs2,bc2)
#' betay<- rbind(betay1,betay2)
#' coa<- c(0,0.4)
#' j<- 56 #sample(1,1:1e+09)
#' sim_data<- sim_medi_cate_data(sample_size=sample_size, exp_cont=exp_cont, med_cont=med_cont, n_cate=n_cate, deltaq=deltaq, alphas=alphas, betay=betay, coa=coa,j)
#' #seeting of mutiple mediation analysis
#' confounders<- median(sim_data[,5])
#' intv<- 4
#' intval<- c(0, 1)
#' nb<- 0
#' n_core<- 2
#' result<- Ord_mediation_analysis_pal(Indata=sim_data, n_cate, confounders, intv, intval, nb, n_core)
#'
#' #if want to calculate the confidence interval by using bootstrapping
#' nb<- 50
#' result_boot<- Ord_mediation_analysis_pal(Indata=sim_data, n_cate, confounders, intv, intval, nb, n_core)

#' @return PSE: a table with the following values in the row; columns store the values for different categories of the outcome
#' @return p0000: the probability of the counterfactual outcome of the outcome under the intervention=(0,0,0,0)
#' @return p1000: the probability of the counterfactual outcome of the outcome under the intervention=(1,0,0,0)
#' @return p1100: the probability of the counterfactual outcome of the outcome under the intervention=(1,1,0,0)
#' @return p1110: the probability of the counterfactual outcome of the outcome under the intervention=(1,1,1,0)
#' @return p1111: the probability of the counterfactual outcome of the outcome under the intervention=(1,1,1,1)
#' @return total RD: the total effect under the risk difference scale
#' @return total RR: the total effect under the risk ratio scale
#' @return RD W>Y: the PSE of path from W to Y under the risk difference scale
#' @return lower(RDWY): the lower bound of the confidence interval the PSE of path from W to Y under the risk difference scale
#' @return upper(RDWY): the upper bound of the confidence interval the PSE of path from W to Y under the risk difference scale
#' @return pv(RDWY): the p-value of the PSE of path from W to Y under the risk difference scale
#' @return RR W>Y: the PSE of path from W to Y under the risk ratio scale
#' @return lower(RRWY): the lower bound of the confidence interval the PSE of path from W to Y under the risk ratio scale
#' @return upper(RRWY): the upper bound of the confidence interval the PSE of path from W to Y under the risk ratio scale
#' @return pv(RRWY): the p-value of the PSE of path from W to Y under the risk ratio scale
#' @return Note: other return values with the similar names are for the path from W through Q to Y  or through  S, or through both
Ord_mediation_analysis_pal=function(Indata, n_cate, confounders, intv=4, intval=c(0,1), nb=0, n_core=1){
  library(ordinal)
  library(parallel)
  colnames(Indata)[1:4]<-c("Y","W","Q","S")
  if(n_cate!=length(names(summary(Indata$Y)))){stop("please make sure 'n_cate' is equal to the number of categories of the outcome!")}
  n_cate1<-n_cate-1
  n_para<-4+length(confounders)
  if(n_para!=dim(Indata)[2]){stop("please make sure the length of 'confounders' equal to the number of covariates of the model of outcome!")}

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
    databoot<- mclapply(1:nb_r, "create_data_boot_one", Indata=Indata, mc.cores=n_core, mc.cleanup = TRUE)
  }

  if(nb>0){
    for(cate.count in 1:(n_cate)){
      # for each Y=a=cate.count calculate the PSEs
      PSE_a<-PSEa_wiboot_c(data_boot=databoot, Para=Para, n_para=n_para, n_cate=n_cate, ITE=ITE, theta.hat=theta.hat, sig.hat=sig.hat, Vcov.matrix=Vcov.matrix, intv=intv, intval=intval, confounders=confounders, a=cate.count, nb=nb, n_core=n_core)
      PSE<-cbind(PSE,PSE_a)#;print(PSE)
      coln_PSE<-c(coln_PSE,paste0("Outcome=",cate.count))#;print(coln_PSE)
    }
  }else{
   for(cate.count in 1:(n_cate)){
   	PSE_a<-PSEa_woboot(Indata, Para=Para, n_para=n_para, n_cate=n_cate, ITE=ITE, theta.hat=theta.hat, sig.hat=sig.hat, Vcov.matrix=Vcov.matrix, intv=intv, intval=intval, confounders=confounders, a=cate.count)
  	PSE<-cbind(PSE,PSE_a)#;print(PSE)
  	coln_PSE<-c(coln_PSE,paste0("Outcome=",cate.count))#;print(coln_PSE)
  	}
  }

  colnames(PSE)<-coln_PSE
  return(PSE)
}
