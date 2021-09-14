#This function is for sensitvity analysis
#' @param eta_msdatap the measured confounder without esp
mean_PSE_U<- function(eta_msdatap, msdata, confounders_ms, n_cate, n_cate1, sc, n_coaU, n_paraU, intv, intval, nb, j, etasd, n_core){
  ndata<- cbind(msdata, eta_msdatap)
  if(dim(msdata)[1]==length(eta_msdatap)){
    n_sample<-length(eta_msdatap)
  }else{cat("Error")}
  wc<- seq(1,n.cate,by=1); wc<- wc[-sc]
  wc_count<-length(wc)
  Outcome<-paste0("Outcome=",wc)
  path<- c("W>Y","W>S>Y","W>Q>Y","W>Q>S>Y")
  #repeat 100 times
  i<- 1
  PSE_100_RD<- 0
  PSE_100_RR<- 0
  repeat{
    if(i>101){break}
    ndata$Ums<- ndata$eta_msdatap+rnorm(n_sample,mean=0,sd=etasd)
    confoundersU<- c(confounders_ms, median(ndata$Ums))
    PSE_ms<-Ord_mediation_analysis_pal_sc(Indata=ndata[,-which(colnames(ndata)=="eta_msdatap")], n.cate=n_cate, n.coa=n_coaU, n.para=n_paraU, n.cate1=n_cate1, confounders=confoundersU, intv=intv, intval=intval, nb=nb, j=j, w_c=wc, n.core=n_core)
    #matrix of PSE under RD
    PSE_RD_vector<- NULL
    for(pathcount in 1:4){
      for(catecount in 1:wc_count){PSE_RD_vector<-c(PSE_RD_vector,PSE_ms[paste0("RD ",path[pathcount]), paste0("Outcome=",wc[catecount])])}
    }
    PSE_RD<- matrix(PSE_RD_vector, ncol = wc_count, nrow = 4, byrow = T)
    PSE_100_RD<- PSE_100_RD+PSE_RD
    #matrix of PSE under RD

    #matrix of PSE under RR
    PSE_RR_vector<- NULL
    for(pathcount in 1:4){
      for(catecount in 1:wc_count){PSE_RR_vector<-c(PSE_RR_vector,PSE_ms[paste0("RR ",path[pathcount]), paste0("Outcome=",wc[catecount])])}
    }
    PSE_RR<- matrix(PSE_RR_vector, ncol = wc_count, nrow = 4, byrow = T)
    PSE_100_RR<- PSE_100_RR+log(PSE_RR)
    #matrix of PSE under RR
    i<- i+1
  }
  #repeat 100 times
  PSE_mean_RD<- PSE_100_RD/100
  PSE_mean_RR<- PSE_100_RR/100
  colnames(PSE_mean_RD)<-Outcome; colnames(PSE_mean_RR)<-Outcome
  return(list(PSE_mean_RD,PSE_mean_RR))
}
