#' This function is for sensitivity analysis
#' It will calculate the measured confounder without the disturbance(epsilon)
#' reference: 'Sensitivity analysis.pdf' in 'ref' file folder
#' @title measured confounder
#' @param eta eta
#' @param ratio eta_tilde
#' @param changevariable the names of variables that used to create a measured confounder; ex: c("Q", "W")
#' @param msdatap the data only contain the exposure, mediators and outcome; ex: cbind(Y, Q, S, W)
#' @return the measured confounder without the disturbance(epsilon)
etamsdatap_de<- function(eta, ratio, changevariable, msdatap){
  #n_smaple<- dim(msdatap)[1]
  etavector<-seq(0,0,length.out = 4)
  etavector[which(colnames(msdatap)==changevariable[1])]<- eta
  etavector[which(colnames(msdatap)==changevariable[2])]<- eta*ratio
  eta_msdatap_de_1<-msdatap%*%etavector
  return(eta_msdatap_de_1)
}
