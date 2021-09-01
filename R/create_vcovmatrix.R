#' create the variance-covariance matrix of theta
#' @title create_vcovmatrix
#' @param out.reg the output of the R function 'clm' i.e. the regression output of the outcome
#' @param s.reg the output of the R functiom 'lm' under the model of 'S' (the 2nd medistor)
#' @param q.red the output of the R function 'lm' under the model of 'Q' (the 1st medistor)

#' @export create_vcovmatrix

create_vcovmatrix=function(out.reg,s.reg,q.reg){
   n.beta<-length(out.reg$coefficients)
   n.alpha<-length(s.reg$coefficients)
   n.delta<-length(q.reg$coefficients)
   dimensions<-sum(n.beta+n.alpha+n.delta)
   vcovm<-matrix(0, nr=dimensions, nc=dimensions)
   vcovm[c(1:n.beta),c(1:n.beta)]<-vcov(out.reg)
   vcovm[c((1+n.beta):(n.beta+n.alpha)),c((1+n.beta):(n.beta+n.alpha))]<-vcov(s.reg)
   vcovm[c((1+n.beta+n.alpha):(n.beta+n.alpha+n.delta)),c((1+n.beta+n.alpha):(n.beta+n.alpha+n.delta))]<-vcov(q.reg)
   return(vcovm)
}
