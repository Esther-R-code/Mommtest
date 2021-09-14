#' the probability of the counterfactual outcome of the outcome
#'
#' This function allows you to calculate the probability of the counterfactual outcome when assume outcome=a and intervention=w.
#' Before using this function, please make sure the R packages "ordinal" and "parallel" have been install and libraries.
#' @title rho
#' @param TH theta; all of the parameters of the model #list(beta.hat, alpha.hat, delta.hat)
#' @param SH sigma: all of the estimates of the variance of the mediators
#' @param Name names of exposure and mediators, def:c("W","Q","S")
#' @param n_cate1: the number of categories of the outcome minus 1
#' @param a the value of outcome
#' @param w the values of intervention
#' @param confounders the values of confounders, def:NULL
#' @export rho
#' @keywords counterfactual
#' @return the probability of the counterfactual outcome
#' @examples
#' deltaq<-c(0.1,0.5) #Q c(d0,dw)
#' alphas<-c(0.2,0.5,0.5) #S c(al0,alw,alq)
#'
#' betay1<-c(-1.8, 1.0, 1.2, 1.3) #c(b01,bw1,bq1,bs1)
#' betay2<-c(0.8, 0.8, 1.1, 1.2) #c(b02,bw2,bq2,bs2)
#' betay<-rbind(betay1,betay2)
#' Tbeta<- NULL
#' for(beta.count in 1:4){Tbeta<-c(Tbeta,betay[,beta.count])}
#' Ttheta<- list(Tbeta, alphas, deltaq)
#' TH= Ttheta
#' SH= c(0.23,0.23)
#' n_cate1= 2
#' a= 2
#' w0000 <- c(0, 0, 0, 0)
#' confounders<- c()
#'
#' rho(TH, SH, Name=c("W","Q","S"), n_cate1, a, w0000, confounders)

rho<-function(TH,SH,Name=c("W","Q","S"),n_cate1,a,w,confounders=NULL){
  a0<-TH[[2]][1]
  aw<-TH[[2]][2]
  aq<-TH[[2]][3]
  d0<-TH[[3]][1]
  dw<-TH[[3]][2]
  prhob<-rep(0,length(TH[[1]]))
  n.M<-length(Name)
  b0a<-TH[[1]][paste0(a,"|",(a+1),".","(Intercept)")]; b0a_1<-TH[[1]][paste0((a-1),"|",a,".","(Intercept)")]
  for(m in 1:n.M){
     assign(paste0("b",Name[m],"a"),TH[[1]][paste0(a,"|",(a+1),".",Name[m])])
     assign(paste0("b",Name[m],"a_1"),TH[[1]][paste0((a-1),"|",a,".",Name[m])])
  }
  bCa<-NULL; bCa_1 <- NULL
  aC <- NULL; dC <- NULL
  n.C <- length(confounders)

  for(cc in 0:n.C){
     bCa <- c(bCa,TH[[1]][paste0(a,"|",(a+1),".",paste0("Con",cc))])
     bCa_1 <- c(bCa_1,TH[[1]][paste0((a-1),"|",a,".",paste0("Con",cc))])
     aC <- c(aC,TH[[2]][paste0("Con",cc)])
     dC <- c(dC,TH[[3]][paste0("Con",cc)])
  }
  bCa <- bCa[-1]; bCa_1 <- bCa_1[-1]
  aC <- aC[-1];   dC <- dC[-1]
  ss<-SH[1]; sq<-SH[2]

  #intervation 4
  if(length(w)==4){
  pmu.b.a <-b0a+bWa*w[1]+sum(bCa*confounders,na.rm=TRUE)   #b0a+bwa*wi+bcaT*C
  pmu.b.a_1 <-b0a_1+bWa_1*w[1]+sum(bCa_1*confounders,na.rm=TRUE)
  mu.q.l <- d0+dw*w[4]+sum(dC*confounders,na.rm=TRUE)  #d0+dw*wl+dCT*C
  mu.s <- a0+aw*w[2]+sum(aC*confounders,na.rm=TRUE)+aq*mu.q.l #a0+aw*wj+aq*(d0+dw*wl+dCT*C)+aCT*C
  mu.q.k <- d0+dw*w[3]+sum(dC*confounders,na.rm=TRUE)
  #b0a+bwa*wi+bcaT*C+bSa(a0+aw*wj+aq*(d0+dw*wl+dCT*C)+aCT*C)+bQa*(d0+dw*wl+dCT*C)=b0a+bwa*wi+bcaT*C+bSa*mu.s+bQa*pmu.s.k
  gtha <- pmu.b.a+bSa*mu.s+bQa*mu.q.k
  gtha_1 <- pmu.b.a_1+bSa_1*mu.s+bQa_1*mu.q.k
  ha <- sqrt(bSa^2*((aq^2)*(sq^2)+ss^2)+(bQa^2)*(sq^2)+1)
  ha_1 <- sqrt(bSa_1^2*(aq^2*sq^2+ss^2)+bQa_1^2*sq^2+1)
  Ga <- gtha/ha; Ga_1 <- gtha_1/ha_1
  ###############derivative########################
  pb0 <- (exp(-(Ga^2)/2)/sqrt(2*pi))/ha #prhoa/pb0a
  pb0_1 <- -(exp(-(Ga_1^2)/2)/sqrt(2*pi))/ha_1 #prhoa/pb0a_1
  pbw <- pb0*w[1]; pbw_1 <- pb0_1*w[1]
  pbq <- pb0*(mu.q.k-Ga*bQa*sq^2/ha)
  pbq_1 <- pb0_1*(mu.q.k-Ga_1*bQa_1*sq^2/ha_1)
  pbs <- pb0*(mu.s-Ga*bSa*(ss^2+aq^2*sq^2)/ha)
  pbs_1 <- pb0_1*(mu.s-Ga_1*bSa_1*(ss^2+aq^2*sq^2)/ha_1)

  pa0 <- sum(c(pb0,pb0_1)*c(bSa,bSa_1),na.rm=T)
  paw <- pa0*w[2]
  paq <- sum(c(pb0,pb0_1)*c(bSa*mu.q.l-(Ga*bSa^2*sq^2*aq/ha),bSa_1*mu.q.l-(Ga_1*bSa_1^2*sq^2*aq/ha_1)),na.rm=T)

  pd0 <- sum(c(pb0,pb0_1)*c((bSa*aq+bQa),(bSa_1*aq+bQa_1)),na.rm=T)
  pdw <- sum(c(pb0,pb0_1)*c((bSa*aq*w[4]+bQa*w[3]),(bSa_1*aq*w[4]+bQa_1*w[3])),na.rm=T)
  }
  prhobC <- rep(0,n.C*n_cate1)
  prhob0 <- rep(0,n_cate1); prhobW <- rep(0,n_cate1)
  prhobQ <- rep(0,n_cate1); prhobS <- rep(0,n_cate1)
  #prhob0: the derivative of rho with respect to intercept   #prhobW: the derivative of rho with respect to W
  #prhobQ: the derivative of rho with respect to Q           #prhobS: the derivative of rho with respect to S
  #prhobC: the derivative of rho with respect to covariates
  for(aa in 1:n_cate1){
    if(aa==(a-1)){
      for(k in 1:n.C){
        prhobC[(k-1)*n_cate1+(a-1)] <- pb0_1*confounders[k]
      }
      prhob0[(a-1)] <- pb0_1; prhobW[(a-1)] <- pbw_1
      prhobQ[(a-1)] <- pbq_1; prhobS[(a-1)] <- pbs_1
    }else if(aa==a){
      for(k in 1:n.C){
        prhobC[(k-1)*n_cate1+(a)] <- pb0*confounders[k]
      }
      prhob0[a] <- pb0; prhobW[a] <- pbw
      prhobQ[a] <- pbq; prhobS[a] <- pbs
    }else{}
  }
  prhob <- c(prhob0,prhobW,prhobQ,prhobS,prhobC)
  paC <- pa0*confounders
  pdC <- pd0*confounders

  derivatives=c(prhob, pa0, paw, paq, paC, pd0, pdw, pdC)
  if(is.na(Ga_1)){
  rhoa <- pnorm(Ga)
  }else if(is.na(Ga)){
  rhoa <- 1-pnorm(Ga_1)
  }else{rhoa <- pnorm(Ga)-pnorm(Ga_1)}

  return(c(rhoa,derivatives))
}
