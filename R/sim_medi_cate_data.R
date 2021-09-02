#'this function to generate simulation data
#'
#'how the mediators and exposures are generated
#'W=exposure, Q=mediator one right after W, S=second mediator closest to outcome, Y=outcome
#'@title simulate a data
#'@param sample_size the number of case wanted to be create
#'@param exp_cont the mean and standard deviations of exposure; only consider continuous exposure
#'@param med_cont vector of standard deviations of the mediators; only consider continous mediators
#'@param n_cate the number of categories of the outcome
#'@param deltaq the true value of the parameter of the model of Q
#'@param alphas the true value of the parameter of the model of S
#'@param betay the true value of the parameter of the model of Y
#'@param coa the true value of mean and variance of the covariates seeting as c(mean(coa1),sd(coa1),mean(coa2),sd(coa2),...); default:coa=C()
#'@param j is the number to set seed; set.seed(j)
#'@export sim_medi_cate_data
#'@examples
#'sample_size<- 500
#'exp_cont<- c(0,1)
#'med_cont<- c(0.23,0.23)
#'n_cate<- 3
#'deltaq<- c(0.1,0.5) #Q c(d0,dw)
#'alphas<- c(0.2,0.5,0.5) #S c(al0,alw,alq)
#'betay1<- c(-1.8, 1.0, 1.2, 1.3, 0.05) #c(b01,bw1,bq1,bs1,bc1)
#'betay2<- c(0.8, 0.8, 1.1, 1.2, 0.03) #c(b02,bw2,bq2,bs2,bc2)
#'betay<- rbind(betay1,betay2)
#'coa<- c(0,0.4)
#'j<- 56 #sample(1,1:1e+09)
#'sim_cover<- sim_medi_cate_data(sample_size=sample_size,exp_cont=exp_cont,med_cont=med_cont,n_cate=n_cate,deltaq=deltaq,alphas=alphas,betay=betay,coa=coa,j)
#'@return Y the outcome
#'@return W the exposure
#'@return Q the 1st mediator
#'@return S the 2nd mediator
#'@return Coa1 the 1st covariate

sim_medi_cate_data=function(sample_size, exp_cont, med_cont, n_cate, deltaq, alphas, betay, coa=c(), j){
  set.seed(j)
  if(length(coa)==0){
      n_coa<- 0
      W<- rnorm(sample_size,mean=exp_cont[1],sd=exp_cont[2])
      Q<- deltaq[1]+deltaq[2]*W+rnorm(sample_size,mean=0,med_cont[1])
      S<- alphas[1]+alphas[2]*W+alphas[3]*Q+rnorm(sample_size,mean=0,med_cont[2])
      X.matrix<- cbind(rep(1,sample_size), W, Q, S)
    }else{
      n_coa<- length(coa)/2
      Coa<- NULL
      Coa_name<- NULL
      for(count_coa in 1:n_coa){
        Coa<-cbind(Coa,rnorm(sample_size,mean=coa[(1+(n_coa-1)*2)],sd=coa[(2+(n_coa-1)*2)]))
        Coa_name<- c(Coa_name,paste0("Coa", count_coa))
        }
      colnames(Coa)<-Coa_name
      W<- rnorm(sample_size, mean=exp_cont[1], sd=exp_cont[2])
      Q<- deltaq[1]+deltaq[2]*W+rowSums(t(t(Coa)*deltaq[-(1:2)]))+rnorm(sample_size, mean=0, med_cont[1])
      S<- alphas[1]+alphas[2]*W+alphas[3]*Q+rowSums(t(t(Coa)*alphas[-(1:3)]))+rnorm(sample_size, mean=0, med_cont[2])
      X.matrix<- cbind(rep(1,sample_size), W, Q, S, Coa)
    }

    eta<- X.matrix%*%t(betay)
    CP<- pnorm(eta)#link function pnorm

    if(n_cate>2){
    CP2<- cbind(CP[,2:(n_cate-1)],rep(1,sample_size))
    P<- cbind(CP[,1],CP2-CP)
    }else{
    P<- cbind(CP,rep(1,sample_size)-CP)
    }

    Y<- as.matrix(rep(0,sample_size))
    for(i in 1:sample_size){Y[i,]<- (which(rmultinom(1, size = 1, prob = P[i,])==1))}
    Y<- factor(Y)

    if(sum(which(Y==0))!=0){
        print("Error! Number of Y not equal to the sample size.")
    }else{
        df<- data.frame(Y,X.matrix[,-1])
        return(df)
    }
}

