

simrun<-function(X,y1,df=5,ratio,alpha=0.05,nlam=10,beta_ll)
  
{ n=length(y1)
  Xs<-scale(X[c(1:n),],scale=F)
  Xs2<-scale(X,scale=F)
  N=nrow(X)-n
  p=ncol(X)
  Y_all<-c()
  bsMat<-c()
  for(i in 1:p){
    bsmat<-bs(X[,i],df=df,degree=2)
    bsMat<-cbind(bsMat,bsmat)
  }
  bsMat1<-scale(bsMat,scale=F)
  index<-c()
  for (i in 1:p)  index<-c(index,rep(i,df))
  XX<-bsMat1[c(1:n),]
  lambda.grp<- lambdamax(XX, y = y1, index = index, penscale = sqrt,center=F,
                         model = LinReg()) * 0.8^(0:nlam)
  fit1 <- grplasso(XX, y = y1, index = index, lambda = lambda.grp, model = LinReg(),
                   penscale = sqrt,center=F,
                   control = grpl.control(update.hess = "lambda", trace = 0))
  
  clambda<-gBIC(fit1,y1)
  beta_cc<-fit1$coefficients[,clambda]
  nonempty<-ceiling(which(fit1$coefficients[,clambda]!=0)/df)
  nonempty<-unique(nonempty)
  Y_all<-bsMat1%*%beta_cc
  
  sum((y1-Y_all[c(1:n)])^2)/(n)
  sum((y1-Xs[c(1:n),]%*%beta_ll)^2)/(n)
  #  sum((Ytrue-Y_all)^2)/(n+N)
  noise<-sum((Y_all[c(1:n)]-y1)^2)/(n-df*length(nonempty))
  
  newy<-as.numeric(y1-Xs[c(1:n),]%*%beta_ll)
  sdxinv = 1/sqrt(colSums(Xs^2)/(n - 1))
  newY<-diag(newy)%*%(Xs[c(1:n),]%*%diag(sdxinv))
  # 
  XX2<-diag(as.numeric(Y_all))%*%(Xs2%*%diag(sdxinv))
  XX_new<-scale(XX2,scale=F)[c(1:n),]
  XXs<-scale(XX2[c(1:n),],scale=F)
  zeta_new<-colMeans(diag(as.numeric(y1))%*%(Xs%*%diag(sdxinv)))
  zeta_original<-zeta_new

  BB<-matrix(0,nrow=p,ncol=p)
  for(i in 1:p)
  {
    y_new=newY[,i]
    fit_new= glmnet(XXs,y_new,lambda=
                      cv.glmnet(XXs,y_new,family="gaussian",nfolds=5)$lambda.1se)
    BB[,i]<-coef(fit_new)[-1]
    zeta_new[i]<-zeta_new[i]-(colMeans(XX_new)%*%coef(fit_new)[-1])
  }
  
  X_n1<-scale(X[c(1:n),],center=T,scale=T)
  
  X_new<-scale(X,center=T,scale=T)
  sigma_n<-(1/(n+N))*(t(X_new)%*%X_new) 
  M <- InverseLinfty(sigma_n, n+N, resol=1.2, maxiter=50, threshold=1e-2, verbose=T)

  y2<-diag(as.numeric(y1))%*%X_new[c(1:n),] 
  XX_2<-diag(as.numeric(Y_all))%*%X_new
  XXnew<-scale(XX_2,center=T,scale=F)[c(1:n),]
  zeta_n<-colSums(y2-XXnew)/n
  
  dantzig_n<-dantz(X,y1,zeta_n,lambda.min.ratio=ratio,nlambda=50,clam=1,verbose=T)
  beta_n<-dantzig_n$beta
  lambdamin<-cv_dantzig(nfolds=5,X=X,Y_all=Y_all,y1=y1,zetafull=zeta_n,nlambda=50,ratio=ratio)
  clam_n=lambdamin$lambda.min
  e1<-beta_n[,clam_n]-beta_t
  error1<-e1^2
  error2<-abs(e1)
  y3<-diag(as.numeric(y1))%*%X_n1
  zeta_n1<-colMeans(y3)
  dantzig_n1<-dantz(X,y1,zeta_n1,Ulasso=T,lambda.min.ratio=0.2,nlambda=50,clam=1,verbose=T)
  beta_n1<-dantzig_n1$beta
  lambdamin1<-cv_dantzig2(nfolds=5,X=X,y1=y1,zetafull=zeta_n1,nlambda=50,ratio=0.2)
  clam_n1=lambdamin1$lambda.min
  e2<-beta_n1[,clam_n1]-beta_t
  error3<-e2^2
  error4<-abs(e2)
  
  
  dantzig_n2<-dantz(X[c(1:n),],y1,zeta_new,lambda.min.ratio=ratio,nlambda=50,clam=1,verbose=T)
  beta_n2<-dantzig_n2$beta
  chooselam<-cv_nonadd(nfolds=5,XX2=XX2,Xs=Xs,y1=y1,zetafull=zeta_new,BB=BB,nlambda=50,ratio=ratio)
  clam<-chooselam$lambda.min
  e4<-beta_n2[,clam]-beta_t
  error5<-e4^2
  error6<-abs(e4)
  
  #   # 
  Gamma<-cov(X_new*as.numeric(Y_all-Xs2%*%beta_ll))
  htheta<-beta_ll/(sdxinv)
  sunbiased.Lasso <- as.numeric((htheta) + (M%*%(zeta_n-sigma_n %*% (htheta))))
  B<-M %*% Gamma %*% t(M)
  A <- M %*% sigma_n %*% t(M)
  interval.size <- qnorm(1-(alpha/2))*sqrt(noise*diag(A)+(n/(n+N))*diag(B))/sqrt(n);
  interval.size<-interval.size * sdxinv
  sunbiased.Lasso<-sunbiased.Lasso *sdxinv
  lowerbd<-sunbiased.Lasso-interval.size
  upperbd<-sunbiased.Lasso+interval.size
  dif1<-interval.size
  coveragel1<- as.numeric(lowerbd<=beta_t+0.05 & upperbd+0.05>=beta_t)
  
  
  htheta2<-beta_ll
  sunbiased.Lasso1 <- as.numeric((htheta2) + (M%*%(zeta_new-(1/n)*(t(X_n1)%*%Xs)%*%htheta2))*sdxinv)
  AA<-diag(as.numeric(y1-Xs%*%beta_ll))%*%X_n1
  
  var2<-M%*%cov(AA)%*%t(M)
  sunbiased.Lasso2 <- as.numeric((beta_ll) + (M%*%(colMeans(AA)))*sdxinv)
  Ar<-(AA-XXs%*%BB)
  var1<-M%*%(N/(n+N)*cov(Ar))%*%t(M)+n/(n+N)*M%*%cov(AA)%*%t(M)
  #
  se2<-sqrt(diag(var2))/sqrt(n)*sdxinv
  se1<-sqrt(diag(var1))/sqrt(n)*sdxinv
  diff1<-qnorm(1-(alpha/2))*se1
  diff2<-qnorm(1-(alpha/2))*se2
  lowerbd1<- sunbiased.Lasso1-qnorm(1-(alpha/2))*se1
  lowerbd2<- sunbiased.Lasso2-qnorm(1-(alpha/2))*se2
  upperbd2<- sunbiased.Lasso2+qnorm(1-(alpha/2))*se2
  upperbd1<- sunbiased.Lasso1+qnorm(1-(alpha/2))*se1
  coverage1<- as.numeric( beta_t-lowerbd1>-0.05 & upperbd1-beta_t>-0.05)
  coverage2<-as.numeric( beta_t-lowerbd2>-0.05 & upperbd2-beta_t>-0.05)
  
  Zfit<-calculate.Z(X_new,verbose=T,Z=NULL,debug.verbose = F)
  Z<-scale(Zfit$Z, center = FALSE, scale = 1/Zfit$scaleZ)
  ZZ<-Zfit$Z
  object2<-lasso.proj(X[c(1:n),],y1,Z=Z[c(1:n),],robust=T)
  y2<-diag(as.numeric(y1))%*%ZZ[c(1:n),]
  XX2<-diag(as.numeric(Y_all))%*%ZZ
  XX_new<-scale(XX2,center=T,scale=F)[c(1:n),]
  zeta_n1<-colSums(y2-XX_new)/n
  htheta<-beta_n[,clam_n]/sdxinv
  sunbiasedL<-as.numeric((htheta) + (zeta_n1-(t(ZZ)%*%X_new%*% (htheta))/(n+N)))
  #
  A1<- t(ZZ[c(1:n),])%*%ZZ[c(1:n),]/n
  B1<- t(ZZ)%*%diag(as.numeric(Y_all-Xs2%*%beta_ll)^2)%*%ZZ/(n+N)
  intervals <- qnorm(1-(alpha/2))*sqrt(noise*diag(A1)+(n/(n+N))*diag(B1))/sqrt(n);
  intervals<-intervals *sdxinv
  sunbiasedL<-sunbiasedL *sdxinv
  lowerbdl1<-sunbiasedL-intervals
  upperbdl1<-sunbiasedL+intervals
  coveragel2<- as.numeric( beta_t-lowerbdl1>-0.05 & upperbdl1-beta_t>-0.05)
  dif2<-intervals
  #
  # dif3<-qnorm(1-(alpha/2))*object2$se
  # objectll2<-object2$bhat-qnorm(1-(alpha/2))*object2$se
  # objectul2<-object2$bhat+object2$se*qnorm(1-(alpha/2))
  # coveragel3<-as.numeric(objectll2< beta_t+0.05 & objectul2+0.05 >=beta_t)
  
  re1<-rbind(error1,error2,error3,error4,error5,error6,coverage1,coverage2,coveragel1,coveragel2,diff1,diff2,dif1,dif2)
  return(re1)
}


##############cross fitting

simrun2<-function(X,y1,df=5,nfolds=4,alpha=0.05,nlam=10,beta_ll)
{
  Xs<-scale(X[c(1:n),],scale=F)
  Xs2<-scale(X,scale=F)
  n=length(y1)
  N=dim(X)[1]-n
  flds <- split(sample(n,n,replace=F),1:nfolds)
  flds2 <- split(sample(N,N,replace=F),1:nfolds)
  
  k=1
  Y_all=rep(0,n+N)
  test<-c(1:n)[-flds[[k]]]
  Ntest<-c(1:N)[-flds2[[k]]]
  Ntest=Ntest+n
  
  cross<-getzeta(X,y1,test,Ntest,beta_ll=beta_ll)
  Y_all[-c(test,Ntest)]<-cross$Y_all
  
  
  for(k in 2:nfolds){
    test<-c(1:n)[-flds[[k]]]
    Ntest<-c(1:N)[-flds2[[k]]]
    Ntest=Ntest+n
    val=paste('folds',k,sep = "")
    crossfit<-getzeta(X,y1,test,Ntest,beta_ll=beta_ll)
    cross<- lapply(seq_along(crossfit),function(i)
      (crossfit[[i]]+cross[[i]]) )
    Y_all[-c(test,Ntest)]<-crossfit$Y_all
  }
  
  
  zeta_n=cross[[2]]/nfolds
  zeta_new=cross[[4]]/nfolds
  noise=cross[[3]]/nfolds
  BB= cross[[7]]/nfolds
  Gamma= cross[[5]]/nfolds
  
  sdxinv = 1/sqrt(colSums(Xs^2)/(n - 1))
  
  X_n1<-scale(X[c(1:n),],center=T,scale=T)
  X_new<-scale(X,center=T,scale=T)
  sigma_n<-(1/(n+N))*(t(X_new)%*%X_new) 
  M <- InverseLinfty(sigma_n, n+N, resol=1.2, maxiter=50, threshold=1e-2, verbose=T)
  #   M<-solve(sigma_n)
  #   M[abs(M)<0.05]<-0
  ###### if wantTheta is True, the method performed is Oldschool, which invovles the inaccurate prediction.
  
  #   
  #   dantzig_n<-dantz(X,y1,zeta_n,lambda.min.ratio=ratio,nlambda=50,clam=1,verbose=T)
  #   beta_n<-dantzig_n$beta
  #   lambdamin<-cv_dantzig(nfolds=5,X=X,Y_all=Y_all,y1=y1,zetafull=zeta_n,nlambda=50,ratio=ratio)
  #   clam_n=lambdamin$lambda.min
  #   e1<-beta_n[,clam_n]-beta_t
  #   error1<-e1^2
  #   error2<-abs(e1)
  #   y3<-diag(as.numeric(y1))%*%X_n1
  #   zeta_n1<-colMeans(y3)
  #   #########replicate zeta_n1
  #   dantzig_n1<-dantz(X,y1,zeta_n1,Ulasso=T,lambda.min.ratio=0.15,nlambda=50,clam=1,verbose=T)
  #   beta_n1<-dantzig_n1$beta
  #   lambdamin1<-cv_dantzig2(nfolds=5,X=X,y1=y1,zetafull=zeta_n1,nlambda=50,ratio=0.15)
  #   clam_n1=lambdamin1$lambda.min
  #   e2<-beta_n1[,clam_n1]-beta_t
  #   error3<-e2^2
  #   error4<-abs(e2)
  #   fitl2= glmnet(Xs,y1,family = "gaussian",lambda=cv.glmnet(Xs,y1,family="gaussian")$lambda.min)
  #   beta1<-coef(fitl2)[-1]
  #   e3<-beta1-beta_t
  #   error5<-e3^2
  #   error6<-abs(e3)
  #   
  # #  zeta_new=(cross1$zeta_new+cross2$zeta_new)/2
  #   dantzig_n2<-dantz(X[c(1:n),],y1,zeta_new,lambda.min.ratio=ratio,nlambda=50,clam=1,verbose=T)
  #   beta_n2<-dantzig_n2$beta
  # #  BB=(cross1$BB+cross2$BB)/2
  #   chooselam<-cv_nonadd(nfolds=5,XX2=XX2,Xs=Xs,y1=y1,zetafull=zeta_new,BB=BB,nlambda=50,ratio=ratio)
  #   clam<-chooselam$lambda.min
  #   e4<-beta_n2[,clam]-beta_t
  #   error7<-e4^2
  #   error8<-abs(e4)


  htheta<-beta_ll/(sdxinv)
  sunbiased.Lasso <- as.numeric((htheta) + (M%*%(zeta_n-sigma_n %*% (htheta))))
  B<-M %*% Gamma %*% t(M)
  A <- M %*% sigma_n %*% t(M)
  interval.size <- qnorm(1-(alpha/2))*sqrt(noise*diag(A)+(n/(n+N))*diag(B))/sqrt(n);
  interval.size<-interval.size * sdxinv
  sunbiased.Lasso<-sunbiased.Lasso *sdxinv
  lowerbd<-sunbiased.Lasso-interval.size
  upperbd<-sunbiased.Lasso+interval.size
  dif1<-interval.size
  coveragel1<- as.numeric(lowerbd<=beta_t+0.05 & upperbd+0.05>=beta_t)
  
  
  htheta2<-beta_ll
  sunbiased.Lasso1 <- as.numeric((htheta2) + (M%*%(zeta_new-(1/n)*(t(X_n1)%*%Xs)%*%htheta2))*sdxinv)
  
  XX2<-diag(as.numeric(Y_all))%*%(Xs2%*%diag(sdxinv))
  XX_new<-scale(XX2,scale=F)[c(1:n),]
  XXs<-scale(XX2[c(1:n),],scale=F)
  
  AA<-diag(as.numeric(y1-Xs%*%beta_ll))%*%X_n1
  Ar=AA-XXs%*%BB  
  var2<-M%*%cov(AA)%*%t(M)
  sunbiased.Lasso2 <- as.numeric((beta_ll) + (M%*%(colMeans(AA)))*sdxinv)
  
  #    cov_Ar<-(cov(cross1$Ar)+cov(cross2$Ar))/2
  var1<-M%*%(N/(n+N)*cov(Ar))%*%t(M)+n/(n+N)*M%*%cov(AA)%*%t(M)
  
  
  se2<-sqrt(diag(var2))/sqrt(n)*sdxinv
  se1<-sqrt(diag(var1))/sqrt(n)*sdxinv
  diff1<-qnorm(1-(alpha/2))*se1
  diff2<-qnorm(1-(alpha/2))*se2
  lowerbd1<- sunbiased.Lasso1-qnorm(1-(alpha/2))*se1
  lowerbd2<- sunbiased.Lasso2-qnorm(1-(alpha/2))*se2
  upperbd2<- sunbiased.Lasso2+qnorm(1-(alpha/2))*se2
  upperbd1<- sunbiased.Lasso1+qnorm(1-(alpha/2))*se1
  coverage1<- as.numeric( beta_t-lowerbd1>-0.05 & upperbd1-beta_t>-0.05)
  coverage2<-as.numeric( beta_t-lowerbd2>-0.05 & upperbd2-beta_t>-0.05)
  
  
  Zfit<-calculate.Z(X_new,verbose=T,Z=NULL,debug.verbose = F)
  Z<-scale(Zfit$Z, center = FALSE, scale = 1/Zfit$scaleZ)
  ZZ<-Zfit$Z
  object2<-lasso.proj(X[c(1:n),],y1,Z=Z[c(1:n),],robust=T)
  y2<-diag(as.numeric(y1))%*%ZZ[c(1:n),]
  XX2<-diag(as.numeric(Y_all))%*%ZZ
  XX_new<-scale(XX2,center=T,scale=F)[c(1:n),]
  zeta_n1<-colSums(y2-XX_new)/n
  htheta<-beta_n[,clam_n]/sdxinv
  sunbiasedL<-as.numeric((htheta) + (zeta_n1-(t(ZZ)%*%X_new%*% (htheta))/(n+N)))
  #
  A1<- t(ZZ[c(1:n),])%*%ZZ[c(1:n),]/n
  B1<- t(ZZ)%*%diag(as.numeric(Y_all-Xs2%*%beta_ll)^2)%*%ZZ/(n+N)
  intervals <- qnorm(1-(alpha/2))*sqrt(noise*diag(A1)+(n/(n+N))*diag(B1))/sqrt(n);
  intervals<-intervals *sdxinv
  sunbiasedL<-sunbiasedL *sdxinv
  lowerbdl1<-sunbiasedL-intervals
  upperbdl1<-sunbiasedL+intervals
  coveragel2<- as.numeric( beta_t-lowerbdl1>-0.05 & upperbdl1-beta_t>-0.05)
  dif2<-intervals
  #
  dif3<-qnorm(1-(alpha/2))*object2$se
  objectll2<-object2$bhat-qnorm(1-(alpha/2))*object2$se
  objectul2<-object2$bhat+object2$se*qnorm(1-(alpha/2))
  coveragel3<-as.numeric(objectll2<beta_t+0.05 & objectul2+0.05 >=beta_t)
  re1<-rbind(coverage1,coverage2,coveragel1,coveragel2,diff1,diff2,dif1,dif2)
  return(re1)
}
