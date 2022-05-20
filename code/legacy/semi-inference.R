######
# Simulation for optimal and safe semi-supervised inference
#
# 1.Two ways of deriving \hat Omega: the inverse sample covariance matrix are considered:
#   a. from debiased lasso by Andrea Montanari
#   b. with nodewise lasso in the r package 'hdi'
# 2.Additive model is performed to impute the unlabeled data.
# 3.The main function for deriving optimal and safe estimators is 'dantz' contained in 'dantzig.R', which is based on 
#   the 'slim' function in 'flare' and takes zeta as an input to adapt to our estimation problem.
#   The input ratio in the function is the ratio of minimal lambda over maximal lambda (tuning parameter). Don't set 
#   it too small, since when it is too small, 'flare' takes too long time to run.
# 4.Change the generative model in data_generate4 contained in nec.R.
# 5.output 
#    1) error1,3,5,7 are l2 errors for Optimal, U-dantzig, Lasso, and safe estimators.
#    2) error2,4,6,8 are l1 errors for Optimal, U-dantzig, Lasso, and safe estimators.
#    3) dif1, dif2 are lengths of intervials for optimal estimator with method a and b.
#    4) diff1 is the length of intervial for safe estimator with method a.
#    5) diff2, diff3 and dif3 are length of intervial for lasso using different methods.

source("nec.R")

n=200
p=500
alpha=0.05
N=200
ratio=0.08
nlam=10
df=5
an<-30000
set.seed(7)
data_t<-data_generate4(r=0.3,n,N=an,p)
ay<-data_t$y
ay<-scale(ay,center=T,scale=F)
ax<-data_t$X
ax<-scale(ax,center=T,scale=F)
beta_t<-solve(t(ax)%*%ax)%*%t(ax)%*%ay
beta_t[which(abs(beta_t)<0.1)]=0
beta_t<-rep(0,p)
beta_t[1:6]<-c(1.1,0,2.4,4,4,2)

registerDoParallel(cores=4)

result1<-foreach (repe=1:20,.combine = "rbind") %dopar%
  { 
    data1<-data_generate4(r=0.3,n,N,p)
    X<-data1$X
    Ytrue=scale(data1$y,scale=F)
    y1<-scale(data1$y[c(1:n)],center=T,scale=F)
    Xs<-scale(X[c(1:n),],scale=F)
    Xs2<-scale(X,scale=F)
    fit_b= glmnet(Xs,y1,family = "gaussian",lambda=
                    cv.glmnet(Xs,y1,family="gaussian")$lambda.1se)
    beta_ll<-coef(fit_b)[-1]
    ################################
    #  Different methods of imputing unlabled data
    #
    #  Xwork<-cbind(X,X[,1]*X[,2],X[,3]*X[,4],X[,3]^2,X[,4]^2,X[,5]^2,X[,6]^2)
    #  Xws<-scale(Xwork[c(1:n),],scale=F)
    #  fit_w= glmnet(Xws,y1,family = "gaussian",lambda=cv.glmnet(Xws,y1,family="gaussian")$lambda.1se)
    #  beta_w<-coef(fit_w)[-1]
    #  Y_all<-Xws%*%beta_w
    #  Y_all<-scale(Xwork,scale=F)%*%beta_w
    
    ###  Impute with additive model
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
    
    #  sum((y1-Y_all[c(1:n)])^2)/(n)
    #  sum((y1-Xs[c(1:n),]%*%beta_ll)^2)/(n)
    #  sum((Ytrue-Y_all)^2)/(n+N)
    
    noise<-sum((Y_all[c(1:n)]-y1)^2)/(n-df*length(nonempty))
    
    ####zeta_new associates safe debiased estimator    
    newy<-as.numeric(y1-Xs[c(1:n),]%*%beta_ll)
    sdxinv = 1/sqrt(colSums(Xs^2)/(n - 1))
    newY<-diag(newy)%*%(Xs[c(1:n),]%*%diag(sdxinv))
    #
    XX2<-diag(as.numeric(Y_all))%*%(Xs2%*%diag(sdxinv))
    XX_new<-scale(XX2,scale=F)[c(1:n),]
    XXs<-scale(XX2[c(1:n),],scale=F)
    
    zeta_new<-colMeans(diag(as.numeric(y1))%*%(Xs%*%diag(sdxinv)))
    zeta_original<-zeta_new
    p=ncol(X)
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
    #   M<-solve(sigma_n)
    #   M[abs(M)<0.05]<-0
    ###### zeta_n associates with optimal estimator
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
    
    y3<-diag(as.numeric(y1))%*%X_new[c(1:n),] 
    zeta_n1<-colMeans(y3)
    dantzig_n1<-dantz(X,y1,zeta_n1,Ulasso=T,lambda.min.ratio=0.15,nlambda=50,clam=1,verbose=T)
    beta_n1<-dantzig_n1$beta
    lambdamin1<-cv_dantzig(nfolds=5,X=X,Y_all=rep(0,n+N),y1=y1,zetafull=zeta_n1,nlambda=50,ratio=0.15)
    clam_n1=lambdamin1$lambda.min
    e2<-beta_n1[,clam_n1]-beta_t
    error3<-e2^2
    error4<-abs(e2)
    
    fitl2= glmnet(Xs,y1,family = "gaussian",lambda=cv.glmnet(Xs,y1,family="gaussian")$lambda.min)
    beta1<-coef(fitl2)[-1]
    e3<-beta1-beta_t
    error5<-e3^2
    sum(error5)
    error6<-abs(e3)
    
    dantzig_n2<-dantz(X[c(1:n),],y1,zeta_new,lambda.min.ratio=ratio,nlambda=50,clam=1,verbose=T)
    beta_n2<-dantzig_n2$beta
    chooselam<-cv_nonadd(nfolds=5,XX2=XX2,Xs=Xs,y1=y1,zetafull=zeta_n,BB=BB,nlambda=50,ratio=ratio)
    clam<-chooselam$lambda.min
    e4<-beta_n2[,clam]-beta_t
    error7<-e4^2
    error8<-abs(e4)
    
    ####
    #### CI for optimal     
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
    
    ## CI for safe estimator and original lasso estimator
    ## diff1 is for safe and diff2 is for lasso with hat Omega estimated with all data.
    htheta2<-beta_n2[,clam]
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
    
    ## diff3 is the result given by hdi using only labeled data with robust option.
    object<-lasso.proj(X[c(1:n),],y1,robust=T)
    diff3<-qnorm(1-(alpha/2))*object$se
    objectll<-object$bhat-qnorm(1-(alpha/2))*object$se
    objectul<-object$bhat+object$se*qnorm(1-(alpha/2))
    coverage3<-as.numeric(objectll<beta_t+0.05 & objectul+0.05 >=beta_t)
    
    ## The other way to get hat Omega with function in hdi. The output Z is hat Omega%*%X.
    ## dif2 is for the optimal estimator with Z and dif3 is for the original lasso.
    ## 
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
    
    re1<-rbind(error1,error2,error3,error4,error5,error6,error7,error8,coverage1,coverage2,coverage3,coveragel1,coveragel2,coveragel3,diff1,diff2,diff3,dif1,dif2,dif3)
    
    datt<-data_generate4(r=0.3,200,400,p)
    X<-rbind(X,datt$X)
    re2<-simrun(X=X,y1=y1,ratio=0.1,beta_ll=beta_ll)
    datt1<-data_generate4(r=0.3,400,400,p)
    X<-rbind(X,datt1$X)
    re3<-simrun(X=X,y1=y1,ratio=0.1,beta_ll=beta_ll)
    rbind(re1,re2,re3)
  }  
################################
################################


  



