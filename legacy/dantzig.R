library(flare)
dantz<-function (X, Y, zeta,lambda = NULL,Ulasso=FALSE, nlambda = NULL, nlabel= NULL,clam=1,bsdf=9,lambda.min.value = NULL, 
                 lambda.min.ratio = NULL, rho = 1, method = "dantzig", q = 2, res.sd = FALSE, 
                 prec = 1e-05, max.ite = 1e+05, verbose = F) 
{
  
  if (verbose)
    cat("Sparse Linear Regression with L1 Regularization.\\n")
  
  n = nrow(X)
  d = ncol(X)
  n1=length(Y)
  
  maxdf = max(n, d)
  xm = matrix(rep(colMeans(X), n), nrow = n, ncol = d, byrow = T)
  x1 = X - xm
  sdxinv = 1/sqrt(colSums(x1^2)/(n - 1))
  xx = x1 * matrix(rep(sdxinv, n), nrow = n, ncol = d, byrow = T)
  
  ym = mean(Y)
  y1 = Y - ym
  if (res.sd == TRUE) {
    sdy = sqrt(sum(y1^2)/(n1 - 1))
    yy = y1/sdy
  }
  else {
    sdy = 1
    yy = y1
  }
  
  
  intercept = FALSE
  if (intercept) {
    xx = cbind(rep(1, nrow(xx)), xx)
    X = cbind(rep(1, nrow(X)), X)
    d = d + 1
  }
  if (!is.null(lambda)) 
    nlambda = length(lambda)
  if (is.null(lambda)) {
    if (is.null(nlambda)) 
      nlambda = 8
    if (method == "dantzig") {
      if (intercept) 
        
        lambda.max = max(abs(crossprod(xx[, 2:d], yy/n)))
      else
      {
        # if (d>=n)
        #  lambda.max = max(abs(crossprod(xx[c(1:n1),], yy/n)))
        #else 
        lambda.max =clam* max(abs(zeta))
      }
    }
    
    if (method == "dantzig") {
      if (is.null(lambda.min.ratio)) {
        {if (d>=n) 
          lambda.min.ratio = 0.06
        else
          lambda.min.ratio = 0.04
        }
      }
      if (is.null(lambda.min.value)) {
        lambda.min.value = lambda.min.ratio * lambda.max
      }
    }
    else {
      if (is.null(lambda.min.value)) {
        lambda.min.value = sqrt(log(d)/n)
      }
      else {
        if (is.null(lambda.min.ratio)) {
          lambda.min.ratio = lambda.min.value/lambda.max
        }
      }
    }
    if (lambda.max < lambda.min.value) {
      lambda.max = 1
      lambda.min.value = 0.4
    }
    lambda = exp(seq(log(lambda.max), log(lambda.min.value), 
                     length = nlambda))
    rm(lambda.max, lambda.min.value, lambda.min.ratio)
    gc()
  }
  if (is.null(rho)) 
    rho = 1
  begt = Sys.time()
  if (method == "dantzig") {
    if (d >= n) 
      out =modified.dantzig.ladm.scr(zeta, xx, lambda, nlambda, 
                                     n, d , maxdf)
    else out = modified.dantzig.ladm.scr2(zeta, xx, lambda, nlambda, 
                                          n, d, maxdf)
    q = "infty"
  }
  
  runt = Sys.time() - begt
  df = rep(0, nlambda)
  if (intercept) {
    for (i in 1:nlambda) df[i] = sum(out$beta[[i]][2:d] != 0)
  }
  else {
    for (i in 1:nlambda) df[i] = sum(out$beta[[i]] != 0)
  }
  est = list()
  intcpt0 = matrix(0, nrow = 1, ncol = nlambda)
  intcpt = matrix(0, nrow = 1, ncol = nlambda)
  
  Xlabel=scale(X[c(1:n1),],scale=F)
  sdxinvlabel=1/sqrt(colSums(Xlabel^2)/(n1-1))
  
  if(Ulasso) { sdxinv=sdxinvlabel}
  
  if (intercept) {
    beta1 = matrix(0, nrow = d - 1, ncol = nlambda)
    for (k in 1:nlambda) {
      tmp.beta = out$beta[[k]][2:d]
      beta1[, k] = sdxinv * tmp.beta * sdy
      intcpt[k] = ym - as.numeric(xm[1, ] %*% beta1[, k]) + 
        out$beta[[k]][1] * sdy
      intcpt0[k] = intcpt[k]
    }
  }
  else {
    beta1 = matrix(0, nrow = d, ncol = nlambda)
    for (k in 1:nlambda) {
      tmp.beta = out$beta[[k]]
      intcpt0[k] = 0
      beta1[, k] = sdxinv * tmp.beta * sdy
      intcpt[k] = ym - as.numeric(xm[1, ] %*% beta1[, k])
    }
  }
  est$beta0 = out$beta
  est$beta = beta1
  est$intercept = intcpt
  est$intercept0 = intcpt0
  est$Y = Y
  est$X = X
  est$lambda = lambda
  est$nlambda = nlambda
  est$df = df
  est$method = method
  est$q = q
  est$ite = out$ite
  est$verbose = verbose
  est$runtime = runt
  class(est) = "slim"
  if (verbose) 
    print(est)
  return(est)
}

##X is the total data, n IS # OF total data 
##d is the number of covariates

modified.dantzig.ladm.scr<-function (zeta, X, lambda, nlambda, n, d, maxdf, rho=1,prec = 1e-05, max.ite = 1e+05,
                                     intercept=F, verbose=F) 
{
  if (verbose == TRUE) 
    cat("Dantzig selector with screening.\\n")
  XY = zeta
  XX = crossprod(X, X/n)
  beta = matrix(0, nrow = d, ncol = nlambda)
  ite.int = rep(0, nlambda)
  ite.int1 = rep(0, nlambda)
  ite.int2 = rep(0, nlambda)
  if (intercept) {
    intcep = 1
  }
  else {
    intcep = 0
  }
  if (n <= 3) {
    num.scr1 = n
    num.scr2 = n
  }
  else {
    num.scr1 = ceiling(n/log(n))
    num.scr2 = n - 1
  }
  order0 = order(abs(XY), decreasing = TRUE)
  idx.scr = order0
  num.scr = length(idx.scr)
  idx.scr1 = order0[1:num.scr1]
  idx.scr2 = order0[1:num.scr2]
  X1 = X[, idx.scr]
  XXX = crossprod(X1, crossprod(tcrossprod(X1, X1), X1))/(n^2)
  gamma = max(colSums(abs(XXX)))
  str = .C("slim_dantzig_ladm_scr", as.double(XY), as.double(XX), 
           as.double(XXX), as.double(beta), as.integer(n), as.integer(d), 
           as.double(rho), as.integer(ite.int), as.integer(ite.int1), 
           as.integer(ite.int2), as.integer(num.scr1), as.integer(num.scr2), 
           as.integer(idx.scr), as.integer(idx.scr1), as.integer(idx.scr2), 
           as.double(gamma), as.double(lambda), as.integer(nlambda), 
           as.integer(max.ite), as.double(prec), as.integer(intcep), 
           PACKAGE = "flare")
  beta.list = vector("list", nlambda)
  for (i in 1:nlambda) {
    beta.i = unlist(str[4])[((i - 1) * d + 1):(i * d)]
    beta.list[[i]] = beta.i
  }
  ite.int = unlist(str[8])
  ite.int1 = unlist(str[9])
  ite.int2 = unlist(str[10])
  ite = list()
  ite[[1]] = ite.int1
  ite[[2]] = ite.int2
  ite[[3]] = ite.int
  return(list(beta = beta.list, ite = ite))
}


modified.dantzig.ladm.scr2<-function (zeta, X, lambda, nlambda, n, d,  maxdf, rho=1,prec = 1e-05, max.ite = 1e+05,
                                      intercept=F, verbose=F) 
{
  if (verbose == TRUE) 
    cat("Dantzig selector with screening.\\n")
  XY = zeta
  XX = crossprod(X)/n
  beta = matrix(0, nrow = d, ncol = nlambda)
  ite.int = rep(0, nlambda)
  ite.int1 = rep(0, nlambda)
  if (intercept) {
    intcep = 1
  }
  else {
    intcep = 0
  }
  if (d <= 3) {
    num.scr1 = d
  }
  else {
    num.scr1 = ceiling(d/log(d))
  }
  order0 = order(abs(XY), decreasing = TRUE)
  idx.scr = order0
  num.scr = length(idx.scr)
  idx.scr1 = order0[1:num.scr1]
  XX1 = XX[idx.scr, idx.scr]
  XXX = crossprod(XX1, XX1)
  gamma = max(colSums(abs(X)))
  begt = Sys.time()
  str = .C("slim_dantzig_ladm_scr2", as.double(XY), as.double(XX), 
           as.double(XXX), as.double(beta), as.integer(n), as.integer(d), 
           as.double(rho), as.integer(ite.int), as.integer(ite.int1), 
           as.integer(num.scr1), as.integer(idx.scr), as.integer(idx.scr1), 
           as.double(gamma), as.double(lambda), as.integer(nlambda), 
           as.integer(max.ite), as.double(prec), as.integer(intcep), 
           PACKAGE = "flare")
  runt = Sys.time() - begt
  beta.list = vector("list", nlambda)
  for (i in 1:nlambda) {
    beta.i = unlist(str[4])[((i - 1) * d + 1):(i * d)]
    beta.list[[i]] = beta.i
  }
  ite.int = matrix(unlist(str[8]), byrow = FALSE, ncol = nlambda)
  ite.int1 = matrix(unlist(str[9]), byrow = FALSE, ncol = nlambda)
  ite = list()
  ite[[1]] = ite.int1
  ite[[2]] = ite.int
  return(list(beta = beta.list, ite = ite))
}
################

CreateF<-function (y, k = 10, list = TRUE, returnTrain = FALSE) 
{
  if (class(y)[1] == "Surv") 
    y <- y[, "time"]
  if (is.numeric(y)) {
    cuts <- floor(length(y)/k)
    if (cuts < 2) 
      cuts <- 2
    if (cuts > 5) 
      cuts <- 5
    breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
    y <- cut(y, breaks, include.lowest = TRUE)
  }
  if (k < length(y)) {
    y <- factor(as.character(y))
    numInClass <- table(y)
    foldVector <- vector(mode = "integer", length(y))
    for (i in 1:length(numInClass)) {
      min_reps <- numInClass[i]%/%k
      if (min_reps > 0) {
        spares <- numInClass[i]%%k
        seqVector <- rep(1:k, min_reps)
        if (spares > 0) 
          seqVector <- c(seqVector, sample(1:k, spares))
        foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
      }
      else {
        foldVector[which(y == names(numInClass)[i])] <- sample(1:k, 
                                                               size = numInClass[i])
      }
    }
  }
  else foldVector <- seq(along = y)
  if (list) {
    out <- split(seq(along = y), foldVector)
    names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))), 
                        sep = "")
    if (returnTrain) 
      out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
  }
  else out <- foldVector
  out
}

cv_dantzig<-function(nfolds=8,X,Y_all,y1,zetafull,nlambda=30,ratio=0.04,lambda.min.value=NULL)
{
  lambda.max =max(abs(zetafull))
  if(is.null(lambda.min.value)) lambda.min.value = ratio* lambda.max
  lambda = exp(seq(log(lambda.max), log(lambda.min.value), 
                   length = nlambda))
  loss_lambda<-matrix(0,nrow=nfolds,ncol=nlambda)
  Xcenter<-scale(X,center=T,scale=F)
  flds <- CreateF(y1, k = nfolds, list = TRUE, returnTrain = FALSE)
  
  for(i in 1:nfolds){
    
    test<-flds[[i]]
    n<-dim(y1)[1]
    subnumber<-n-length(test)
    Xtrain<-X[-test,]
    Ytrain<-Y_all[-test]
    
    X_new<-scale(Xtrain,center=T,scale=T)
    
    ytrain=scale(y1[-test,],center=T,scale=F)
    y2<-diag(as.numeric(ytrain))%*%X_new[c(1:subnumber),] 
    XX2<-diag(as.numeric(Ytrain))%*%X_new
    XX_new<-scale(XX2,center=T,scale=F)[c(1:subnumber),]
    
    zeta_n1<-colSums(y2-XX_new)/subnumber
    
    Xs<-scale(Xtrain,center=T,scale=F)
    
    dantzig_n<-dantz(Xs,ytrain,zeta_n1,lambda=lambda,verbose=F)
    beta_n<-dantzig_n$beta
    # intercept_n<-dantzig_n$intercept
    for(j in 1:nlambda)
    {
      y_predict=Xcenter[test,]%*%beta_n
      loss_lambda[i,j]<-sum((y1[test,]-y_predict[,j])^2)/length(test)
    }
  }  
  loss=colSums(loss_lambda)
  lambda.min=which.min(loss)
  return(list(lambda.min=lambda.min,lambda_val=lambda[lambda.min],loss=loss))
}


cv_dantzig2<-function(nfolds=8,X,y1,zetafull,nlambda=30,ratio=0.04)
{
  lambda.max =max(abs(zetafull))
  lambda.min.value = ratio* lambda.max
  lambda = exp(seq(log(lambda.max), log(lambda.min.value), 
                   length = nlambda))
  loss_lambda<-matrix(0,nrow=nfolds,ncol=nlambda)
  Xcenter<-scale(X,center=T,scale=F)
  flds <- CreateF(y1, k = nfolds, list = TRUE, returnTrain = FALSE)
  
  for(i in 1:nfolds){
    
    test<-flds[[i]]
    n<-dim(y1)[1]
    subnumber<-n-length(test)
    Xtrain<-X[-test,]
    X_new<-scale(Xtrain[c(1:subnumber),],center=T,scale=T)
    ytrain=scale(y1[-test,],center=T,scale=F)
    y2<-diag(as.numeric(ytrain))%*%X_new
    
    zeta_n1<-colSums(y2)/subnumber
    Xs<-scale(Xtrain,center=T,scale=F)
    
    dantzig_n<-dantz(Xs,ytrain,zeta_n1,lambda=lambda,verbose=F)
    beta_n<-dantzig_n$beta
    # intercept_n<-dantzig_n$intercept
    for(j in 1:nlambda)
    {
      y_predict=Xcenter[test,]%*%beta_n
      loss_lambda[i,j]<-sum((y1[test,]-y_predict[,j])^2)/n
    }
  }  
  loss=colSums(loss_lambda)
  lambda.min=which.min(loss)
  return(list(lambda.min=lambda.min,loss=loss))
}



cv_nonadd<-function(nfolds=8,XX2,Xs,y1,zetafull,BB,nlambda=50,ratio)
{
  lambda.max =max(abs(zetafull))
  lambda.min.value = ratio* lambda.max
  lambda = exp(seq(log(lambda.max), log(lambda.min.value), 
                   length = nlambda))
  loss_lambda<-matrix(0,nrow=nfolds,ncol=nlambda)
  
  flds <- CreateF(y1, k = nfolds, list = TRUE, returnTrain = FALSE)
  
  for(k in 1:nfolds){
    test<-flds[[k]]
    n<-dim(y1)[1]
    subnumber<-n-length(test)
    Xtrain<-XX2[-test,]
    Ytrain<-y1[-test]
    Xstrain<-Xs[-test,]
    XX_new<-scale(Xtrain,scale=F)[c(1:subnumber),]
    zeta_new<-colMeans(diag(as.numeric(Ytrain))%*%Xstrain)
    
    for(i in 1:p)
    {
      zeta_new[i]<-zeta_new[i]-(colMeans(XX_new)%*% BB[,i])
    }
    sdxinv = 1/sqrt(colSums(Xstrain^2)/(subnumber - 1))
    zeta_n<-zeta_new*sdxinv
    dantzig_n<-dantz(Xstrain,Ytrain,zeta_n,lambda=lambda,clam=1,verbose=F)
    beta_n<-dantzig_n$beta
    
    for(j in 1:nlambda)
    {
      y_predict=Xs[test,]%*%beta_n
      loss_lambda[k,j]<-sum((y1[test,]-y_predict[,j])^2)/n
    }
  }
  loss=colSums(loss_lambda)
  lambda.min=which.min(loss)
  return(list(lambda.min=lambda.min,loss=loss))
}

cv_lasso<-function(nfolds=5,X,Y_all,y1,lambdaminratio=0.05)
{ n=length(y1)
N=length(Y_all)-n
ylabel=(n+N)/n*(y1-Y_all[1:n])+Y_all[1:n]
imputey=c(ylabel,Y_all[(n+1):(n+N)])
Xs2<-scale(X,scale=F)
object=glmnet(Xs2,imputey,family = "gaussian",standardize = F,intercept = F,lambda.min.ratio=lambdaminratio)
lambda =object$lambda
nlambda=length(lambda)
loss_lambda<-matrix(0,nrow=nfolds,ncol=nlambda)
Xcenter<-scale(X,center=T,scale=F)
flds <- CreateF(y1, k = nfolds, list = TRUE, returnTrain = FALSE)

for(i in 1:nfolds){
  
  test<-flds[[i]]
  n<-dim(y1)[1]
  subnumber<-n-length(test)
  Xtrain<-X[-test,]
  Ytrain<-Y_all[-test]
  ylabel=(subnumber+N)/subnumber*(y1[-test]-Ytrain[1:subnumber])+Ytrain[1: subnumber]
  imputey=c(ylabel,Ytrain[(subnumber+1):(subnumber+N)])
  Xs_train<-scale(Xtrain,scale=F)
  object_train=glmnet(Xs_train,imputey,family = "gaussian",standardize = F,intercept = F,lambda=lambda)
  beta_n=object_train$beta
  # intercept_n<-dantzig_n$intercept
  for(j in 1:nlambda)
  {
    y_predict=Xcenter[test,]%*%beta_n
    loss_lambda[i,j]<-sum((y1[test,]-y_predict[,j])^2)/length(test)
  }
}  
loss=colSums(loss_lambda)
lambda.min=which.min(loss)
return(list(lambda.min=lambda.min,lambda_val=lambda[lambda.min],loss=loss))
}


