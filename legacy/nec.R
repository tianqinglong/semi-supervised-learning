###
###
###
library(grplasso)
library(flare)
library(spam)
require(stats)
require(splines)
library(glmnet)
library(hdi)
library(foreach)
library(doParallel)

source("getZ.R")
source("lasso_inf.R")
source("dantzig.R")
source("simrun.R")

gBIC<-function(fitlasso,y){
  bic<-c()
  m=ncol(fitlasso$fitted)
  n=nrow(fitlasso$fitted)
  con<-log(n)/n
  for (i in 1:m){
    dff<-length(which(fitlasso$coefficients[,i]!=0))
    bbic<-log(sum((fitlasso$fitted[,i]-y)^2))+dff*con
    if(dff>0) bic<-c(bic,bbic)
    else bic<-c(bic,Inf)
  }
  return(clambda=which.min(bic))
}

data_generate4<-function(r=0.8,n=100,N,p)
{
  row1<-rep(0,p)
  for(i in 1:p){
    row1[i]<-r^(i-1)
  }
  Sigma1<-circulant.spam(row1)
  Sigma1[lower.tri(Sigma1)]<-0
  Sigma<-Sigma1+t(Sigma1)-diag(diag(Sigma1))
  X<-mvrnorm((n+N),rep(0,p),Sigma)
  X[,1]<-abs(X[,1])
  y<-c()
  for (i in 1:(n+N))
  { yy<-0.5*X[i,1]^2+4*X[i,3]^3/5-(X[i,4]-2)^2+2*(X[i,5]+1)^2+2*X[i,6]+rnorm(1,sd=1)
    #yy<-0.6*(X[i,1]+X[i,2])^2+2*X[i,4]^3/5-X[i,5]+2*X[i,6]+rnorm(1,sd=1)
    #yy<-0.5*X[i,1]^2+4*X[i,3]^3/5-(X[i,4]-2)^2+2*(X[i,5]+1)^2+2*X[i,6]+sum(2*X[i,supp_t2])+0.05*(sum(X[i,supp_t1]))^3+rnorm(1,sd=1)
    y<-c(y,yy)
  }  
  return(data=list(y=y,X=X))
}


#beta_t[1:6]<-c(1.45,1.04,0,1.2,-1,2)
#supp_t1=sample(c(10:(p/2)),5)
#supp_t2=sample(c((p/2+1):p),5)
#beta_t[supp_t1]=rep(0.75,5)
#beta_t[supp_t2]=rep(2,5)
