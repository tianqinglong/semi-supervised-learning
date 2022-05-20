
result <- read.csv("inf_M2p5n2.csv")
result = result[,-1]

esresult<-matrix(0,nrow=8,ncol=4)
for(i in 1:8)
{  
  er1<-result[i+(0:rep)*52,]
  err1<-rowSums(er1) ##l2 or l1
  sqe1<-mean(err1) 
  sqsd1<-sd(err1)
  mer1<-apply(er1, 1,function(x) {max(x)})
  mmer1<-mean(mer1)
  sdm1<-sd(mer1)
  esresult[i,]<-c(sqe1,sqsd1,mmer1,sdm1)
}

infresult<-matrix(0,nrow=6,ncol=6)
for(i in 9:14)
{  
  er1<-result[i+(0:rep)*52,]
  err1<-colMeans(er1[,c(1,2,4:6)])
  err2<-mean(colMeans(er1[,-c(1,2,4:6)]))
  infresult[i-8,]<-c(err1,err2)
}


infresult1<-matrix(0,nrow=6,ncol=6)
for(i in 15:20)
{  
  er1<-result[i+(0:rep)*52,]
  err1<-colMeans(er1[,c(1,2,4:6)])
  err2<-mean(colMeans(er1[,-c(1,2,4:6)]))
  infresult1[i-14,]<-c(err1,err2)
}

esresult1<-matrix(0,nrow=6,ncol=4)
for(i in 21:26)
{  
  er1<-result[i+(0:rep)*52,]
  err1<-rowSums(er1) ##l2 or l1
  sqe1<-mean(err1) 
  sqsd1<-sd(err1)
  mer1<-apply(er1, 1,function(x) {max(x)})
  mmer1<-mean(mer1)
  sdm1<-sd(mer1)
  esresult1[i-20,]<-c(sqe1,sqsd1,mmer1,sdm1)
}

infresult2<-matrix(0,nrow=5,ncol=6)
for(i in 27:31)
{  
  er1<-result[i+(0:rep)*52,]
  err1<-colMeans(er1[,c(1,2,4:6)])
  err2<-mean(colMeans(er1[,-c(1,2,4:6)]))
  infresult2[i-26,]<-c(err1,err2)
}


infresult3<-matrix(0,nrow=5,ncol=6)
for(i in 32:36)
{  
  er1<-result[i+(0:rep)*52,]
  err1<-colMeans(er1[,c(1,2,4:6)])
  err2<-mean(colMeans(er1[,-c(1,2,4:6)]))
  infresult3[i-31,]<-c(err1,err2)
}

infresult4<-matrix(0,nrow=5,ncol=6)


esresult2<-matrix(0,nrow=6,ncol=4)
for(i in 37:42)
{  
  er1<-result[i+(0:rep)*52,]
  err1<-rowSums(er1) ##l2 or l1
  sqe1<-mean(err1) 
  sqsd1<-sd(err1)
  mer1<-apply(er1, 1,function(x) {max(x)})
  mmer1<-mean(mer1)
  sdm1<-sd(mer1)
  esresult2[i-36,]<-c(sqe1,sqsd1,mmer1,sdm1)
}

infresult4<-matrix(0,nrow=5,ncol=6)

for(i in 43:47)
{  
  er1<-result[i+(0:rep)*52,]
  err1<-colMeans(er1[,c(1,2,4:6)])
  err2<-mean(colMeans(er1[,-c(1,2,4:6)]))
  infresult4[i-42,]<-c(err1,err2)
}

infresult5<-matrix(0,nrow=5,ncol=6)
for(i in 48:52)
{  
  er1<-result[i+(0:rep)*52,]
  err1<-colMeans(er1[,c(1,2,4:6)])
  err2<-mean(colMeans(er1[,-c(1,2,4:6)]))
  infresult5[i-47,]<-c(err1,err2)
}

