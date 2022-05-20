#### Change the color in p_ab, p_sq for estimation and in line 60 for inference.
#### Change the file name in data source to plot results of different setting. 
#### The code for inference parts does not work for model1. 
#### Plots for estimation results can be generated using the same code.
############################################################################################
#### Estimation

source("data source.R")
plotdata<-rbind(esresult[c(2,8),c(1:2)],esresult1[c(2,6),c(1:2)],
                esresult2[c(2,6),c(1:2)],esresult[6,c(1:2)],esresult[6,c(1:2)],esresult[6,c(1:2)])
method=c(rep(c("EFF","SAFE"),3),rep("NAIVE",3))
nNratio<-c(rep(1,2),rep(4,2),rep(8,2),1,4,8)
colnames(plotdata)<- c("square error","sd")
ppd<-data.frame(method,nNratio,plotdata)
ppd$sd<-ppd$sd/10

p_ab<-ggplot(ppd, aes(x=nNratio,y=square.error, colour=method, group=method)) + 
  geom_errorbar(aes(ymin=square.error-sd, ymax=square.error+sd), width=.1) +
  geom_line() +
  geom_point()+
  theme(legend.position = "none")+
  xlab("N/n ratio") +
  ylab("absolute error")


plotdata<-rbind(esresult[c(1,7),c(1:2)],esresult1[c(1,5),c(1:2)],
                esresult2[c(1,5),c(1:2)],esresult[5,c(1:2)],esresult[5,c(1:2)],esresult[5,c(1:2)])
method=c(rep(c("EFF","SAFE"),3),rep("NAIVE",3))
nNratio<-c(rep(1,2),rep(4,2),rep(8,2),1,4,8)
colnames(plotdata)<- c("square error","sd")
ppd<-data.frame(method,nNratio,plotdata)
ppd$sd<-ppd$sd/10

p_sq<-ggplot(ppd, aes(x=nNratio,y=square.error, colour=method, group=method)) + 
  geom_errorbar(aes(ymin=square.error-sd, ymax=square.error+sd), width=.1) +
  geom_line() +
  geom_point()+
  theme(legend.position = "none")+
  xlab("N/n ratio") +
  ylab("squared error") 

figure1<-ggarrange(p_sq,p_ab,common.legend = T,
                   ncol = 2, nrow = 1)

annotate_figure(figure1,top = text_grob("Estimation Results for Model 2 (p=500,n=200)", color = "black", face = "bold", size = 12))

########################################################################################################
# Plots for Inference.
# To draw the plots of x1, x2, x5, x6 used as the input of line 69 ggarrange function,
# change m to be 1,2,5,6 in line 54 and the names of variables x6 and coverage6 respectively in line 60.
# For example, to get subplot x1, change m=1, x6 to x1, and the label from coverage6 to coverage1.
#######################################################################################################

m=6
plotdata<-c(infresult1[c(5,1),m],infresult3[c(4,1),m],infresult5[c(4,1),m],infresult1[3,m],infresult1[3,m],infresult1[3,m])
nNratio<-c(rep(1,2),rep(4,2),rep(8,2),1,4,8)
coverage=c(infresult[c(5,1),m],infresult2[c(4,1),m],infresult4[c(4,1),m],infresult[3,m],infresult[3,m],infresult[3,m])
ppd<-data.frame(method,nNratio,plotdata)
assign(paste("coverage", m, sep=""),as.factor(round(coverage,2)))
x6<-ggplot(ppd, aes(x=nNratio,y=plotdata, colour=method,label= coverage6, group=method)) + 
  geom_line() +
  geom_point(size=1)+
  # geom_text(size=3.3,hjust = 0,vjust=0,nudge_x=-0.1,nudge_y = c(rep(c(-0.02,0.01,-0.02),3),rep(0.01,3)))+
  geom_text_repel(size=3.3)+
  #nudge_y = c(rep(c(0.01,0.02,-0.01),3),rep(0.01,3)))+
  xlab("N/n ratio") +
  ylab("length")+theme(legend.position = "none")

figure<-ggarrange(x6, x1, x2, x5, 
                  labels = c("X0","X1", "X2", "X6"),common.legend = TRUE,
                  ncol = 2, nrow = 2)

annotate_figure(figure,
                top = text_grob("Inference Results for Model 2 (p=500,n=200)", color = "black", face = "bold", size = 12))





