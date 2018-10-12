#sDiagnosis codes', 0- carriage, 1- pneumonia, 2- OM, AOM, 3- sinusitis, 4- conjuctivitis (1 strain ceratitis), 5- sepsis, 6-abscess
#Analerobic codes: Ambient+ Catalase is 2, oxyrase is 1, ambient air=0

#Clear workspace to get rid of old junk
rm(list=ls(all=TRUE))
#setwd("growth curve analyses\\")
library(grofit)
library(reshape)
library(RColorBrewer)
library(dplyr)
library(fdapace)
library(EMCluster)
library( aplpack )
library(ks)
library(nlme)


##IMPORT AND FORMAT DATA
d1a<- read.csv("master 5_24_2018.csv")
d1a$st[which(d1a$st=="11")]<-"11A"
d1a$st[which(d1a$st=="15")]<-"15B"
d1a$st[d1a$st=='5' &d1a$Diagnosis==5]<-"14"  #Our ST 5 strain from CDC is actually ST 14
#d1a$dx<-NA
#d1a$dx[d1a$Diagnosis ==0]<-0
#d1a$dx[d1a$Diagnosis %in% c(5)]<-5
#d1a$dx[d1a$Diagnosis %in% c(2)]<-2
d1a$st=as.character(d1a$st)
d1a$Diagnosis=as.numeric(d1a$Diagnosis)
d1a$X<-NULL
d1a$Age<-NULL
d1a$Agegroup=NULL
d1a$Ery<-NULL
d1a$Pen<-NULL
d1a$Pathogen<-NULL
d1a$Pilus<-NULL
d1a$Variant<-NULL
#index<-1:nrow(d1a)
ID2<-paste(d1a$Group,d1a$st,d1a$Diagnosis,sep="_")
d1a<-cbind.data.frame(ID2, d1a)
d1a$ID<-as.character(d1a$ID)
d1a$ID[d1a$ID==""]<-as.character(d1a$ID2[d1a$ID==""])
d1a$ID<-as.factor(d1a$ID)
d1a$Group<-NULL
d1a$Diagnosis[substr(d1a$st,1,3)=="603"]<-7 #isogenic
d1a$Diagnosis[substr(d1a$st,1,4)=="TIGR"]<-8 #isogenic
d1a$Diagnosis[d1a$st=="15B-" | d1a$st=='10A-']<-9 #cps- knockout
d1a<-d1a[substr(d1a$st,1,1) != "P",]
d1a<-d1a[d1a$st != "SAUR",]
d1a<-d1a[d1a$st != "Morax",]
d1a<-d1a[d1a$st != "Blank",]
d1a<-d1a[d1a$st != "",]
d1a$st[d1a$st %in% c("15B", "15C", '15B/C')] <-"15BC"
unique(d1a$st)
## blank based on t=0
d1a$X1[d1a$X1>d1a$X2]<- d1a$X2[d1a$X1>d1a$X2]
d1a[,7:ncol(d1a)]<-d1a[,7:ncol(d1a)] - (d1a$X1) #
d1a[,7:ncol(d1a)][d1a[,7:ncol(d1a)]<0]<-0  #Values must be >=0
d1a<-d1a[d1a$temp>=30 & d1a$temp<=39,]
d1a<-d1a[d1a$st !="unknown",]
d1a<-d1a[d1a$st !="GBS",]
d1a<-d1a[which(!is.na(d1a$Diagnosis)),]
d1a$st[d1a$st=="NT"]<-"cps-"
d1a$st[d1a$st=="TIGR4 Janus"]<-"cps-"
d1a$st[d1a$st=="603 Janus"]<-"cps-"
d1a$st[d1a$st=="603"]<-"6B"
d1a$st[d1a$st=="TIGR4"]<-"4"
d1a$st[d1a$st=="TIGR5"]<-"5"
d1a$st[d1a$st=="TIGR14"]<-"14"
d1a$st[d1a$st=="TIGR19F"]<-"19F"
d1a$st[d1a$st=="R6"]<-"cps-"
d1a<-d1a[d1a$temp>=30 & d1a$temp<=39,] #Restrict to 30-39C temperatures (readings at lower temps not reliable)

#identify which STs have isolates for both carriage AND IPD??
st.dx.table<- as.data.frame.matrix(table(d1a$st, d1a$Diagnosis))
st.keep.dx<- row.names(st.dx.table)[st.dx.table$'5'!=0 & st.dx.table$'0'!=0] 

##Function PCA
# Principal Modes of Variation for Processes With Continuous Sample Curves Castro, Lawton, SYLVESTRE
##strightforwrd review: https://www.tandfonline.com/doi/full/10.1080/14763141.2017.1392594?src=recsys
str(d1a)
#Y.list<- split(log(d1a[,-c(1:8, (48+8+1):ncol(d1a))] +0.01), as.factor(1:nrow(d1a)))
Y.list<- split(d1a[,-c(1:8, (48+8+1):ncol(d1a))], as.factor(1:nrow(d1a)))
Y.list<-lapply(Y.list, function(x) as.numeric(x[1,])) #List of OD600 values
T.list<-lapply(Y.list, function(x) (1:length(x))/2 ) #List of times (h)
fpca.growth<-FPCA(Y.list, T.list, optns= list(dataType='Dense', plot=TRUE))
plot(fpca.growth)
pct.var.pca<-fpca.growth$lambda/sum(fpca.growth$lambda)*100 #percent of variation  associated with each component
matplot(fpca.growth$phi[,c(1:3)]) #First 3 eigenfunctions

fit1<-fitted(fpca.growth, K=4)
#Plot a random selection of 10
matplot(t(fit1[1:10,]), type='l')

#Cluster the curves based on the FPCA
A <- FClust(Y.list,T.list, optnsFPCA = list(dataType='DenseWithMV'), k = 4)
CreatePathPlot( fpca.growth, K=4, showObs=FALSE, lty=1, col= A$cluster, xlab = 'Time (h)', ylab = 'OD600')
clusters<-A$cluster
gc.labels<-d1a[,1:8]
gc.labels$cluster<-clusters

#Extract fitted, first, and second derivatives
fitted.all<-fitted(fpca.growth,K = 3,derOptns = list(p = 0, bw = 1.01 , kernelType = 'epan') )
derivs.all<-fitted(fpca.growth,K = 3,derOptns = list(p = 1, bw = 1.01 , kernelType = 'epan') )
#2nd deriv
derivs2.all<-fitted(fpca.growth,K = 3,derOptns = list(p = 2, bw = 1.01 , kernelType = 'epan') )

plot.select<-c(ana.select=2, temp.select=37,diag.select=5)
par(mfrow=c(1,3))
matplot(t(fitted.all[d1a$anaerobic==plot.select['ana.select'] & d1a$temp==plot.select['temp.select'] & d1a$Diagnosis==plot.select['diag.select'],]), type='l',col=trans.black, lty=1, bty='l')
title("Smooth Observed")
matplot(t(derivs.all[d1a$anaerobic==plot.select['ana.select'] & d1a$temp==plot.select['temp.select'] & d1a$Diagnosis==plot.select['diag.select'],]), type='l',col=trans.black, lty=1, bty='l')
  abline(h=0, col='red')
  title("1st Derivative")
matplot(t(derivs2.all[d1a$anaerobic==plot.select['ana.select'] & d1a$temp==plot.select['temp.select'] & d1a$Diagnosis==plot.select['diag.select'],]), type='l',col=trans.black, lty=1, bty='l')
  abline(h=0, col='red')
title("2nd Derivative")
################################################


##################################################
#Try to pull out time series that increase, then flatten, then increase more
#Which have 1st derivs that stay >0 (never decline)
max.deriv1<-apply(derivs.all[,c(1:30)],1,function(x) max(x, na.rm=TRUE))
max.deriv2<-apply(derivs2.all[,c(1:30)],1,function(x) max(x, na.rm=TRUE))
min.deriv1<-apply(derivs.all[,c(1:30)],1,function(x) min(x, na.rm=TRUE))
min.deriv2<-apply(derivs2.all[,c(1:30)],1,function(x) min(x, na.rm=TRUE))
deriv1.start.increase<-apply(derivs.all[,c(1:30)],1,function(x){
  which(lag(x)<0 & x>0)
} )
min.deriv1.0<-min.deriv1
min.deriv1.0[min.deriv1>0]<-0 #if never declines, set minimum of deriv1 to 0
#Ratio of max rate of increase vs max rate of decrease 9if decreases
ratio_min_max_rate<-abs(min.deriv1.0)/max.deriv1 #Bigger values=faster decline relative to increase

#What is time of max 1st derivative?
t.max.deriv1<-rep(NA, length=length(max.deriv1))
for(i in 1:length(max.deriv1)){
  t.max.deriv1[i]<-which(derivs.all[i,c(1:30)]==max.deriv1[i])[1]
}
#What is time of min 2nd derivative?
t.min.deriv2<-rep(NA, length=length(min.deriv2))
for(i in 1:length(min.deriv2)){
  t.min.deriv2[i]<-which(derivs2.all[i,c(1:30)]==min.deriv2[i])[1]
}
#Plot based on 
curve.cuts<-list(c(0,0.05),c(0.05,0.2),c(0.2,10))
par(mfrow=c(length(curve.cuts),3))
for(i in 1:length(curve.cuts)){
increase.indices<-which(ratio_min_max_rate<curve.cuts[[i]][2] & ratio_min_max_rate>=curve.cuts[[i]][1]) #What is cutoff--lower values give ones that decrease less
matplot(t(fitted.all[increase.indices,]), type='l',col=trans.black, lty=1, bty='l')
title(paste0("Smooth Observed ",i))
matplot(t(derivs.all[increase.indices,]), type='l',col=trans.black, lty=1, bty='l')
abline(h=0, col='red')
title("1st Derivative")
matplot(t(derivs2.all[increase.indices,]), type='l',col=trans.black, lty=1, bty='l')
abline(h=0, col='red')
title("2nd Derivative")
increasing.samples<-d1a[increase.indices, c(1:6)]
increasing.samples<-increasing.samples[order(increasing.samples$st),]
increasing.samples
}

#Predictors of having bigger ratio
mod.df<-cbind.data.frame(ratio_min_max_rate=ratio_min_max_rate,id=d1a$ID, st=d1a$st, anaerobic=as.factor(d1a$anaerobic), t.max.deriv1=t.max.deriv1, temp= as.factor(d1a$temp))
mod1<-lme(sqrt(ratio_min_max_rate) ~ st+anaerobic +t.max.deriv1 + temp , random=~1|id, data=mod.df[mod.df$anaerobic !='0',] )  
summary(mod1)
coef1<-mod1$coefficients$fixed
coef.st<-coef1[substr(names(coef1),1,2)=='st']
sort(coef.st)

#Focus only on ST3 and 8
mod.df<-cbind.data.frame(ratio_min_max_rate=ratio_min_max_rate,id=d1a$ID, st=d1a$st, anaerobic=as.factor(d1a$anaerobic), t.max.deriv1=t.max.deriv1, temp= as.factor(d1a$temp))
mod2<-lme(sqrt(ratio_min_max_rate) ~ st+anaerobic +t.max.deriv1 + temp , random=~1|id,
          data=mod.df[mod.df$anaerobic !='0' & mod.df$ st %in% c('3','8'),] )  
summary(mod2)
coef1<-mod1$coefficients$fixed
coef.st<-coef1[substr(names(coef1),1,2)=='st']
sort(coef.st)


#Look at ST3 only: at 38C, see the phenotype in some strains temps, see an
#Seen with both analerobic and aerobic; also see it at 30,33; not really 35,37
subset1<- (d1a$st=='8' &  d1a$anaerobic!=0 &d1a$anaerobic==2 &d1a$temp==35)
par(mfrow=c(1,3))
matplot(t(fitted.all[subset1,]), type='l',col=trans.black, lty=1, bty='l')
title("Smooth Observed")
matplot(t(derivs.all[subset1,]), type='l',col=trans.black, lty=1, bty='l')
abline(h=0, col='red')
title("1st Derivative")
matplot(t(derivs2.all[subset1,]), type='l',col=trans.black, lty=1, bty='l')
abline(h=0, col='red')
title("2nd Derivative")


#matplot(t(fitted.all[d1a$ID2=='CDC_3_5' &d1a$temp==33 &d1a$Diagnosis==5,]), type='l',col=d1a$anaerobic+1, lty=1, bty='l')


#Select subset based on PC3 variation:
#Select subset based on PC3 variation:
#3rd PC seems to show evidence of carbon switching:
pc3<-fpcs$PC3
pc3.select<- sort(pc3)[c(50, round(length(pc3)/2), length(pc3)-50)]
pc3.select.index<- which(pc3 %in% pc3.select)
CreatePathPlot(fpca.growth, subset = c(pc3.select.index), K = 3, main = 'K = 3', showObs = FALSE) ; grid()
CreatePathPlot(fpca.growth, subset = c(pc3.select.index), K = 3, main = 'K = 3: Derivative', showObs = FALSE, derOptns = list(p = 1, bw = 1.01 , kernelType = 'epan') ) ; grid()





#####################
#####################
#VISUALIZE HOW GROWTH CURVES VARY ALONG PC3
plot_func<-function(fpcaObj, plot.quants){
  #Eigenfunctions
  eigenfunc<-fpcaObj$phi  #Eigenfunction
  xiEst2<-fpcaObj$xiEst     
  est.meanfunc<-fpcaObj$mu  #Mean function
  xi.est2.quantiles<- apply(xiEst2,2, quantile, probs=plot.quants)
  yEst.q<-array(NA, dim=c(length(est.meanfunc),nrow(xiEst2),ncol(xiEst2) ))
  #Just reconstruct based on PCN
  for(pcn in 1:ncol(xiEst2) ){
    yEst.q[,,pcn]<- t(xiEst2[,pcn, drop=FALSE] %*%  t(eigenfunc[,pcn]) + t(matrix(rep(est.meanfunc,nrow(xiEst2)), nrow=length(est.meanfunc))) )
  }
  return(yEst.q)
}
pca.all<-plot_func(fpca.growth,plot.quants=c(0.1,0.5,0.9))
pc3<-fpcs$PC3
pc2<-fpcs$PC2

pc2.q<-quantile(pc2, probs=c(0.05, 0.1,0.2,0.4,0.6,0.8))
pc3.q<-quantile(pc3, probs=c(0.05, 0.1,0.2,0.4,0.6,0.8))
pc3.grp1<- which(pc3<pc3.q[1])
pc3.grp2<-which(pc3>=pc3.q[1] & pc3<pc3.q[2])
pc3.grp3<-which(pc3>=pc3.q[2] & pc3<pc3.q[3])
pc3.grp4<-which(pc3>=pc3.q[3] & pc3<pc3.q[4])
pc3.grp5<-which(pc3>=pc3.q[4] & pc3<pc3.q[5])
pc3.grp6<-which(pc3>=pc3.q[5] & pc3<pc3.q[6])
pc3.grp7<-which(pc3>=pc3.q[6])
#PC2 score is high
pc2.grp6<-which(pc2>=pc2.q[5] & pc2<pc2.q[6])
pc2.grp7<-which(pc2>=pc2.q[6])

#Plot this variability
par(mfrow=c(1,7))
matplot(pca.all[,pc3.grp1,2], type='l', col=trans.black, bty='l')
matplot(pca.all[,pc3.grp2,2], type='l', col=trans.black, bty='l')
matplot(pca.all[,pc3.grp3,2], type='l', col=trans.black, bty='l')
matplot(pca.all[,pc3.grp4,2], type='l', col=trans.black, bty='l')
matplot(pca.all[,pc3.grp5,2], type='l', col=trans.black, bty='l')
matplot(pca.all[,pc3.grp6,2], type='l', col=trans.black, bty='l')
matplot(pca.all[,pc3.grp7,2], type='l', col=trans.black, bty='l')

# #Correlates of PC3
# #Low on PC3 and High on PC2 should give desire phenotype?
# Reduce(intersect, list(c(pc3.grp1,pc3.grp2),c(pc2.grp6,pc2.grp7)))
#Observed time series for low PC3 scores--don't bear much resemblance to PC plots
par(mfrow=c(1,1))
matplot(t(d1a[pc3.grp1, -c(1:6)]), col=trans.black, type='l')
matplot(t(d1a[increase.indices, -c(1:6)]), col=trans.black, type='l')




#Or first take mean of PCs by ST, then do regression
sub.fpc2<-fpc2[fpc2$Diagnosis %in% c('0','5') & fpc2$anaerobic %in% c('1','2') , c('st','Diagnosis','PC1','PC2','PC3','PC4')]
mean.fpcs<-aggregate( cbind(sub.fpc2$PC1, sub.fpc2$PC2, sub.fpc2$PC3, sub.fpc2$PC4)~sub.fpc2$st+ sub.fpc2$Diagnosis, FUN=mean)
names(mean.fpcs)<-names(sub.fpc2)
mean.fpcs2<- merge(cov1, mean.fpcs, by='st')
#fpc2<-mean.fpcs2[fpc2$Diagnosis %in% c(0,5),]
reg.mean.pc<-glm(sqrt(mean.fpcs2$ave_carr_pre)~ PC1+PC2+Diagnosis, data=mean.fpcs2)
summary(reg.mean.pc)

reg.mean.pc<-glm(sqrt(mean.fpcs2$ave_carr_pre)~ PC2+Diagnosis + PC2*Diagnosis, data=mean.fpcs2)
summary(reg.mean.pc)

pred.mean.pc<-predict(reg.mean.pc)
plot(pred.mean.pc, mean.fpcs2$ave_carr_pre, col='white', bty='l')
text(pred.mean.pc, mean.fpcs2$ave_carr_pre, mean.fpcs2$st, col=as.numeric(factor(mean.fpcs2$Diagnosis)))

#Regression where PC1 or PC2 are the outcome:
#carriage prev associated with PC2 (timing) but not PC1 (intenisty)
fpc3<-merge(fpc2,ps1,by='st', all=TRUE)
fpc3$NAC<-0
fpc3$NAC[fpc3$GlcNAc==1 | fpc3$GalNAc==1 | fpc3$ManNAc==1 |  fpc3$ManNAcA==1 |fpc3$FucNAc==1 |fpc3$PneNAc ==1] <-1
fpc3$uronic<-0
fpc3$uronic[fpc3$GalA==1 |fpc3$GlcA==1 ]<-1  #GlcA and GalA derived from same pathway

#3rd PC seems to show evidence of carbon switching--low values of PC3  have a carbon switching phenotype (plateau then increase, high values peak then decline)
pc3<-fpcs$PC3
plot(as.factor(fpcs$Diagnosis) , fpcs$PC3)
plot(as.factor(fpcs$anaerobic) , fpcs$PC3)
plot(as.factor(fpcs$st) , fpcs$PC3)
fpcs.model<-fpcs
fpcs.model$Diagnosis<-as.factor(fpcs.model$Diagnosis)
fpcs.model$anaerobic<-as.factor(fpcs.model$anaerobic)
fpcs.model$temp<-as.factor(fpcs.model$temp)
fpcs.model$st<-as.factor(fpcs.model$st)
mod.pc3<-lm(PC3~ st+Diagnosis+anaerobic +temp, data=fpcs.model)
summary(mod.pc3)

mod.pc3b<-lm(PC3~ st+anaerobic , data=fpcs.model[fpcs.model$temp=='33' &fpcs.model$Diagnosis=='0',])
summary(mod.pc3b)
mod.pc3c<-lm(PC3~ st+anaerobic , data=fpcs.model[fpcs.model$temp=='37' &fpcs.model$Diagnosis=='5',])
summary(mod.pc3c)

#EFFECT OF OXYGEN BY SEROTYPE FOR PC! AND 2
reg1<-lme(PC1 ~ ST + anaerobic +temp +Diagnosis
          +ST*anaerobic
          , random=~ 1|ID, data=fpc3[fpc3$anaerobic %in% c('1','2'),], na.action=na.omit )
summary(reg1)
reg2<-lme(PC2 ~ ST + anaerobic +temp +Diagnosis
          +ST*anaerobic
          , random=~ 1|ID, data=fpc3[fpc3$anaerobic %in% c('1,2'),], na.action=na.omit )
summary(reg2)

reg1<-lme(PC1 ~ ST + anaerobic +temp +Diagnosis, random=~ 1|ID, data=fpc3[fpc3$Diagnosis %in% c('0','5') & fpc3$anaerobic %in% c('1','2'),], na.action=na.omit )
summary(reg1)
reg1a<-lme(PC1 ~ sqrt(ave_carr_pre) + anaerobic +temp +Diagnosis, random=~ 1|ID, data=fpc3[fpc3$Diagnosis %in% c('0','5') & fpc3$anaerobic %in% c('1','2'),], na.action=na.omit )
summary(reg1a)
reg2<-lme(PC2 ~ ST + anaerobic +temp +Diagnosis, random=~ 1|ID, data=fpc3[fpc3$Diagnosis %in% c('0','5') & fpc3$anaerobic %in% c('1','2'),], na.action=na.omit )
summary(reg2)
reg2a<-lme(PC2 ~ sqrt(ave_carr_pre) + anaerobic +temp +Diagnosis, random=~ 1|ID, data=fpc3[fpc3$Diagnosis %in% c('0','5') & fpc3$anaerobic %in% c('1','2'),], na.action=na.omit )
summary(reg2a)

#cARBON ASSOCIATED WITH PC1(EHIGH) NOT pc2 (TIMING)
reg4<-lme(PC1 ~ NAC + anaerobic +temp +Diagnosis,random=~ 1|ID, data=fpc3[fpc3$Diagnosis %in% c('0','5') & fpc3$anaerobic %in% c('1','2'),], na.action=na.omit )
summary(reg4)
reg4<-lme(PC2 ~ NAC + anaerobic +temp +Diagnosis, random=~ 1|ID, data=fpc3[fpc3$Diagnosis %in% c('0','5') & fpc3$anaerobic %in% c('1','2'),], na.action=na.omit )
summary(reg4)
reg5<-lme(PC1 ~ uronic + anaerobic +temp +Diagnosis, random=~ 1|ID, data=fpc3[fpc3$Diagnosis %in% c('0','5') & fpc3$anaerobic %in% c('1','2'),], na.action=na.omit )
summary(reg5)
reg6<-lme(PC2 ~ uronic + anaerobic +temp +Diagnosis, random=~ 1|ID, data=fpc3[fpc3$Diagnosis %in% c('0','5') & fpc3$anaerobic %in% c('1','2'),], na.action=na.omit )
summary(reg6)



#####
#FUNCTION TO EXTRACT fitted value based on 1st or second PCs
# based on " Displaying the Important Features of Large Collections of Similar Curves" Jones and Rice 1989
plot_func<-function(fpcaObj, plot.quants){
  #Eigenfunctions
  eigenfunc<-fpcaObj$phi  #Eigenfunction
  xiEst2<-fpcaObj$xiEst     
  est.meanfunc<-fpcaObj$mu  #Mean function
  xi.est2.quantiles<- apply(xiEst2,2, quantile, probs=plot.quants)
  yEst.q<-array(NA, dim=c(length(est.meanfunc),3,ncol(xi.est2.quantiles) ))
  #Just reconstruct based on PCN
  for(pcn in 1:ncol(xi.est2.quantiles) ){
    yEst.q[,,pcn]<- t(xi.est2.quantiles[,pcn] %*%  t(eigenfunc[,pcn]) + t(matrix(rep(est.meanfunc,3), nrow=length(est.meanfunc))) )
  }
  return(yEst.q)
}
quants.pca<-plot_func(fpca.growth,plot.quants=c(0.1,0.5,0.9))
par(mfrow=c(2,2))
for(i in 1:4){
matplot(quants.pca[,,i], type='l', bty='l', col=c('#1b9e77','#d95f02','#7570b3' ), lty=c(2,1,3), lwd=c(1,2,1)) #Plot 1st principal component range
  title(paste0("Variation associated with PCA ", i,"; ", round(pct.var.pca[i],0),"%"))
}
###


###Anne suggestion:
#PC1 vs PC2--overlay st, temp, etc info
par(mfrow=c(1,2))
subs<-fpcs[fpcs$anaerobic %in% c('1','2') ,]
subs$temp<-relevel(subs$temp,ref='30')
plot(subs$PC1, subs$PC2, col=subs$anaerobic , bty='l')
plot(subs$PC1, subs$PC2, col=ts.col[subs$temp ], bty='l')
#plot(subs$PC1, subs$PC2, col=ts.col[subs$st ])

subs.cl<-scale(cbind(subs$PC1, subs$PC2))
fit <- kmeans(subs.cl, 3) # 5 cluster solution
plot(subs$PC1, subs$PC2, col=fit$cluster, bty='l')
abline(v=subs$PC1[subs$PC2==min(subs$PC2)])

#quadrants based on whether above or below median for PC1 and PC2
subs$grp<-1
subs$grp[subs$PC1>=median(subs$PC1) & subs$PC2<median(subs$PC2)]<-1
subs$grp[subs$PC1>=median(subs$PC1) & subs$PC2>=median(subs$PC2)]<-2
subs$grp[subs$PC1<median(subs$PC1) & subs$PC2>=median(subs$PC2)]<-3
subs$grp[subs$PC1<median(subs$PC1) & subs$PC2<median(subs$PC2)]<-4
plot(subs$PC1, subs$PC2, col=subs$grp, bty='l')
tab1<-table(subs$st,subs$grp )
 prop.table(tab1,1)
 library(dendextend)
 d <- dist(subs[,c('PC1','PC2')], method = "euclidean") # distance matrix
 fit <- hclust(d, method="ward.D") 
 hcl.groups <- cutree(fit, k=4) # cut tree into 4 clusters
 hcd <- as.dendrogram(fit)
 labels(hcd)<-paste(subs$st,subs$anaerobic, subs$temp, subs$Diagnosis,  sep="_")
 par(mfrow=c(1,1), mar=c(6,1,1,1))
 plot(hcd, xlim=c(1,40), ylim=c(1,30)) # display dendogram
 
 #HCL GROUPINGS
 tab2<-table(subs$st,groups )
 prop.table(tab2,1)
 plot(subs$PC1, subs$PC2, col='white', bty='l')
 text(subs$PC1, subs$PC2,hcl.groups, col=hcl.groups)
  tab2<-table(subs$temp,groups )
 prop.table(tab2,1)
 
 #distance matrix--after controling for shared temp or O2 level, are some ST pairs found together more frequently on PC1, PC2?#
 #USE INVERSE DISTANCE IN OUTCOME OF REGRESSION, CONTROL FOR WHETHER SHARE CONDITIONS< THEN RAND EFECT FOR PAIRING...
#length of d2 is (nrow(subs)^2-nrow(subs))/2
 d2 <- as.matrix(dist(subs[,c('PC1','PC2')], method = "euclidean", diag=FALSE, upper=FALSE)) # distance matrix
 dimnames(d2)[[1]]<-paste(subs$st,subs$anaerobic, subs$temp, subs$Diagnosis, 1:nrow(subs), sep="_")
 dimnames(d2)[[2]]<-paste(subs$st,subs$anaerobic, subs$temp, subs$Diagnosis, 1:nrow(subs), sep="_")
 d3<-melt(d2)  #pairwise distances
 d3<-d3[d3$Var1 != d3$Var2,]
 lab.var1<-as.data.frame(matrix(unlist(strsplit( as.character(d3$Var1), "_", fixed = TRUE)), ncol=5, byrow=TRUE))
 names(lab.var1)<-c('st.var1','anaerobic.var1', 'temp.var1', 'Diagnosis.var1', 'index.var1')
 lab.var2<-as.data.frame(matrix(unlist(strsplit( as.character(d3$Var2), "_", fixed = TRUE)), ncol=5, byrow=TRUE))
 names(lab.var2)<-c('st.var2','anaerobic.var2', 'temp.var2', 'Diagnosis.var2', 'index.var2')
 d4<-cbind.data.frame(d3, lab.var1, lab.var2)
 d4$same.st<-0
 d4$same.st[d4$st.var1==d4$st.var2]<-1
 d4$same.o2<-0
 d4$same.o2[d4$anaerobic.var1==d4$anaerobic.var2]<-1
 d4$same.temp<-0
 d4$same.temp[d4$temp.var1==d4$temp.var2]<-1
 d4$same.diag<-0
 d4$same.diag[d4$Diagnosis.var1==d4$Diagnosis.var2]<-1
 d4$inv.distance<-1/d4$value
 d4$st.pair<-as.factor(paste0(d4$st.var1, '_',d4$st.var2))
 mod1<-lme(inv.distance~ same.o2 + same.temp +same.st +same.diag ,random=~1|st.pair, data=d4[!is.infinite(d4$inv.distance),], na.action=na.exclude)
 re.coefs<-summary(mod1)$coefficients$random[[1]][,1]
 re.coefs<-sort(re.coefs, decreasing=TRUE)
 #Higher values==smaller distance in pc1, PC2 after controling for covariates
 re.coefs[1:1000]
 re.coefs.labs<- matrix(unlist(strsplit( as.character(names(re.coefs)), "_", fixed = TRUE)), ncol=2, byrow=TRUE)
 re.coefs2<-cbind.data.frame(re.coefs.labs, re.coefs)
 names(re.coefs2)<-c('stA','stB', 'dist')
 re.coefs.m<-melt(re.coefs2, id.vars=c('stA','stB'))
re.coefs.c<-as.data.frame(cast(re.coefs.m, stA~stB))
re.coef.mat<-as.matrix(re.coefs.c[,-1])
dimnames(re.coef.mat)[[1]]<-re.coefs.c$stA
dimnames(re.coef.mat)[[2]]<-re.coefs.c$stA
hmcols<-colorRampPalette(c("white","blue"))(256)
hm1<-heatmap(re.coef.mat, scale='none', col=hmcols)
hm1<-heatmap(re.coef.mat, Rowv=NA, Colv=NA,scale='none' , revC = TRUE, col=hmcols)
title('Inverse Distance between serotypes based on PC1 and PC2')

dist<- dist(-re.coef.mat) #re is for model fit to inverse distance, so negative is distance?
 hc <- hclust(dist)                # apply hirarchical clustering 
plot(hc)
title("Serotypes Relatedness of variation along PC1 and PC2 after adjusting for culture and strain characteristics")

