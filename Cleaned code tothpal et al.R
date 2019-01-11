#################################################################################################
#Code to analyze growth curve data from Tothpal et al. manuscript
#Key for Diagnosis vairable: 0- carriage, 1- pneumonia, 2- OM, AOM, 3- sinusitis, 4- conjuctivitis (1 strain ceratitis), 5- sepsis, 6-abscess
#Key for anaerobic variable: #Ambient air+Catalase is 2, oxyrase (anaerobic) = 1, ambient air (no catalase)=0
###################################################################################################

#Clear workspace to get rid of old junk
rm(list=ls(all=TRUE))
setwd("C:\\Users\\dmw63\\Desktop\\My documents h\\LAB\\growth curve analyses\\")
library(grofit)
library(reshape)
library(corrplot)
library(RColorBrewer)
library(nlme)
library(phia)
library(data.table)
require(caret)

##IMPORT AND CLEAN DATA, get rid of unneeded variables; fix some serotype coding issues
d1a<- read.csv("master 5_24_2018.csv")
d1a$st[which(d1a$st=="11")]<-"11A"
d1a$st[which(d1a$st=="15")]<-"15B"
d1a$st[d1a$st=='5' &d1a$Diagnosis==5]<-"14"  #Our ST 5 strain from CDC is actually ST 14
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
d1a<-d1a[d1a$temp>=30 & d1a$temp<=39,] #Restrict analysis to 30-39C (higher and lower temps not reliable)
d1a<-d1a[d1a$st !="unknown",]
d1a<-d1a[d1a$st !="GBS",]
d1a<-d1a[which(!is.na(d1a$Diagnosis)),]

## blank based on t=0
d1a$X1[d1a$X1>d1a$X2]<- d1a$X2[d1a$X1>d1a$X2]
d1a[,7:ncol(d1a)]<-d1a[,7:ncol(d1a)] - (d1a$X1) #
d1a[,7:ncol(d1a)][d1a[,7:ncol(d1a)]<0]<-0  #Values must be >=0

#Fix coding for isogenic strains
d1a$st[d1a$st=="NT"]<-"cps-"
d1a$st[d1a$st=="TIGR4 Janus"]<-"cps-"
d1a$st[d1a$st=="603 Janus"]<-"cps-"
d1a$st[d1a$st=="603"]<-"6B"
d1a$st[d1a$st=="TIGR4"]<-"4"
d1a$st[d1a$st=="TIGR5"]<-"5"
d1a$st[d1a$st=="TIGR14"]<-"14"
d1a$st[d1a$st=="TIGR19F"]<-"19F"
d1a$st[d1a$st=="R6"]<-"cps-"

#identify which STs have isolates for both carriage AND IPD??
st.dx.table<- as.data.frame.matrix(table(d1a$st, d1a$Diagnosis))
st.keep.dx<- row.names(st.dx.table)[st.dx.table$'5'!=0 & st.dx.table$'0'!=0] 


#tabulate # replicates for each condition
d1.tab<-d1a[,c('ID','st','anaerobic','Diagnosis','temp' )]
tab1<-setDT(d1.tab)[,list(Count=.N) ,names(d1.tab)]

####Generate average curves for use in plotting (average over replictes of same isolate)
d1a.agg<-d1a
d1a.agg$max.od<-NULL
d1_agg2<-aggregate(d1a.agg[,-c(1:6)], list(d1a.agg$ID2,d1a.agg$st, d1a.agg$anaerobic, d1a.agg$temp, d1a.agg$ID,d1a.agg$Diagnosis), mean, na.rm=TRUE)
names(d1_agg2)[1:6]<-names(d1a.agg)[1:6]

##################################################################
#PLOTS FOR PAPER: Figure 1 and supplementary figures with isogenic and capsule KO strains
#########FIG 1:
####PLOTS FOR PAPER--Fig 1 
ts.col<- rev(brewer.pal(6,"RdYlBu")) #define color palette
####CHANGE OPTIONS HERE FOR PLOT SELECTION
st.select<-'19F'  #Which serotype do you want to plot?
id.select.19F<-c('H133','H12','CDC_19F','H111') #H111 is a pneumonia isolate, H133 is conjuncitivitis; 
####
sub1<-par(mfrow=c(1,1))
#names(d1_agg2)[1]<-'temp.lab'
d1a.agg$anaerobic<-as.numeric(as.character(d1a.agg$anaerobic))
d1a.agg$temp<-as.numeric(as.character(d1a.agg$temp))
d1a.agg$st<-as.character(d1a.agg$st)
tiff('plot1.tiff',width=11, height=5.3, units='in', res=400)
# par(mfcol=c(3,4),mar=c(2,3,1,1))
par(mfrow=c(1,3),mar=c(2,3,1,1))
# for(dx.select in c(0,4,1,5)){ #Body sites
for(dx.select in c(0)){ #Body sites
  for(ana.select in c(0,2,1)){  #Catalase, anaerobic, aertobic
    id.composite<-paste0(d1a.agg$ID,d1a.agg$st,d1a.agg$temp,d1a.agg$anaerobic)
    d1.df<-d1a.agg[d1a.agg$st==st.select & duplicated(id.composite)==FALSE & d1a.agg$anaerobic %in% c(ana.select) & d1a.agg$Diagnosis %in% c(dx.select) &d1a.agg$ID %in% id.select.19F, ]
    d1.df<-d1.df[ order(-d1.df$anaerobic, -d1.df$temp, d1.df$st), ]
    hm.df<-(as.matrix(d1.df[ , 7:52]))
    max.ts.time<-apply(t(hm.df),2,which.max)/2 -0.5 #time at which max is achieved for each TS
    max.ts.od<-apply(t(hm.df),2,max, na.rm=TRUE) #OD at  max  for each TS
    temp.lab<-as.numeric(as.factor(d1.df$temp))
    matplot(seq(from=0, to=(nrow(t(hm.df))-1)/2, by=0.5),t(hm.df),lwd=2,  type="l",lty=1,yaxt='n', col=ts.col[temp.lab],xlab="Time (h)",bty="n", ylab="OD",ylim=c(0,0.5) )
    if(dx.select=='0'){axis(side=2 ,at=c(0,0.1,0.2,0.3,0.4,0.5))    } #if yaxt='n'
    text(max.ts.time, max.ts.od+0.02, d1.df$temp, col=ts.col[temp.lab], cex=0.75 )
  }
}
dev.off()


#######SUPPLEMENTARY FIGURE PLOTTING TIGR ISOLATES
####CHANGE OPTIONS HERE FOR PLOT SELECTION
sub1<-par(mfrow=c(1,1))
#names(d1_agg2)[1]<-'temp.lab'
d1a.agg$anaerobic<-as.numeric(as.character(d1a.agg$anaerobic))
d1a.agg$temp<-as.numeric(as.character(d1a.agg$temp))
d1a.agg$st<-as.character(d1a.agg$st)
ts.col2<-ts.col[2:4]
#TIGR4 strains; stratified by strain and condition
tiff('TIGR cps switch by strain.tiff',width=11, height=5.3, units='in', res=400)
par(mfcol=c(3,4),mar=c(2,3,1,1))
for(st.select in c('5','14','19F','cps-')){ #Body sites
  for(ana.select in c(0,2,1)){  #Catalase, anaerobic, aertobic
    id.composite<-paste0(d1a.agg$ID,d1a.agg$st,d1a.agg$temp,d1a.agg$anaerobic)
    d1.df<-d1a.agg[ duplicated(id.composite)==FALSE & d1a.agg$anaerobic %in% c(ana.select) & d1a.agg$Diagnosis %in% c(8)  & d1a.agg$st %in% st.select, ]
    d1.df<-d1.df[ order(-d1.df$anaerobic, -d1.df$temp, d1.df$st), ]
    hm.df<-(as.matrix(d1.df[ , 7:52]))
    max.ts.time<-apply(t(hm.df),2,which.max)/2 -0.5 #time at which max is achieved for each TS
    max.ts.od<-apply(t(hm.df),2,max, na.rm=TRUE) #OD at  max  for each TS
    temp.lab<-as.numeric(as.factor(d1.df$temp))
    matplot(seq(from=0, to=(nrow(t(hm.df))-1)/2, by=0.5),t(hm.df),lwd=2,  type="l",lty=1,yaxt='n', col=ts.col2[temp.lab],xlab="Time (h)",bty="n", ylab="OD",ylim=c(0,0.5) )
    if(st.select=='5'){axis(side=2 ,at=c(0,0.1,0.2,0.3,0.4,0.5))    } #if yaxt='n'
    text(max.ts.time, max.ts.od+0.02, d1.df$temp, col=ts.col2[temp.lab], cex=0.75 )
    if(ana.select==0){title(paste0("TIGR4:",st.select))}
  }
}
dev.off()

#TIGR STRAINS, stratified by temp and oxygen
tiff('TIGR cps switch by temp and o2.tiff',width=11, height=5.3, units='in', res=400)
col3=c('#66c2a5','#fc8d62', '#8da0cb', '#e78ac3' )
par(mfcol=c(3,3),mar=c(2,3,1,1))
for(temp.select in c('33','35','37')){ #Body sites
  for(ana.select in c(0,2,1)){  #Catalase, anaerobic, aertobic
    id.composite<-paste0(d1a.agg$ID,d1a.agg$st,d1a.agg$temp,d1a.agg$anaerobic)
    d1.df<-d1a.agg[  d1a.agg$anaerobic %in% c(ana.select) & d1a.agg$Diagnosis %in% c(8)  & d1a.agg$st %in% c('5','14','19F','cps-') & d1a.agg$temp ==temp.select, ]
    d1.df<-d1.df[ order(-d1.df$anaerobic, -d1.df$temp, d1.df$st), ]
    hm.df<-(as.matrix(d1.df[ , 7:52]))
    max.ts.time<-apply(t(hm.df),2,which.max)/2 -0.5 #time at which max is achieved for each TS
    max.ts.od<-apply(t(hm.df),2,max, na.rm=TRUE) #OD at  max  for each TS
    st.lab<-as.numeric(as.factor(d1.df$st))
    matplot(seq(from=0, to=(nrow(t(hm.df))-1)/2, by=0.5),t(hm.df),lwd=2,  type="l",lty=1,yaxt='n', col=col3[st.lab],xlab="Time (h)",bty="n", ylab="OD",ylim=c(0,0.5) )
    if(temp.select=='33'){axis(side=2 ,at=c(0,0.1,0.2,0.3,0.4,0.5))    } #if yaxt='n'
    text(max.ts.time, max.ts.od+0.02, d1.df$st, col=col3[st.lab], cex=0.75 )
    if(ana.select==0){title(paste0(temp.select))}
  }
}
dev.off()

#15B CPS knockouts
tiff('cpsKO_15B.tiff',width=11, height=5.3, units='in', res=400)
col3=c('#66c2a5','#fc8d62', '#8da0cb', '#e78ac3' )
par(mfcol=c(3,3),mar=c(2,3,1,1))
id.select1<-c('15B-','CDC_15B')
id.select2<-c( '10A-', 'CDC10A')
for(temp.select in c('33','35','37')){ #Body sites
  for(ana.select in c(0,2,1)){  #Catalase, anaerobic, aertobic
    id.keep<- grepl(id.select1[1],d1a.agg$ID) + grepl(id.select1[2],d1a.agg$ID)
    id.composite<-paste0(d1a.agg$ID,d1a.agg$st,d1a.agg$temp,d1a.agg$anaerobic)
    d1.df<-d1a.agg[  id.keep==TRUE & d1a.agg$anaerobic %in% c(ana.select) & d1a.agg$temp ==temp.select, ]
    d1.df<-d1.df[ order(-d1.df$anaerobic, -d1.df$temp, d1.df$st), ]
    d1.df$st[d1.df$st=='15B-'] <-'cps-'
    d1.df$st[d1.df$st=='15BC'] <-'15B'
    d1.df$st[d1.df$st=='10A-'] <-'cps-'
    hm.df<-(as.matrix(d1.df[ , 7:52]))
    max.ts.time<-apply(t(hm.df),2,which.max)/2 -0.5 #time at which max is achieved for each TS
    max.ts.od<-apply(t(hm.df),2,max, na.rm=TRUE) #OD at  max  for each TS
    st.lab<-as.numeric(as.factor(d1.df$st))
    matplot(seq(from=0, to=(nrow(t(hm.df))-1)/2, by=0.5),t(hm.df),lwd=2,  type="l",lty=1,yaxt='n', col=col3[st.lab],xlab="Time (h)",bty="n", ylab="OD",ylim=c(0,0.5) )
    if(temp.select=='33'){axis(side=2 ,at=c(0,0.1,0.2,0.3,0.4,0.5))    } #if yaxt='n'
    text(max.ts.time, max.ts.od+0.02, d1.df$st, col=col3[st.lab], cex=0.75 )
    if(ana.select==0){title(paste0(temp.select))}
  }
}
dev.off()
#10A CPS knockouts
tiff('cpsKO 10A.tiff',width=11, height=5.3, units='in', res=400)
col3=c('#66c2a5','#fc8d62', '#8da0cb', '#e78ac3' )
par(mfcol=c(3,3),mar=c(2,3,1,1))
id.select1<-c('15B-','CDC_15B')
id.select2<-c( '10A-', 'CDC_10A')
for(temp.select in c('33','35','37')){ #Body sites
  for(ana.select in c(0,2,1)){  #Catalase, anaerobic, aertobic
    id.keep<- grepl(id.select2[1],d1a.agg$ID) + grepl(id.select2[2],d1a.agg$ID)
    id.composite<-paste0(d1a.agg$ID,d1a.agg$st,d1a.agg$temp,d1a.agg$anaerobic)
    d1.df<-d1a.agg[  id.keep==TRUE & d1a.agg$anaerobic %in% c(ana.select) & d1a.agg$temp ==temp.select, ]
    #d1.df<-d1_agg2[d1_agg2$st==st.select & d1_agg2$anaerobic %in% c(ana.select) & d1_agg2$Diagnosis %in% c(dx.select) &d1_agg2$ID %in% id.select.11a ,]
    d1.df<-d1.df[ order(-d1.df$anaerobic, -d1.df$temp, d1.df$st), ]
    d1.df$st[d1.df$st=='15B-'] <-'cps-'
    d1.df$st[d1.df$st=='15BC'] <-'15B'
    d1.df$st[d1.df$st=='10A-'] <-'cps-'
    hm.df<-(as.matrix(d1.df[ , 7:52]))
    max.ts.time<-apply(t(hm.df),2,which.max)/2 -0.5 #time at which max is achieved for each TS
    max.ts.od<-apply(t(hm.df),2,max, na.rm=TRUE) #OD at  max  for each TS
    st.lab<-as.numeric(as.factor(d1.df$st))
    matplot(seq(from=0, to=(nrow(t(hm.df))-1)/2, by=0.5),t(hm.df),lwd=2,  type="l",lty=1,yaxt='n', col=col3[st.lab],xlab="Time (h)",bty="n", ylab="OD",ylim=c(0,0.5) )
    if(temp.select=='33'){axis(side=2 ,at=c(0,0.1,0.2,0.3,0.4,0.5))    } #if yaxt='n'
    text(max.ts.time, max.ts.od+0.02, d1.df$st, col=col3[st.lab], cex=0.75 )
    if(ana.select==0){title(paste0(temp.select))}
  }
}
dev.off()

#603 6B CPS knockouts
tiff('cpsKO 603.tiff',width=11, height=5.3, units='in', res=400)
col3=c('#66c2a5','#fc8d62', '#8da0cb', '#e78ac3' )
par(mfcol=c(3,3),mar=c(2,3,1,1))
for(temp.select in c('33','35','37')){ #Body sites
  for(ana.select in c(0,2,1)){  #Catalase, anaerobic, aertobic
    id.keep<- grepl(id.select2[1],d1a.agg$ID) + grepl(id.select2[2],d1a.agg$ID)
    id.composite<-paste0(d1a.agg$ID,d1a.agg$st,d1a.agg$temp,d1a.agg$anaerobic)
    d1.df<-d1a.agg[  d1a.agg$Diagnosis==7 & d1a.agg$anaerobic %in% c(ana.select) & d1a.agg$temp ==temp.select, ]
    d1.df<-d1.df[ order(-d1.df$anaerobic, -d1.df$temp, d1.df$st), ]
    hm.df<-(as.matrix(d1.df[ , 7:52]))
    max.ts.time<-apply(t(hm.df),2,which.max)/2 -0.5 #time at which max is achieved for each TS
    max.ts.od<-apply(t(hm.df),2,max, na.rm=TRUE) #OD at  max  for each TS
    st.lab<-as.numeric(as.factor(d1.df$st))
    matplot(seq(from=0, to=(nrow(t(hm.df))-1)/2, by=0.5),t(hm.df),lwd=2,  type="l",lty=1,yaxt='n', col=col3[st.lab],xlab="Time (h)",bty="n", ylab="OD",ylim=c(0,0.5) )
    if(temp.select=='33'){axis(side=2 ,at=c(0,0.1,0.2,0.3,0.4,0.5))    } #if yaxt='n'
    text(max.ts.time, max.ts.od+0.02, d1.df$st, col=col3[st.lab], cex=0.75 )
    if(ana.select==0){title(paste0(temp.select))}
  }
}
dev.off()

##Calculate max OD600 and time to mid point OD600
d1a$anaerobic<- relevel(as.factor(d1a$anaerobic), ref = '1')
d1a$st<- relevel(as.factor(d1a$st), ref = '14') 
d1a$temp<- relevel(as.factor(d1a$temp), ref = '33') #reference
d1a$Diagnosis<- relevel(as.factor(d1a$Diagnosis), ref = '5') #Sepsis as reference
max.od<-apply(d1a[,-c(1:6)],1,max, na.rm=TRUE)
ts<-d1a[,-c(1:6)] #times series only
t.past.mid.log<- ts >= max.od/2 #is OD past mid point?
t.past.mid.log.cumsum<-apply(t.past.mid.log,1,cumsum)
t.past.mid.log.cumsum2<-apply(t.past.mid.log.cumsum,2,cumsum) #if multiple with value of 1, take first
t.midpoint<-unlist(apply(t.past.mid.log.cumsum2,2,function(x) which(x==1)/2))
d1a<-cbind.data.frame(max.od,t.midpoint,d1a)
d1a<-d1a[d1a$max.od>=0.05,]    #FILTER OUT IF MAX OD <0.05

#Format data for regression
d1b<-d1a[d1a$Diagnosis %in% c('0','5', '1','4'),]
d1b$Diagnosis<-factor(d1b$Diagnosis)
d1b$st<-factor(d1b$st)
d1b$temp<-factor(d1b$temp)
d1b$anaerobic<-factor(d1b$anaerobic)

#Make a plot showing the overall patterns of max OD by temp and site of isolation
library(reshape2)
sub1<-d1b[,c('temp','st','ID', 'Diagnosis','max.od','anaerobic')]
mo1.m<-melt(sub1, id=c('temp','st','ID', 'Diagnosis','anaerobic'))
mo1.c<-acast(mo1.m[mo1.m$variable=='max.od',], temp~ID~Diagnosis~anaerobic, fun.aggregate=mean, drop=FALSE)
mo1.c<-mo1.c[order(as.numeric(dimnames(mo1.c)[[1]])),,,]
trans.black=rgb(0,0,0,alpha=0.2)
par(mfrow=c(2,3))
matplot(as.numeric(dimnames(mo1.c)[[1]]), mo1.c[,,'5','0'], type='l',bty='l', ylim=c(0,0.6), main='IPD, aerobic no catalase', col=trans.black, lty=1) #IPD
matplot(as.numeric(dimnames(mo1.c)[[1]]), mo1.c[,,'5','2'], type='l',bty='l', ylim=c(0,0.6), main='IPD, aerobic with catalase', col=trans.black, lty=1)
matplot(as.numeric(dimnames(mo1.c)[[1]]), mo1.c[,,'5','1'], type='l',bty='l', ylim=c(0,0.6), main='IPD, Anaerobic', col=trans.black, lty=1)

matplot(as.numeric(dimnames(mo1.c)[[1]]), mo1.c[,,'0','0'], type='l',bty='l', ylim=c(0,0.6), main='Carriage, aerobic no catalase', col=trans.black, lty=1) #IPD
matplot(as.numeric(dimnames(mo1.c)[[1]]), mo1.c[,,'0','2'], type='l',bty='l', ylim=c(0,0.6), main='Carriage, aerobic with catalase', col=trans.black, lty=1)
matplot(as.numeric(dimnames(mo1.c)[[1]]), mo1.c[,,'0','1'], type='l',bty='l', ylim=c(0,0.6), main='Carriage, Anaerobic', col=trans.black, lty=1)

#d1b[d1b$ID=='H71',c(1:10)]
#d1b[d1b$ID=='2011M8',c(1:10)] #carriage strain from Dutch collection, only tested at 33C

###SPLINE ANALYSIS TO GET DERIVATIVES Use smooth.spline function, which is same thing they do in the grofit package; use a fixed smoothing parameter for simplicity
#library(fdapace)
Y.list<- split(d1a[,-c(1:8, (48+8+1):ncol(d1a))], as.factor(1:nrow(d1a)))
d1.data<-d1a[,-c(1:8, (48+8+1):ncol(d1a))]
Y.list<-lapply(Y.list, function(x) as.numeric(x[1,])) #List of OD600 values
T.list<-lapply(Y.list, function(x) (1:length(x))/2 ) #List of times (h)

spl1<-lapply(Y.list, function(x) smooth.spline(log(x[!is.na(x)]+0.01), spar=0.5)) #smooth the log-growth curve--gives mch more reliable results than smothing raw OD
  #Spar determined empircally by trying values between 0.2 and 0.9  to find one that was adquate for many curves (ie eyeballing where the estimated
  #start to log-grwth was vs the calculated time when 2nd derivative was at maximum. cross-validation
  #gave unstable 
fitted.all<-t(sapply(spl1, function(x) c(predict(x)$y, rep(NA, 48-length(predict(x)$y ))  )))
derivs.all<-t(sapply(spl1, function(x) c(predict(x, deriv=1)$y, rep(NA, 48-length(predict(x)$y ))  )))
derivs2.all<-t(sapply(spl1, function(x) c(predict(x, deriv=2)$y, rep(NA, 48-length(predict(x)$y ))  )))

#When is 2nd deriv at max? ie start of log phase
max.deriv2.time<-apply(derivs2.all,1,function(x) which(x==max(x[1:24], na.rm=TRUE))[1] ) #capture peak that occures within first 12 h
max.deriv2.time<-cbind.data.frame(d1a[,1:8],max.deriv2.time)

par(mfrow=c(4,1), mar=c(1,1,1,1))
for(i in 1:20){
  plot(log(Y.list[[i]]+0.01), type='l', bty='l', ylim=log(c(0.01, 0.5)))
  title(paste("Observed",i))
  abline(v=max.deriv2.time$max.deriv2.time[i], lty=2, col='gray')
  plot(fitted.all[i,], type='l', bty='l', ylim=log(c(0.01, 0.5)))
  abline(v=max.deriv2.time$max.deriv2.time[i], lty=2, col='gray')
  title("Smoothed")
  plot(derivs.all[i,], type='l', bty='l')
  abline(v=max.deriv2.time$max.deriv2.time[i], lty=2, col='gray')
  title("1st deriv")
  plot(derivs2.all[i,], type='l', bty='l')
  abline(v=max.deriv2.time$max.deriv2.time[i], lty=2, col='gray')
  title('2nd deriv')
  }
#d1a[6, c(1:8)]  #H14
par(mfrow=c(1,2), mar=c(1,2,1,1))
for( i in 1:2){
test1<-d1a[d1a$ID=='H14' & d1a$anaerobic %in% c(i),]
temp.lab<-as.numeric(as.factor(as.numeric(as.character(test1$temp))))
matplot(t(test1[, -c(1:8)]), type='l',col=ts.col[temp.lab], xlim=c(1,48),ylim=c(0,0.5), bty='l')
}
max.deriv<-apply(derivs.all,1,max, na.rm=T)

d1c<-cbind.data.frame(d1a[,c('max.od','ID','st','anaerobic','temp', 'Diagnosis')],max.deriv, max.deriv2.time)
d1c<


#Compare max.d and length of lag phase and max.deriv
comp1.gro<-cbind.data.frame(d1a[,c('max.od','ID','st','anaerobic','temp', 'Diagnosis')],max.deriv, max.deriv2.time)
cor(comp1.gro[,c('max.od','max.deriv','max.deriv2.time')])#growth rate correlated with max.od, not with lag length
#plot(comp1.gro$max.deriv, comp1.gro$max.od)
cor.test(comp1.gro[,c('max.od')],comp1.gro[,c('max.deriv')] )
cor.test(comp1.gro[,c('max.deriv2.time')],comp1.gro[,c('max.deriv')] )
cor.test(comp1.gro[,c('max.deriv2.time')],comp1.gro[,c('max.od')] )
##Does correlation differ with environmental conditions?
# comp1.gro.ls<-split(comp1.gro, list(comp1.gro$anaerobic, comp1.gro$temp))
# corr.grp<-lapply(comp1.gro.ls, function(x) cor(x[,c('max.od','max.deriv','max.deriv2.time')] ))


###############################################
#EVALuate 3 ways to measure serotype effect:
#1) Max OD achieved for a growth curve
#2) Time to max second derivative (when is growth rate increasing most)
#3) OD at average value for #2 for a given temp/anaerobic/Diagnostic level
#Test with aerobic and anaerobic together, and then separate, so have 9 
#measures total
##############################################
#ST Effect  on max.od  
mod.st<-lme(max.od  ~ st +  anaerobic*temp * Diagnosis ,  random = ~ 1|ID,data=d1b[d1b$anaerobic %in% c('1','2') & d1b$Diagnosis %in% c('0','5'),])
coef.fe.mod.st<-fixef(mod.st)
coef.st1.2<-coef.fe.mod.st[substr(names(coef.fe.mod.st),1,2)=='st']

mod.st<-lme(max.od  ~ st +  temp * Diagnosis ,  random = ~ 1|ID,data=d1b[d1b$anaerobic %in% c('1') & d1b$Diagnosis %in% c('0','5'),])
coef.fe.mod.st<-fixef(mod.st)
coef.st1<-coef.fe.mod.st[substr(names(coef.fe.mod.st),1,2)=='st']

mod.st<-lme(max.od  ~ st +  temp * Diagnosis ,  random = ~ 1|ID,data=d1b[d1b$anaerobic %in% c('2') & d1b$Diagnosis %in% c('0','5'),])
coef.fe.mod.st<-fixef(mod.st)
coef.st2<-coef.fe.mod.st[substr(names(coef.fe.mod.st),1,2)=='st']
coef.st.df3<-cbind.data.frame('coef.max.od.2'=coef.st2,'coef.max.od.1'=coef.st1,'coef.max.od.1.2'=coef.st1.2, 'st'=substring(names(coef.st1.2),3))
cor(coef.st.df3[,1:3])

###############################################
##############################################
##############################################
#ST Effect  on max.deriv2.time; ie time to max rate of growth
#This effect is stronger when look at anaerobic vs aerobic+catalase (-0.61 (P=0.02) vs -0.2 (P>0.05))
lag.length.mod<-lme(max.deriv2.time  ~ st +anaerobic+  temp  ,  random = ~ 1|ID,
                    data=max.deriv2.time[max.deriv2.time$anaerobic %in% c('1','2') 
                                         & max.deriv2.time$Diagnosis %in% c('5'),])
summary(lag.length.mod)
mod.st.t<-lme(max.deriv2.time  ~ st +  anaerobic*temp * Diagnosis ,  random = ~ 1|ID,
              data=max.deriv2.time[max.deriv2.time$anaerobic %in% c('1','2') 
                                   & max.deriv2.time$Diagnosis %in% c('0','5'),])
coef.fe.mod.st.time<-fixef(mod.st.t)
coef.st.time.1.2<-coef.fe.mod.st.time[substr(names(coef.fe.mod.st.time),1,2)=='st']

mod.st.t<-lme(max.deriv2.time  ~ st +  temp * Diagnosis ,  random = ~ 1|ID,
              data=max.deriv2.time[max.deriv2.time$anaerobic %in% c('1') 
                                   & max.deriv2.time$Diagnosis %in% c('0','5'),])
coef.fe.mod.st.time<-fixef(mod.st.t)
coef.st.time.1<-coef.fe.mod.st.time[substr(names(coef.fe.mod.st.time),1,2)=='st']

mod.st.t<-lme(max.deriv2.time  ~ st +  temp * Diagnosis ,  random = ~ 1|ID,
              data=max.deriv2.time[max.deriv2.time$anaerobic %in% c('2') 
                                   & max.deriv2.time$Diagnosis %in% c('0','5'),])
coef.fe.mod.st.time<-fixef(mod.st.t)
coef.st.time.2<-coef.fe.mod.st.time[substr(names(coef.fe.mod.st.time),1,2)=='st']
coef.st.time.df<-cbind.data.frame('coef.st.time1.2'=coef.st.time.1.2,'coef.st.time.2'=coef.st.time.2,'coef.st.time.1'=coef.st.time.1, 'st'=substring(names(coef.st.time.1),3))
cor(coef.st.time.df[,1:3])
###############################################
##############################################

#calculate OD at ave max.deriv2.time for a temp.aernarobic/diagnosis combo
ds.sub<-max.deriv2.time[max.deriv2.time$Diagnosis %in% c('0','5'),]
ds.sub$Diagnosis<-factor(ds.sub$Diagnosis)
spl.ds1<-split(ds.sub, list(ds.sub$anaerobic, ds.sub$temp, ds.sub$Diagnosis))
ave.max.deriv2.time<-round(sapply(spl.ds1, function(x) median(x$max.deriv2.time,na.rm=TRUE)))+2 # Look 2 time periods (1 hour) AFTEr the start of log-growth
labs1<-matrix(unlist(strsplit(names(ave.max.deriv2.time), '\\.')), ncol=3, byrow=T)
ave.max.deriv2.time<-cbind.data.frame(ave.max.deriv2.time, labs1)
names(ave.max.deriv2.time)<-c('ave.max.deriv2.time','anaerobic','temp','Diagnosis' )
ds1<-merge(d1b, ave.max.deriv2.time, by=c('anaerobic', 'temp','Diagnosis'))
od.fixed.time<-apply(ds1, 1, function(x) x[(8+as.numeric(x[length(x)] ))]) 
od.fixed.time<-cbind.data.frame(ds1[,1:8], od.fixed.time=as.numeric(as.character((od.fixed.time))))
hist(od.fixed.time$od.fixed.time)
##Models
mod.st.fixed.od<-lme(od.fixed.time  ~ st +  anaerobic*temp * Diagnosis ,  random = ~ 1|ID,
                     data=od.fixed.time[od.fixed.time$anaerobic %in% c('1','2') 
                                        & od.fixed.time$Diagnosis %in% c('0','5'),])
coef.fe.mod.st.fixed.od<-fixef(mod.st.fixed.od)
coef.st.fixed.od.1.2<-coef.fe.mod.st.fixed.od[substr(names(coef.fe.mod.st.fixed.od),1,2)=='st']

mod.st.fixed.od<-lme(od.fixed.time  ~ st +  temp*Diagnosis ,  random = ~ 1|ID,
                     data=od.fixed.time[od.fixed.time$anaerobic %in% c('1') 
                                        & od.fixed.time$Diagnosis %in% c('0','5'),])
coef.fe.mod.st.fixed.od<-fixef(mod.st.fixed.od)
coef.st.fixed.od.1<-coef.fe.mod.st.fixed.od[substr(names(coef.fe.mod.st.fixed.od),1,2)=='st']

mod.st.fixed.od<-lme(od.fixed.time  ~ st +  temp * Diagnosis ,  random = ~ 1|ID,
                     data=od.fixed.time[od.fixed.time$anaerobic %in% c('2') 
                                        & od.fixed.time$Diagnosis %in% c('0','5'),])
coef.fe.mod.st.fixed.od<-fixef(mod.st.fixed.od)
coef.st.fixed.od.2<-coef.fe.mod.st.fixed.od[substr(names(coef.fe.mod.st.fixed.od),1,2)=='st']
coef.st.fixed.od.df<-cbind.data.frame('coef.st.fixed.od.1.2'=coef.st.fixed.od.1.2,'coef.st.fixed.od.1'=coef.st.fixed.od.1, 
                                      'coef.st.fixed.od.2'=coef.st.fixed.od.2,'st'=substring(names(coef.st.fixed.od.2),3))
cor(coef.st.fixed.od.df[,1:3])

#################################################################################
#################################################################################
##REGRESSION MODELS
#Regression model; controls for repeated measurements of same isolate across conditions with random intercept
mod1<-lme(max.od  ~ st + anaerobic + temp + Diagnosis, random = ~ 1|ID, data=d1b)
summary(mod1)

#Test which serotypes grow best without catalase
mod2<-lme(max.od  ~ st +  temp + Diagnosis, random = ~ 1|ID, data=d1b[d1b$anaerobic==0 &d1b$Diagnosis %in% c(0,5),])
summary(mod2)

##Test whether effect of oxygen is different depending on sit eof origin of the strains
#Suggests carriage strains do better than IPD strains
tiff('fig 2 site o2 effect.tiff',width=7.5, height=2.7, units='in', res=200)
par(mfrow=c(1,3))
temp.ranges<-list( c(30,33,35,37,38,39),c(30,33,35),c(37,38,39))
for(i in 1:3){
  d1c<-d1b[d1b$temp %in% temp.ranges[[i]],]
  d1c$temp<-factor(d1c$temp)
  d1c$st<-factor(d1c$st)
  d1c$anaerobic<-factor(d1c$anaerobic)
  d1c$Diagnosis<-factor(d1c$Diagnosis)
  print(nrow(d1c))
  mod2<-lme(max.od  ~ st + anaerobic + temp + Diagnosis +Diagnosis*anaerobic,  random = ~ 1|ID,data=d1c)
  summary(mod2)
  test<- interactionMeans(mod2, factors=c('Diagnosis','anaerobic'))
  test<-test[test$Diagnosis %in% c('0','5', '1','4'),]
  test$Diagnosis<- relevel(as.factor(test$Diagnosis), ref = '5') #Sepsis as reference
  test$plot.index<-NA
  test$plot.index[test$Diagnosis=='5']<-1
  test$plot.index[test$Diagnosis=='0']<-2
  test$plot.index[test$Diagnosis=='1']<-3
  test$plot.index[test$Diagnosis=='4']<-4
  test$plot.index[test$anaerobic==1]<-test$plot.index[test$anaerobic==1]+0.1
  test$plot.index[test$anaerobic==0]<-test$plot.index[test$anaerobic==0]-0.1
  test$lcl<-test$`adjusted mean` -1.96*test$`std. error`
  test$ucl<-test$`adjusted mean` +1.96*test$`std. error`
  
  #Colors for plotting
  col.ana<-c('#1b9e77','#d95f02','#7570b3')
  pch.ana=c(15,16,17)
  
  plot(test$plot.index, test$`adjusted mean`, col=col.ana[test$anaerobic],bty='l',xaxt='n',xlab="", ylab="OD600", pch=pch.ana[test$anaerobic], ylim=c(0,0.4))
  arrows(test$plot.index,test$lcl,test$plot.index,test$ucl, code=3, angle=90, length=0.0,  col=col.ana[test$anaerobic])
  axis(side=1,at=c(1,2,3,4), labels=FALSE)
  text(c(1,2,3,4), par("usr")[3] - 0.05, labels=c('IPD','Carriage','Pneumonia','Conjunctivitis' ), srt=45, pos=1, xpd=TRUE)
  if(i==1){
    legend(0.8,0.43,cex=0.8 ,pch=pch.ana[test$anaerobic],bty='n',  col=col.ana[test$anaerobic],text.font=2,legend=c('Anaerobic', 'Aerobic without catalase','Aerobic with catalase'))
    title(expression(paste("30-39",degree,"C")) )
  }
  if(i==2){
    title(expression(paste("30-35",degree,"C")) )
  }
  if(i==3){
    title(expression(paste("37-39",degree,"C")) )
  }
}
dev.off()

##Test whether effect of temp is different depending on anaerobic, by body site
tiff('fig 3 temp o2 effect.tiff',width=7.5, height=2.7, units='in', res=200)
par(mfrow=c(1,3))
site.select<-list( unique(d1b$Diagnosis), 5,0 )
for(i in c(1:3)){
  sub<-d1b[d1b$Diagnosis %in% site.select[[i]] & d1b$st %in% st.keep.dx,] #0=carriage, 5=IPD
  sub$st<-factor(sub$st)
  sub$anaerobic<-factor(sub$anaerobic)
  sub$Diagnosis<-factor(sub$Diagnosis)
  print(nrow(sub))
  mod3<-lme(max.od  ~ st + anaerobic + temp + anaerobic*temp,  random = ~ 1|ID,data=sub)
  summary(mod3)
  test3<- interactionMeans(mod3, factors=c('anaerobic','temp'))
  #plot(test3)
  test3$temp.index<- as.numeric(relevel(as.factor(test3$temp), ref = '30')) #reference
  test3$temp.index[test3$anaerobic==1]<-test3$temp.index[test3$anaerobic==1]+0.1
  test3$temp.index[test3$anaerobic==0]<-test3$temp.index[test3$anaerobic==0]-0.1
  test3$lcl<-test3$`adjusted mean` -1.96*test3$`std. error`
  test3$ucl<-test3$`adjusted mean` +1.96*test3$`std. error`
  plot(test3$temp.index, test3$`adjusted mean`, col=col.ana[test3$anaerobic],bty='l', xaxt='n',pch=pch.ana[test3$anaerobic],xlab='', ylim=c(0,0.4), ylab='OD600')
  arrows(test3$temp.index,test3$lcl,test3$temp.index,test3$ucl, code=3, angle=90, length=0.0,  col=col.ana[test3$anaerobic])
  axis(side=1,at=c(1,2,3,4,5,6), labels=FALSE)
  temp.labels=c( expression(paste("30",degree,"C")), expression(paste("33",degree,"C")),expression(paste("35",degree,"C")),
                 expression(paste("37",degree,"C")), expression(paste("38",degree,"C")), expression(paste("39",degree,"C")))
  text(c(1,2,3,4,5,6), par("usr")[3] - 0.02, labels=temp.labels, srt=45, pos=1, xpd=TRUE)
  if(i==1){
    legend(0.8,0.43,cex=0.8 ,pch=pch.ana[test$anaerobic],bty='n',  col=col.ana[test$anaerobic],text.font=2,legend=c('Anaerobic', 'Aerobic without catalase','Aerobic with catalase'))
    title("All isolates" )
  }
  if(i==2){
    title("IPD" )
  }
  if(i==3){
    title("Carriage")
  }
}
dev.off()

#test whether effect of ST varies by anaerobic
sub<-d1b[d1b$Diagnosis %in% c(0,1,2,3,4,5,6) ,]
col.ana2<-c('#1b9e77','white','#7570b3')
pch.ana2=c(15,16,17)
sub$st<-factor(sub$st)
sub$anaerobic<-factor(sub$anaerobic)
sub$Diagnosis<-factor(sub$Diagnosis)
mod4<-lme(max.od  ~ st + anaerobic +Diagnosis+ temp +st*anaerobic,  random = ~ 1|ID,data=sub)
summary(mod4)
test4<- interactionMeans(mod4, factors=c('anaerobic','st'))

#For maile:
sub<-d1b[d1b$Diagnosis %in% c(0,1,2,3,4,5,6) &d1b$anaerobic %in% c(0,1,2)
         &d1b$temp %in% c(33,35,37),]
sub$anaerobic<-factor(sub$anaerobic)
sub$temp<-factor(sub$temp)
sub$Diagnosis<-factor(sub$Diagnosis)
mod4a.maile<-lme(max.od  ~ st + anaerobic +Diagnosis+ temp +st*anaerobic,  random = ~ 1|ID,data=sub)
summary(mod4a.maile)
test4.maile<- interactionMeans(mod4a.maile, factors=c('anaerobic','st'))
write.csv(test4.maile, 'st.effect.csv')

test4$st.num<-as.numeric(test4$st)
test4$lcl<-test4$`adjusted mean` -1.96*test4$`std. error`
test4$ucl<-test4$`adjusted mean` +1.96*test4$`std. error`
test4.plot<-test4
par(mfrow=c(1,1))
plot(test4.plot$st.num, test4.plot$`adjusted mean`, col=col.ana2[test4.plot$anaerobic],bty='l', pch=pch.ana2[test4.plot$anaerobic])
arrows(test4.plot$st.num,test4.plot$lcl,test4.plot$st.num,test4.plot$ucl, code=3, angle=90, length=0.0,  col=col.ana2[test4.plot$anaerobic])
text(test4.plot$st.num[test4.plot$anaerobic==2], test4.plot$`adjusted mean`[test4.plot$anaerobic==2]*1.2, test4.plot$st[test4.plot$anaerobic==2]) 
coefs<-t(coef(mod4)[1,])
keep.index.coef<-which(row.names(coefs)=='st1:anaerobic2'):length(coefs)
coef.lab<-dimnames(coefs)[[1]][keep.index.coef]
coefs.ana.interact<-cbind.data.frame(coef.lab,coefs[keep.index.coef])
coefs.ana.interact[,2]<-coefs.ana.interact[,2]-min(coefs.ana.interact[,2])
names(coefs.ana.interact)<-c('st','coef.ana.st')
coefs.ana.interact<-coefs.ana.interact[order(coefs.ana.interact$coef.ana.st),]
coefs.ana.interact$st<-gsub( ":.*$", "", coefs.ana.interact$st ) #exclude text after colon
coefs.ana.interact$st<-substr(coefs.ana.interact$st,3,nchar(coefs.ana.interact$st)) #remove "st" from the variable
plot(1:nrow(coefs.ana.interact), coefs.ana.interact$coef.ana.st, bty='l')
text(1:nrow(coefs.ana.interact), coefs.ana.interact$coef.ana.st+0.01,coefs.ana.interact$st )
#Plot to show effect of O2 by serotype
cis.mod4<-as.data.frame(intervals(mod4)[[1]])
keep.index.coef<-which(row.names(cis.mod4)=='st1:anaerobic2'):nrow(cis.mod4)
cis.mod4.int<-cis.mod4[keep.index.coef,]
cis.mod4.ref<-as.data.frame(t(c(0,0,0))) #add bak in ref serotype
names(cis.mod4.ref)<-names(cis.mod4.int)
cis.mod4.ref$st<-'14'
cis.mod4.int$st<-gsub( ":.*$", "", row.names(cis.mod4.int) ) #exclude text after colon
cis.mod4.int$st<-substr(cis.mod4.int$st,3,nchar(cis.mod4.int$st)) #remove "st" from the variable
cis.mod4.int<-rbind.data.frame(cis.mod4.ref,cis.mod4.int)
remove<-c("43/45/48","6B/D",'9',"25FA/38",'5' ) #our ST5 in CDC collection not actually ST5
cis.mod4.int<-cis.mod4.int[!cis.mod4.int$st %in% remove,]
cis.mod4.int<-cis.mod4.int[order(cis.mod4.int$est.),]
cis.mod4.int$st.num<-1:nrow(cis.mod4.int)
cis.mod4.int[,1:3]<-cis.mod4.int[,1:3] - min(cis.mod4.int[,2]) #subtract to center on 0
tiff('fig 4 st o2 effect.tiff',width=10.6, height=5.3, units='in', res=400)
plot(cis.mod4.int$st.num,cis.mod4.int[,2], bty='l', ylim=c(-0.02,0.18), pch=16, xaxt='n',xlab='',
     ylab="Max density aerobic (with catalase) - anaerobic")
arrows(cis.mod4.int$st.num,cis.mod4.int$lower,cis.mod4.int$st.num,cis.mod4.int$upper, code=3, angle=90, length=0.0,  col='gray')
text(cis.mod4.int$st.num,cis.mod4.int[,2]+0.01,cis.mod4.int$st ,srt=90)
dev.off()

#COmpare O2 effect to invasiveness
inv1<-read.csv('C:/Users/dmw63/Desktop/My documents h/INVASIVENESS/mcmc_invasive_single_stage.csv')
inv1$st=as.character(inv1$st)
inv1$st[inv1$st=='6A/C']<-'6A'
comp1<-merge(cis.mod4.int, inv1, by='st') 
par(mfrow=c(1,1))
symbols(comp1$est., comp1$log.inv.age1, sqrt(comp1$log.inv.prec.age1/pi),inches=0.35, bty='l',fg="white", bg=trans.red, xlab="O2 effect", ylab="log(Invasiveness")
text(comp1$est., comp1$log.inv.age1, comp1$st, cex=0.75)
reg.inv<-lm(comp1$log.inv.age1~ est., data=comp1, weights=log.inv.prec.age1)
summary(reg.inv)


#test whether effect of ST varies by anaerobic
tiff('fig 5 st o2 effect.tiff',width=15, height=5, units='in', res=200)
par(mfrow=c(2,1),mar=c(2,3,1,1))
st.position<-sort(as.character(unique(d1b$st)))
'%!in%' <- function(x,y)!('%in%'(x,y))
st.position<-st.position[st.position %!in% c('9', '6B/D', '43/45/48')]
serogrp<-as.numeric(gsub("([0-9]+).*$", "\\1", st.position))
serogrp[st.position=='cps-']<-length(serogrp)
sero.order<-cbind.data.frame(st.position,serogrp)
sero.order<-sero.order[order(sero.order$serogrp,sero.order$st.position),]
sero.order$sero.plot.N<-1:nrow(sero.order)
for(diag.select in c(0,5)){
  sub<-d1b[d1b$Diagnosis %in% c(diag.select) ,]
  col.ana2<-c('#1b9e77','#d95f02','#7570b3')
  pch.ana2=c(15,16,17)
  if (diag.select==5){sub<-sub[sub$st!='cps-',] }
  if (diag.select==0){ sub<-sub[sub$st!='14',]}
  sub$st<-factor(sub$st)
  sub$anaerobic<-factor(sub$anaerobic)
  sub$temp<-factor(sub$temp)
  sub$Diagnosis<-factor(sub$Diagnosis)
  table(sub$st, sub$anaerobic)
  print(nrow(sub))
  mod4<-lme(max.od  ~ st + anaerobic + temp +st*anaerobic,  random = ~ 1|ID,data=sub)
  #summary(mod4)
  test4<- interactionMeans(mod4, factors=c('anaerobic','st'))
  
  #plot(test4, multiple=FALSE, bty='l')
  test4$lcl<-test4$`adjusted mean` -1.96*test4$`std. error`
  test4$ucl<-test4$`adjusted mean` +1.96*test4$`std. error`
  #test4.plot<-test4[test4$anaerobic %in% c(1,2),]
  test4.plot<-test4
  test4.plot$st<-as.character(test4.plot$st)
  test4.plot<-merge(test4.plot,sero.order, by.x='st', by.y='st.position')
  test4.plot$sero.plot.N
  test4$sero.plot.N[test4$anaerobic==1]<-test4$sero.plot.N[test4$anaerobic==1]+0.1
  test4$sero.plot.N[test4$anaerobic==0]<-test4$sero.plot.N[test4$anaerobic==0]-0.1
  #cols=c('black','red')
  #pchs<-c(16,17)
  test4.plot$sero.plot.N[test4.plot$anaerobic==1]<-test4.plot$sero.plot.N[test4.plot$anaerobic==1]+0.1
  test4.plot$sero.plot.N[test4.plot$anaerobic==0]<-test4.plot$sero.plot.N[test4.plot$anaerobic==0]-0.1
  plot(test4.plot$sero.plot.N, test4.plot$`adjusted mean`, col=col.ana2[test4.plot$anaerobic],bty='l',xaxt='n', xlab='', ylab='OD600',pch=pch.ana2[test4.plot$anaerobic], xlim=c(1,nrow(sero.order)),ylim=c(0,0.5))
  arrows(test4.plot$sero.plot.N,test4.plot$lcl,test4.plot$sero.plot.N,test4.plot$ucl, code=3, angle=90, length=0.0,  col=col.ana2[test4.plot$anaerobic])
  text(test4.plot$sero.plot.N[test4.plot$anaerobic==2], test4.plot$`adjusted mean`[test4.plot$anaerobic==2]*1.2, test4.plot$st[test4.plot$anaerobic==2], xpd=TRUE) 
  if(diag.select==0){title(main='Carriage')}
  if(diag.select==5){title(main='IPD')}
  
}
dev.off()

#effect of serotype by temp
mod5<-lm(max.od  ~ st + anaerobic + temp + Diagnosis +Diagnosis*anaerobic +st*temp,data=d1b)
summary(mod5)
mod5.r2<-summary(mod5)$adj.r.squared
coef.temp<-coef(mod5)
keep.coef.temp<-grep(":temp38", names(coef.temp),fixed=TRUE)
coef.temp.38<-coef.temp[keep.coef.temp] -min(coef.temp[keep.coef.temp], na.rm=TRUE )
confint(mod5)

#combine estimates for 33 vs 39; anaerobic vs catalase in one dataframe
coefs.interact.st<-cbind.data.frame(coefs.ana.interact,coef.temp.38)

#EFFECT OF AMBIENT AIR BY BODY SITE OF ISOLATE
mod6<-lme(max.od  ~ st + anaerobic + Diagnosis+ Diagnosis*anaerobic,  random = ~ 1|ID,data=d1b[d1b$temp==33,])

##restrict to those serotypes present in both carriage and IPD
mod6a<-lme(max.od  ~ st + anaerobic + Diagnosis+ Diagnosis*anaerobic,  random = ~ 1|ID,data=d1b[d1b$temp==38 & d1b$st  %in% st.keep.dx,])
summary(mod6a)
###############################################
##############################################  
###############################################
##############################################  
###############################################
##############################################


##what isolates are tested under what conditions?
# check1<-d1b[,c('ID', 'temp', 'anaerobic')]
# check1$one<-1
# check1.m<-melt(check1, id=c('ID', 'temp', 'anaerobic'))
# check1.c<-acast(check1.m, ID~anaerobic+temp)
# 


##Coefficients for Maile
#MERGE together all coefficients
maile.coef1<-merge(coef.st.time.df,coef.st.fixed.od.df, by='st', all=T )
maile.coef2<-merge(coef.st.df3,maile.coef1, by='st', all=T )
maile.coef3<-maile.coef2[,c('st', "coef.max.od.1","coef.max.od.2", "coef.st.time1.2", "coef.st.fixed.od.1.2")]
write.csv(maile.coef3, 'st.reg.coeffs.csv')

######CORRELATE ST COEFFFICIENTS VS OTHER VALUES
cov1<-read.csv('C:/Users/dmw63/Desktop/My documents h/OLD LIPSITCH LAB DATA/pneumo_master v2.csv')
cov1$st<-as.character(cov1$ST)
cov1$MASSCARR01[is.na(cov1$MASSCARR01)]<-0
cov1$OSMANCARRN[is.na(cov1$OSMANCARRN)]<-0
cov1$SLEEMANCARRN[is.na(cov1$SLEEMANCARRN)]<-0
cov1$NORWAYCARR_PRE[is.na(cov1$NORWAYCARR_PRE)]<-0
cov1$TROTTERCARRN[is.na(cov1$TROTTERCARRN)]<-0
cov1$GREECECARRN[is.na(cov1$GREECECARRN)]<-0
cov1$ave_carr_pre<-apply( cbind(cov1$MASSCARR01, cov1$OSMANCARRN, cov1$SLEEMANCARRN, cov1$NORWAYCARR_PRE, cov1$TROTTERCARRN,
                                cov1$GREECECARRN), 1, mean, na.rm=TRUE )
cov1$st[cov1$st=='6A/C']<-'6A'


#MERGE together all coefficients
carr.comp1<-merge(coef.st.time.df,cov1, by='st', all.y=T )
carr.comp2<-merge(coef.st.fixed.od.df,carr.comp1, by='st', all.y=T )
carr.comp<-merge(coef.st.df3,carr.comp2, by='st', all.y=T )

#Corr with carriage data
carr.comp$coef.st[carr.comp$st=='14']<-0 #ref
carr.comp$coef.st.time[carr.comp$st=='14']<-0 #ref
reg.pc1<-glm(sqrt(carr.comp$ave_carr_pre)~ coef.st.fixed.od.2, data=carr.comp)
summary(reg.pc1)
reg.pc2<-glm(sqrt(carr.comp$ave_carr_pre)~ coef.st.fixed.od.1, data=carr.comp)
summary(reg.pc2)
plot(carr.comp$coef.st.fixed.od.1, carr.comp$ave_carr_pre)
reg.pc3<-glm(sqrt(carr.comp$ave_carr_pre)~ coef.max.od.2, data=carr.comp)
summary(reg.pc3)
reg.pc4<-glm(sqrt(carr.comp$ave_carr_pre)~ coef.max.od.1, data=carr.comp)
summary(reg.pc4)

reg.pc5<-glm(sqrt(carr.comp$ave_carr_pre)~ coef.st.time.2 , data=carr.comp)
summary(reg.pc5)
reg.pc6<-glm(sqrt(carr.comp$ave_carr_pre)~ coef.st.time.1 , data=carr.comp)
summary(reg.pc6)
par(mfrow=c(1,1))
plot(carr.comp$coef.st.fixed.od.1, carr.comp$ave_carr_pre, col='white', bty='l')
text(carr.comp$coef.st.fixed.od.1, carr.comp$ave_carr_pre, carr.comp$st)   


#CFR Harboe
reg.pc.cfr1<-glm(carr.comp$CFRHARBOE~ coef.max.od.1 , data=carr.comp)
summary(reg.pc.cfr1)
reg.pc.cfr2<-glm(carr.comp$CFRHARBOE~ coef.max.od.2 , data=carr.comp)
summary(reg.pc.cfr2)
reg.pc.cfr3<-glm(carr.comp$CFRHARBOE~ coef.st.fixed.od.1 , data=carr.comp)
summary(reg.pc.cfr3)
reg.pc.cfr4<-glm(carr.comp$CFRHARBOE~ coef.st.fixed.od.2 , data=carr.comp)
summary(reg.pc.cfr4)
reg.pc.cfr5<-glm(carr.comp$CFRHARBOE~ coef.st.time.1 , data=carr.comp)
summary(reg.pc.cfr5)
reg.pc.cfr6<-glm(carr.comp$CFRHARBOE~ coef.st.time.2 , data=carr.comp)
summary(reg.pc.cfr6)

plot(carr.comp$coef.st.fixed.od.1, carr.comp$CFRHARBOE, col='white')
text(carr.comp$coef.st.fixed.od.1, carr.comp$CFRHARBOE, carr.comp$st)

#Total carbon
reg.pc.carbon1<-glm(carr.comp$TOTALCARBONREPEAT~ coef.max.od.1, data=carr.comp)
summary(reg.pc.carbon1)
reg.pc.carbon2<-glm(carr.comp$TOTALCARBONREPEAT~ coef.max.od.2, data=carr.comp)
summary(reg.pc.carbon2)
plot(carr.comp$coef.max.od.2, carr.comp$TOTALCARBONREPEAT, col='white')
text(carr.comp$coef.max.od.2, carr.comp$TOTALCARBONREPEAT, carr.comp$st)

#Invasiveness AJE paper--inverse variance weights
inv1<-read.csv('C:/Users/dmw63/Desktop/My documents h/INVASIVENESS/mcmc_invasive_single_stage.csv')
inv1$st=as.character(inv1$st)
inv1$st[inv1$st=='6A/C']<-'6A'
inv2<-merge(carr.comp, inv1, by='st' , all=TRUE)
trans.red<-rgb(1,0,0, alpha=0.5)
symbols(inv2$coef.st.fixed.od.1, inv2$log.inv.age1, sqrt(inv2$log.inv.prec.age1/pi),inches=0.35, bty='l',fg="white", bg=trans.red, xlab="Early od", ylab="log(Invasiveness")
text(inv2$coef.st.fixed.od.1, inv2$log.inv.age1, inv2$st, cex=0.75)
reg.inv<-lm(inv2$log.inv.age1~ coef.max.od.2, data=inv2, weights=log.inv.prec.age1)
summary(reg.inv)
reg.inv2<-lm(inv2$log.inv.age1~ coef.max.od.1, data=inv2, weights=log.inv.prec.age1)
summary(reg.inv2)
reg.inv3<-lm(inv2$log.inv.age1~ coef.st.fixed.od.1, data=inv2, weights=log.inv.prec.age1)
summary(reg.inv3)
reg.inv4<-lm(inv2$log.inv.age1~ coef.st.fixed.od.2, data=inv2, weights=log.inv.prec.age1)
summary(reg.inv4)
reg.inv5<-lm(inv2$log.inv.age1~ coef.st.time.2, data=inv2, weights=log.inv.prec.age1)
summary(reg.inv5)

#PS components
ps1<-read.csv('C:/Users/dmw63/Desktop/My documents h/LAB/pneumo metabolic genes/PS Composition_SS_final.csv')
ps1$st<-as.character(ps1$Serotype)
ps2<-merge(carr.comp, ps1, by='st' , all=TRUE)
ps2$NAC<-0
ps2$NAC[ps2$GlcNAc==1 | ps2$GalNAc==1 | ps2$ManNAc==1 |  ps2$ManNAcA==1 |ps2$FucNAc==1 |ps2$PneNAc ==1] <-1
ps2$uronic<-0
ps2$uronic[ps2$GalA==1 |ps2$GlcA==1 ]<-1  #GlcA and GalA derived from same pathway
reg.nac1<-glm(ps2$NAC~ coef.max.od.1, data=ps2, family='binomial')
summary(reg.nac1)
reg.nac2<-glm(ps2$NAC~ coef.max.od.2, data=ps2, family='binomial')
summary(reg.nac2)
reg.nac3<-glm(ps2$NAC~ coef.st.fixed.od.1, data=ps2, family='binomial')
summary(reg.nac3)
reg.nac4<-glm(ps2$NAC~ coef.st.fixed.od.2, data=ps2, family='binomial')
summary(reg.nac4)
reg.nac5<-glm(ps2$NAC~ coef.st.time.2, data=ps2, family='binomial')
summary(reg.nac5)

plot(ps2$NAC, ps2$coef.st)
wilcox.test(coef.st ~ NAC, data=ps2[!is.na(ps2$coef.st), ]) 
wilcox.test(coef.st ~ GlcNAc, data=ps2[!is.na(ps2$coef.st), ]) 

reg.uronic1<-glm(ps2$uronic~ coef.max.od.1, data=ps2, family='binomial')
summary(reg.uronic1)
reg.uronic2<-glm(ps2$uronic~ coef.max.od.2, data=ps2, family='binomial')
summary(reg.uronic2)
reg.uronic3<-glm(ps2$uronic~ coef.st.fixed.od.1, data=ps2, family='binomial')
summary(reg.uronic3)    
reg.uronic4<-glm(ps2$uronic~ coef.st.time.1 , data=ps2, family='binomial')
summary(reg.uronic4)      

wilcox.test(coef.st ~ uronic, data=ps2[!is.na(ps2$coef.st), ]) 
table(ps2$uronic, ps2$st)
plot(ps2$uronic, ps2$coef.st)

#multivariate with both uronic and coef.st as predictors
reg.both<-lm( coef.max.od.1~ NAC+uronic, data=ps2)
summary(reg.both)

table(ps2$uronic, ps2$NAC)


#Figure 6 OF PC2 vs components
tiff('fig 6.tiff',width=8, height=4, units='in', res=300)

par(mfrow=c(1,2), mar=c(4,4,1,1))
plot(carr.comp$coef.st, carr.comp$ave_carr_pre, col='white', bty='l', xlab='Serotype effect on growth (Density at early time point)' ,ylab='Pre-vaccine carriage')
text(carr.comp$coef.st, carr.comp$ave_carr_pre, carr.comp$st, xpd=NA)

symbols(inv2$coef.st, inv2$log.inv.age1, sqrt(inv2$log.inv.prec.age1/pi),inches=0.1, bty='l',fg="white", bg=trans.red, xlab="PC2 (Density at early time point)", ylab="log(Invasiveness)", bty='l')
text(inv2$coef.st, inv2$log.inv.age1, inv2$st, cex=0.75, col='gray', adj=c(1,1), xpd=NA)

dev.off()


############################################################
############################################################
#OTHER ANALYSES
#############################################################
##SIMPLE version--calculaye max OD by serotype, body site, ambient vs cat
mean1<-aggregate(d1b$max.od, list(d1b$st, d1b$Diagnosis, d1b$anaerobic, d1b$temp), mean)
names(mean1)<-c('st','Diagnosis','anaerobic','temp','max.od')
m.mean1<-melt(mean1, id=c("st",'Diagnosis','anaerobic','temp'))
c.mean1<-as.data.frame(cast(m.mean1, st+temp+Diagnosis~anaerobic))
c.mean1<-c.mean1[c.mean1$`2`>=0.05,]
c.mean1$ratio_1_2<-c.mean1$`1`/c.mean1$`2`  #Anerobic vs aerobi with catalase
c.mean1$ratio_1_0<-c.mean1$`1`/c.mean1$`0`  #Anerobic vs aerobi without catalase
plot(c.mean1$ratio_1_2,c.mean1$ratio_1_0)
c.mean1<-c.mean1[order(-c.mean1$ratio_1_2),]
c.mean1[c.mean1$temp==33 &c.mean1$Diagnosis==0,] #carriage isolates
c.mean1[c.mean1$temp==33 &c.mean1$Diagnosis==5,] #IPD isolates

c.mean1$log.ratio<-log(c.mean1$ratio_1_2)
t1<-lm(log.ratio~st+temp+Diagnosis, data=c.mean1[c.mean1$log.ratio<=10 & c.mean1$log.ratio>=-10,])
summary(t1)

#Test if ratio of ST at high temps is associated with invasiveness
cov1<-read.csv('C:/Users/dmw63/Desktop/My documents h/OLD LIPSITCH LAB DATA/pneumo_master v2.csv')
cov1$st<-as.character(cov1$ST)
c.mean1$st<-as.character(c.mean1$st)
c.mean2<-merge(c.mean1,cov1, by='st', all=TRUE)
c.mean.sub<-c.mean2[c.mean2$Diagnosis %in% c(0)& c.mean2$temp %in% c(33), ]
plot(c.mean.sub$CFRHARBOE, c.mean.sub$log.ratio, col='white')
text(c.mean.sub$CFRHARBOE, c.mean.sub$log.ratio, c.mean.sub$st)
cor(c.mean.sub$CFRHARBOE, c.mean.sub$log.ratio, method='spearman', use='pairwise.complete.obs')
plot(sqrt(c.mean.sub$ATTACKRATESLEEMAN), c.mean.sub$log.ratio, col='white')
text(sqrt(c.mean.sub$ATTACKRATESLEEMAN), c.mean.sub$log.ratio, c.mean.sub$st)
cor( c.mean.sub[,c('log.ratio','CFRHARBOE','TOTALCARBONREPEAT', 'GREECECARRN','BRUEGCARRN','ATTACKRATESLEEMAN','DUTCHCARRN','OSMANCARRN','SLEEMANCARRN')], method='spearman', use='pairwise.complete.obs')
plot(sqrt(c.mean.sub$TOTALCARBONREPEAT), c.mean.sub$log.ratio, col='white')
text(sqrt(c.mean.sub$TOTALCARBONREPEAT), c.mean.sub$log.ratio, c.mean.sub$st)

#Compare to regression coefficients, which averages over and adjust fors all disease states and temperatures
coef1<-summary(t1)$coefficients[,'Estimate']
coef.st<-coef1[substr(names(coef1),1,2)=='st']
coef.st.label<-substring(names(coef.st), 3)
coef.st<-cbind.data.frame(coef.st,'st'=coef.st.label)
coef.cor<-merge(coef.st,cov1, by='st', all=TRUE)
cor(coef.cor[,c('coef.st','CFRHARBOE','TOTALCARBONREPEAT', 'GREECECARRN','BRUEGCARRN','ATTACKRATESLEEMAN','DUTCHCARRN','OSMANCARRN','SLEEMANCARRN')], method='spearman', use='pairwise.complete.obs')
plot(sqrt(coef.cor$CFRHARBOE), coef.cor$coef.st, col='white')
text(sqrt(coef.cor$CFRHARBOE), coef.cor$coef.st, c.mean.sub$st)

ps1<-read.csv('C:/Users/dmw63/Desktop/My documents h/LAB/pneumo metabolic genes/PS Composition_SS_final.csv')
ps1$st<-as.character(ps1$Serotype)
ps2<-merge(coef.cor, ps1, by='st' , all=TRUE)
ps2$NAC<-0
ps2$NAC[ps2$GlcNAc==1 | ps2$GalNAc==1 | ps2$ManNAc==1 |  ps2$ManNAcA==1 |ps2$FucNAc==1 ] <-1
wilcox.test(coef.st ~ NAC, data=ps2) 
wilcox.test(coef.st ~ GlcNAc, data=ps2) 
wilcox.test(coef.st ~ GlcA, data=ps2) 

################################################
#######corr with carriage
d1.sub<-d1a[,c(1,2,3,4,5,7:48)]
d1.sub<-d1.sub[d1.sub$Diagnosis %in% c(0,5),]
d1.sub$conditions<-paste(d1.sub$anaerobic ,d1.sub$temp, d1.sub$Diagnosis, sep="_" )
d1.sub<-d1.sub[,-which(names(d1.sub) %in% c('anaerobic','temp','t.midpoint','Diagnosis','ID2','ID'))]
table(d1.sub$conditions, d1.sub$st)
test<-melt(d1.sub, id=c('st','conditions'), na.rm=TRUE)
test2.1<-cast(test, st~variable+conditions, mean)
test2a<-as.data.frame(test2.1)
miss.prop<-apply(test2a,2, function(x) sum(is.nan(x)) )/nrow(test2a)
test2<-test2a[,miss.prop<0.5] #remove columns for which there are bservations for < half serotypes
mat.text2<-t(test2[,-1])
#SIMPLE IMPUTE MISSING VALUES FROM MATRIX...by row (first transpose, so we are taking mean of time point and filling that in)
impute.matrix <- function (matrix) {  ##https://gist.github.com/Jfortin1/d4888d68359a36fbda60
  missing <- which(is.na(matrix) | !is.finite(matrix), arr.ind=TRUE)
  if (length(missing)!=0){
    for (j in 1:nrow(missing)){
      mean <- mean(matrix[missing[j,1],][is.finite(matrix[missing[j,1],])], na.rm=TRUE) #Mean of row
      matrix[missing[j,1],missing[j,2]] <- mean
    }
  }
  matrix
}
test3<-as.matrix(t(impute.matrix(mat.text2)))
test3.df<-as.data.frame(test3)
names(test3.df)<-names(test2)[-1]
test4<-cbind.data.frame(factor(test2$st), test3.df)
names(test4)[1]<-'st'
#Compile carriage data into average
cov1$MASSCARR01[is.na(cov1$MASSCARR01)]<-0
cov1$OSMANCARRN[is.na(cov1$OSMANCARRN)]<-0
cov1$SLEEMANCARRN[is.na(cov1$SLEEMANCARRN)]<-0
cov1$NORWAYCARR_PRE[is.na(cov1$NORWAYCARR_PRE)]<-0
cov1$TROTTERCARRN[is.na(cov1$TROTTERCARRN)]<-0
cov1$GREECECARRN[is.na(cov1$GREECECARRN)]<-0
cov1$ave_carr_pre<-apply( cbind(cov1$MASSCARR01, cov1$OSMANCARRN, cov1$SLEEMANCARRN, cov1$NORWAYCARR_PRE, cov1$TROTTERCARRN,
                                cov1$GREECECARRN), 1, mean, na.rm=TRUE )
cov1$st[cov1$st=='6A/C']<-'6A'
carr.comp<-merge(test4,cov1, by='st' )
#Test all time points, temps, and O2 levels, extract aic scores