### PREPARATION ###

#Set working directory, load packages
setwd("~/Mazama gouazoubira/Data files")
mazdata<-read.table("TimesData_MazamaJanv3.txt",header=T)
bosdata<-read.table("Bos.txt",header=T)
library(activity)
library(overlap)
library(beepr)
library(circular)

#Prep
radian=2*pi*mazdata$General 
radianbos=2*pi*bosdata$General
mazama=mazdata$Species=="Mazama"& mazdata$Excluded=="NO"
bos=bosdata$Species=="Bos"&bosdata$Excluded=="NO"
bos2=bosdata$Species=="Bos"&bosdata$Excluded=="No"&bosdata$Mazama=="Yes"

summary(bos)
summary(bos2)

bosdf=bosdata$Vegetation=="DryForest"& bosdata$Excluded=="NO"&bosdata$Mazama=="Yes"
bossv=bosdata$Vegetation=="Savanna"& bosdata$Excluded=="NO"&bosdata$Mazama=="Yes"
bososv=bosdata$Vegetation=="OpenSavanna"& bosdata$Excluded=="NO"&bosdata$Mazama=="Yes"
dryforest=mazdata$Vegetation=="DryForest"& mazdata$Excluded=="NO"

savanna=mazdata$Vegetation=="Savanna"& mazdata$Excluded=="NO"
opensavanna=mazdata$Vegetation=="OpenSavanna"& mazdata$Excluded=="NO"
cattlepresent=mazdata$Cattle=="Yes"& mazdata$Excluded=="NO"
cattleabsent=mazdata$Cattle=="No"& mazdata$Excluded=="NO"


### COMPLETE MAZAMA ACTIVITY PATTERN ###

mazrad<-radian[mazama] #Convert time to radians
length(mazrad) #Find number of records
mazmod<-fitact(mazrad,reps=10000,sample="model",show=TRUE);beep() #Fit distribution to time-of-detection data
show(mazmod) #Show activity level estimate 
plot(mazmod) #Plot activity pattern

jpeg(file="mazmod.jpeg")
plot(mazmod)
abline(v=c(6,17+59/60),lty=3)
dev.off()

### SUBSETTING MAZAMA ACTIVITY PATTERN BY VEGETATION TYPE ###

# Mazama activity: dry forest
mazrad_df<-radian[mazama&dryforest] 
length(mazrad_df)
mazmod_df<-fitact(mazrad_df,reps=10000,sample="model",show=TRUE);beep()
show(mazmod_df)

jpeg(file="mazmod_df.jpeg")
plot(mazmod_df)
abline(v=c(6,17+59/60),lty=3)
dev.off()


# Mazama activity: savanna 
mazrad_sv<-radian[mazama&savanna]
length(mazrad_sv)
mazmod_sv<-fitact(mazrad_sv,reps=10000,sample="model",show=TRUE);beep()
show(mazmod_sv)
jpeg(file="mazmod_sv.jpeg")
plot(mazmod_sv)
abline(v=c(6,17+59/60),lty=3)
dev.off()


# Mazama activity : open savanna 
mazrad_osv<-radian[mazama&opensavanna]
length(mazrad_osv)
mazmod_osv<-fitact(mazrad_osv,reps=10000,sample="model",show=TRUE);beep()
show(mazmod_osv)
jpeg(file="mazmod_osv.jpeg")
plot(mazmod_osv)
abline(v=c(6,17+59/60),lty=3)
dev.off()

#Overlap coefficients preparation
bsmazrad_sv<-resample(mazrad_sv,10000)
bsmazrad_df<-resample(mazrad_df,10000)
bsmazrad_osv<-resample(mazrad_osv,10000)
dim(bsmazrad_sv)
dim(bsmazrad_df)
dim(bsmazrad_osv)

#COMPARE ACTIVITY LEVELS AND OVERLAP COEFFICIENTS

#Compare savanna and dry forest
compareAct(list(mazmod_sv,mazmod_df))
mazsv_v_mazdf<-overlapPlot(mazrad_sv,mazrad_df, extend='lightsteelblue4',rug=TRUE)
overlapPlot(mazrad_sv,mazrad_df, extend='lightsteelblue4',rug=TRUE,lty=c(1,5),linecol=c("black","darkgreen"))
legend('topright', c("Savanna","Dry Forest"),lty=c(1,3,5),y.intersp=1,col=c("black","darkgreen"))
legend('bottom', c("???=0.811(0.736-0.890),U???=0.4261,p<0.001"), y.intersp = 2, x.intersp = 0.01, bty="n")
abline(v=c(9,16+59/60),lty=3)
dev.off()
(dhats_svdf<-overlapEst(mazrad_sv,mazrad_df))
bsout_svdf<-bootEst(bsmazrad_sv,bsmazrad_df,adjust=c(NA,1,NA))
dim(bsout_svdf)
colMeans(bsout_svdf)
bs_svdf<-as.vector(bsout_svdf[,2])
bootCI(dhats_svdf[2],bs_svdf)['basic0',]

watson.two.test(mazrad_df,mazrad_osv)
watson.two.test(mazrad_sv,mazrad_df,alpha=0.05)


watson.two.test(mazrad_osv,mazrad_df)
watson.two.test(mazrad_osv,mazrad_sv)
watson.two.test(mazrad_dfca,mazrad_osvca,alpha=0.05)
#Compare savanna and open savanna
compareAct(list(mazmod_sv,mazmod_osv))
overlapPlot(mazrad_sv,mazrad_osv,extend='lightsteelblue4',linetype=c(1,3),linecol=c("black","blue"),rug=TRUE)
legend('topright', c("Savanna","Open Savanna"),lty=c(1,3),y.intersp=0.8,col=c("black","blue"))
abline(v=c(6,17+59/60),lty=3)
(dhats_svosv<-overlapEst(mazrad_sv,mazrad_osv))
bsout_svosv<-bootEst(bsmazrad_sv,bsmazrad_osv,adjust=c(NA,1,NA))
dim(bsout_svosv)
colMeans(bsout_svosv)
bs_svosv<-as.vector(bsout_svosv[,2])
bootCI(dhats_svosv[2],bs_svosv)['basic0',]

#Compare open savanna and dry forest
compareAct(list(mazmod_osv,mazmod_df))
overlapPlot(mazrad_osv,mazrad_df, extend='lightsteelblue4',rug=TRUE, linetype=c(3,5), linecol=c("blue","darkgreen"))
legend('topright', c("Open Savanna","Dry Forest"),lty=c(3,5),y.intersp=0.8,col=c("blue","darkgreen"))
abline(v=c(6,17+59/60),lty=3)
(dhats_osvdf<-overlapEst(mazrad_osv,mazrad_df))
bsout_osvdf<-bootEst(bsmazrad_osv,bsmazrad_df,adjust=c(NA,1,NA))
dim(bsout_osvdf)
colMeans(bsout_osvdf)
bs_osvdf<-as.vector(bsout_osvdf[,2])
bootCI(dhats_osvdf[2],bs_osvdf)['basic0',]

### SUBSETTING MAZAMA ACTIVITY BY CATTLE ###

# Mazama activity: cattle present
mazrad_cp<-radian[mazama&cattlepresent]
length(mazrad_cp)
mazmod_cp<-fitact(mazrad_cp,reps=10000,sample="model",show=TRUE);beep()
show(mazmod_cp)
plot(mazmod_cp)
abline(v=c(6,17+59/60),lty=3)

# Mazama activity: cattle absent
mazrad_ca<-radian[mazama&cattleabsent]
length(mazrad_ca)
mazmod_ca<-fitact(mazrad_ca,reps=10000,sample="model",show=TRUE);beep()
plot(mazmod_ca)
show(mazmod_ca)
abline(v=c(6,17+59/60),lty=3)

#COMPARE ACTIVITY LEVELS AND OVERLAP COEFFICIENTS

compareAct(list(mazmod_cp,mazmod_ca))
overlapPlot(mazrad_cp,mazrad_ca,extend='lightsteelblue4')
legend('topright',c("Cattle present","Cattle absent"),lty=c(1,2),col=c("black","blue"),bty='n')
abline(v=c(6,17+59/60),lty=3)
(dhats_cpca<-overlapEst(mazrad_cp,mazrad_ca))
bsmazrad_cp<-resample(mazrad_cp,100)
bsmazrad_ca<-resample(mazrad_ca,100)
dim(bsmazrad_cp)
dim(bsmazrad_ca)
bsout_cpca<-bootEst(bsmazrad_cp,bsmazrad_ca,adjust=c(NA,1,NA))
dim(bsout_cpca)
colMeans(bsout_cpca)
bs_cpca<-as.vector(bsout_cpca[,2])
bootCI(dhats_cpca[2],bs_cpca)['basic0',]

### SUBSETTING MAZAMA ACTIVITY BY VEGETATION & CATTLE ###

#Convert times to radian
mazrad_dfcp<-radian[mazama&dryforest&cattlepresent]
mazrad_dfca<-radian[mazama&dryforest&cattleabsent]
mazrad_svcp<-radian[mazama&savanna&cattlepresent]
mazrad_svca<-radian[mazama&savanna&cattleabsent]
mazrad_osvcp<-radian[mazama&opensavanna&cattlepresent]
mazrad_osvca<-radian[mazama&opensavanna&cattleabsent]

# Mazama activity: dry forest, cattle present
length(mazrad_dfcp)

mazmod_dfcp<-fitact(mazrad_dfcp,reps=10000,sample="data",show=TRUE)
show(mazmod_dfcp)
plot(mazmod_dfcp)
abline(v=c(6,17+59/60),lty=3)

# Mazama activity: dry forest, cattle absent
length(mazrad_dfca)
mazmod_dfca<-fitact(mazrad_dfca,reps=10000,sample="model",show=TRUE)
show(mazmod_dfca)
plot(mazmod_dfca)
abline(v=c(6,17+59/60),lty=3)

# Mazama activity: savanna, cattle present
length(mazrad_svcp)
mazmod_svcp<-fitact(mazrad_svcp,reps=10000,sample="data",show=TRUE);beep()
show(mazmod_svcp)
plot(mazmod_svcp)
abline(v=c(6,17+59/60),lty=3)

# Mazama activity: savanna, cattle absent
length(mazrad_svca)
mazmod_svca<-fitact(mazrad_svca,reps=10000,sample="model",show=TRUE);beep()
show(mazmod_svca)
plot(mazmod_svca)
abline(v=c(6,17+59/60),lty=3)

# Mazama activity: open savanna, cattle present
length(mazrad_osvcp)
mazmod_osvcp<-fitact(mazrad_osvcp,reps=100,sample="data",show=TRUE);beep()
show(mazmod_osvcp)
plot(mazmod_osvcp)
abline(v=c(6,17+59/60),lty=3)

# Mazama activity: open savanna, cattle absent
length(mazrad_osvca)
mazmod_osvca<-fitact(mazrad_osvca,reps=10000,sample="data",show=TRUE);beep()
show(mazmod_osvca)
plot(mazmod_osvca)
abline(v=c(6,17+59/60),lty=3)

#COMPARE ACTIVITY LEVELS AND OVERLAP COEFFICIENTS

#Preparation
bsmazrad_dfcp<-resample(mazrad_dfcp,10000)
bsmazrad_dfca<-resample(mazrad_dfca,10000)
bsmazrad_svcp<-resample(mazrad_svcp,10000)
bsmazrad_svca<-resample(mazrad_svca,10000)
bsmazrad_osvcp<-resample(mazrad_osvcp,10000)
bsmazrad_osvca<-resample(mazrad_osvca,10000)
dim(bsmazrad_dfcp)
dim(bsmazrad_dfca)
dim(bsmazrad_svcp)
dim(bsmazrad_svca)
dim(bsmazrad_osvcp)
dim(bsmazrad_osvca)

#Compare cattle present with cattle absent: dry forest
compareAct(list(mazmod_dfcp,mazmod_svcp))

#Overlap
overlapPlot(mazrad_dfcp,mazrad_dfca, extend='lightsteelblue4',rug=TRUE)
legend('topright', c("Cattle present","Cattle absent"),lty=c(1,2),y.intersp=1,col=c("black","blue"),bty='white')
abline(v=c(6,17+59/60),lty=3)
(dhats_dfcpca<-overlapEst(mazrad_dfcp,mazrad_dfca))
bsout_dfcpca<-bootEst(bsmazrad_dfcp,bsmazrad_dfca,adjust=c(0.8,NA,NA))
dim(bsout_dfcpca)
colMeans(bsout_dfcpca)
bs_dfcpca<-as.vector(bsout_dfcpca[,1])
bootCI(dhats_dfcpca[1],bs_dfcpca)['basic0',]
watson.two.test(mazrad_dfcp,mazrad_dfca, alpha=0.05)
#Compare cattle present with cattle absent: savanna
compareAct(list(mazmod_svcp,mazmod_svca))

#Overlap
overlapPlot(mazrad_svcp,mazrad_svca, extend='lightsteelblue4',rug=TRUE)
legend('topright', c("Cattle present","Cattle absent"),lty=c(1,2),y.intersp=1,col=c("black","blue"),bty='white')
abline(v=c(6,17+59/60),lty=3)
length(mazrad_svcp)
length(mazrad_svca)
(dhats_svcpca<-overlapEst(mazrad_svcp,mazrad_svca))
bsout_svcpca<-bootEst(bsmazrad_svcp,bsmazrad_svca,adjust=c(NA,1,NA))
dim(bsout_svcpca)
colMeans(bsout_svcpca)
bs_svcpca<-as.vector(bsout_svcpca[,2])
bootCI(dhats_svcpca[2],bs_svcpca)['basic0',]
watson.two.test(mazrad_svcp,mazrad_svca)
#Compare cattle present with cattle absent: opensavanna
compareAct(list(mazmod_osvcp,mazmod_osvca))

#Overlap
overlapPlot(mazrad_osvcp,mazrad_osvca, extend='lightsteelblue4',rug=TRUE)
legend('topright', c("Cattle present","Cattle absent"),lty=c(1,2),y.intersp=2,col=c("black","blue"),bty='white')
abline(v=c(6,17+59/60),lty=3)
length(mazrad_osvcp)
length(mazrad_osvca)
(dhats_osvcpca<-overlapEst(mazrad_osvcp,mazrad_osvca))
bsout_osvcpca<-bootEst(bsmazrad_osvcp,bsmazrad_osvca,adjust=c(NA,1,NA))
dim(bsout_osvcpca)
colMeans(bsout_osvcpca)
bs_osvcpca<-as.vector(bsout_osvcpca[,2])
bootCI(dhats_osvcpca[2],bs_osvcpca)['basic0',]
watson.two.test(mazrad_osvcp,mazrad_osvca)
###OVERALL OVERLAP BETWEEN DEER AND CATTLE###

#Preparation

bosrad<-radianbos[bos]
bosrad_sv<-radianbos[bossv]
bosrad_osv<-radianbos[bososv]
bosrad_df<-radianbos[bosdf]

#All habitats
overlapPlot(mazrad_cp,bosrad,extend='lightsteelblue4') #use mazrad_cp so that only sites where cattle are present are included
legend('topleft',c("Mazama","Bos"),lty=c(1,1),col=c("black","blue"),bty='n')
length(bosrad)
(dhats_mazbos<-overlapEst(mazrad_cp,bosrad))
bsbosrad<-resample(bosrad,10000)
dim(bsbosrad)
bsout_mazbos<-bootEst(bsmazrad_cp,bsbosrad,adjust=c(NA,1,NA))
colMeans(bsout_mazbos)
bs_mazbos<-as.vector(bsout_mazbos[,2])
bootCI(dhats_mazbos[2],bs_mazbos)['norm0',]

#Dry forest
overlapPlot(mazrad_dfcp,bosrad_df,extend='lightsteelblue4',rug=TRUE) 
legend('topright',c("Mazama","Bos"),lty=c(1,2),y.intersp=1,col=c("black","blue"),bty='white')
 length(bosrad_df)
length(mazrad_dfcp)
(dhats_mazbos_df<-overlapEst(mazrad_dfcp,bosrad_df))
bsbosrad_df<-resample(bosrad_df,10000)
dim(bsbosrad_df)
bsout_mazbos_df<-bootEst(bsmazrad_dfcp,bsbosrad_df,adjust=c(1,NA,NA))

colMeans(bsout_mazbos_df)
bs_mazbos_df<-as.vector(bsout_mazbos_df[,1])
bootCI(dhats_mazbos_df[1],bs_mazbos_df)['basic0',]
watson.two.test(mazrad_dfcp,bosrad_df)
#Savanna
overlapPlot(mazrad_svcp,bosrad_sv,extend='lightsteelblue4',rug=TRUE) 
legend('topright',c("Mazama","Bos"),lty=c(1,2),y.intersp=1,col=c("black","blue"),bty='white')
length(bosrad_sv)
length(mazrad_svcp)
(dhats_mazbos_sv<-overlapEst(mazrad_svcp,bosrad_sv))
bsbosrad_sv<-resample(bosrad_sv,10000)
dim(bsbosrad_sv)
bsout_mazbos_sv<-bootEst(bsmazrad_svcp,bsbosrad_sv,adjust=c(NA,1,NA))
colMeans(bsout_mazbos_sv)
bs_mazbos_sv<-as.vector(bsout_mazbos_sv[,2])
bootCI(dhats_mazbos_sv[2],bs_mazbos_sv)['basic0',]
watson.two.test(mazrad_svcp,bosrad_sv)

#Open Savanna
overlapPlot(mazrad_osvcp,bosrad_osv,extend='lightsteelblue4',rug=T) 
legend('topright',c("Mazama","Bos"),lty=c(1,1),col=c("black","blue"),bty='white')
length(bosrad_osv)
length(mazrad_osvcp)
(dhats_mazbos_osv<-overlapEst(mazrad_osvcp,bosrad_osv))
bsbosrad_osv<-resample(bosrad_osv,10000)
dim(bsbosrad_osv)
bsout_mazbos_osv<-bootEst(bsmazrad_osvcp,bsbosrad_osv,adjust=c(NA,1,NA))
colMeans(bsout_mazbos_osv)
bs_mazbos_osv<-as.vector(bsout_mazbos_osv[,2])
bootCI(dhats_mazbos_osv[2],bs_mazbos_osv)['basic0',]
watson.two.test(mazrad_osvcp,bosrad_osv)
