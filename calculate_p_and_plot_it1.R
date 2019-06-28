library(plyr)
setwd("R:/SET/PlantSci/Staff/Greg Jordan/in progress/prots/prot directional evolution/final simulation results")

#SI then loggcl, then Log_epi, then logSD, then logGmax

correls=data.frame(c(rep(1,5),rep(2,5),rep(4,5)),
                   c(rep("Beerling.CO2",5),rep("Anag.CO2",5),rep("Hansen",5)),
                   (rep(c("SI","gcl","EA","SD","GS"),3)),
                   (c(-0.143083,-0.201605,-0.039522,0.021195,-0.072471,0.196907,-0.424963,
                      -0.06725,0.215795,0.078833,0.22372,-0.49242,-0.1312,0.3369,0.1726)),
                   (rep(c(2.4308543, 1.18148, 1.14046, 0.69834, 1.71374),3)))
names(correls)=c("climate","climate_name","stom_par","r","real.sigma2")

#create function to calculate P

p=data.frame(mu_on_lambda=double(),nreps=double(),theta=double(),climate = double(),
             P_gcl=double(),P_EA=double(),P_SD=double(),P_SI=double(),P_GS=double(),
             P_s2_ACE=double())

files = list.files(pattern=".txt")
#filename=files[[3]]

for(i in 1:length(files)){
  j=i
  filename=files[[j]]
  p2= data.frame(as.numeric(sapply(strsplit(filename,"_"),"[[",2)))
  names(p2)="theta"
  p2$mu_on_lambda = 1/as.numeric(sapply(strsplit(filename,"_"),"[[",3))
  p2$nreps = as.numeric(sapply(strsplit(filename,"_"),"[[",4))
  climate=as.numeric(sapply(strsplit(filename,"_"),"[[",6))
  p2$climate = climate
  
  correls1=correls[which(correls$climate==climate),]
  
  p1= read.csv(filename, sep="")
  names(p1)=c("rand_SI","rand_gcl","rand_EA","rand_SD","rand_GS","sd_SI","sd_gcl","sd_EA","sd_SD","sd_GS",
              "rand.sigma2","rand.ace.sigma2","rand.ace.se")

  p0=p1[which(abs(p1$rand_SI)>abs(correls1[1,4])),]
  p2$P_SI=nrow(p0)/p2$nreps
  p0=p1[which(abs(p1$rand.ace.sigma2)>abs(correls1[1,5])),]
  p2$P_SI_sigma=nrow(p0)/p2$nreps
  
  p0=p1[which(abs(p1$rand_gcl)>abs(correls1[2,4])),]
  p2$P_gcl=nrow(p0)/p2$nreps
  p0=p1[which(abs(p1$rand.ace.sigma2)>abs(correls1[2,5])),]
  p2$P_gcl_sigma=nrow(p0)/p2$nreps
  
  p0=p1[which(abs(p1$rand_EA)>abs(correls1[3,4])),]
  p2$P_EA=nrow(p0)/p2$nreps
  p0=p1[which(abs(p1$rand.ace.sigma2)>abs(correls1[3,5])),]
  p2$P_EA_sigma=nrow(p0)/p2$nreps
  
  p0=p1[which(abs(p1$rand_SD)>abs(correls1[4,4])),]
  p2$P_SD=nrow(p0)/p2$nreps
  p0=p1[which(abs(p1$rand.ace.sigma2)>abs(correls1[4,5])),]
  p2$P_SD_sigma=nrow(p0)/p2$nreps
  
  p0=p1[which(abs(p1$rand_GS)>abs(correls1[5,4])),]
  p2$P_GS=nrow(p0)/p2$nreps
  p0=p1[which(abs(p1$rand.ace.sigma2)>abs(correls1[5,5])),]
  p2$P_GS_sigma=nrow(p0)/p2$nreps
  
  p=rbind(p,p2)
  p=p[order(p$mu_on_lambda),]
}

p$climate_name=ifelse(p$climate==1,"Beerling.CO2",p$climate)
p$climate_name=ifelse(p$climate_name==2,"Anag.CO2",p$climate_name)
p$climate_name=ifelse(p$climate_name==4,"Hansen",p$climate_name)

write.csv(p,"probs.csv")


p$P_gcl=log10(p$P_gcl+0.00005)
p$P_gcl_sigma=log10(p$P_gcl_sigma+0.00005)
p$P_SD=log10(p$P_SD+0.00005)
p$P_SD_sigma=log10(p$P_SD_sigma+0.00005)
p$P_EA=log10(p$P_EA+0.00005)
p$P_EA_sigma=log10(p$P_EA_sigma+0.00005)
p$P_GS=log10(p$P_GS+0.00005)
p$P_GS_sigma=log10(p$P_GS_sigma+0.00005)
p$P_SI=log10(p$P_SI+0.00005)
p$P_SI_sigma=log10(p$P_SI_sigma+0.00005)

plims=c(-4.1,0.2)

#create function to plot the things  
plotz=function(yvar,yvar1){
  plot(x$mu_on_lambda,yvar,ylim=plims,pch=16,ylab="",xlab="",cex=0.8,col="red",yaxt="n",xaxt="n",bg=NULL)
  title(ylab="probability",line=2.2,cex.lab=0.85)
  title(xlab="mu/lambda",line=0.8,cex.lab=0.85)
  axis(side=1,at=c(0.0,0.2,0.4,0.6,0.8),labels=c("0.0","0.2","0.4","0.6","0.8"),
       las=1,mgp=c(2,0.1,0),tck=-.015,cex.axis=0.85,bg=NULL)
  axis(side=2,at=ats,labels=labs,las=2,mgp=c(2,0.4,0),tck=-.015,cex.axis=0.85,bg=NULL)
  #text(x=0.15,y=0.2,"log of guard cell length")
  fitP=lm(yvar~x$mu_on_lambda)
  abline(fitP,col="red",bg=NULL)
  points(x$mu_on_lambda,yvar1,pch=17,cex=0.8,col="blue",bg=NULL)
  fitP_ACE=lm(yvar1~x$mu_on_lambda)
  abline(fitP_ACE,col="blue",bg=NULL)
  abline(h=log10(0.05+0.0001),lty=2,bg=NULL)
  abline(h=log10(0.01+0.0001),lty=3,bg=NULL)
  abline(h=log10(0.001+0.0001),lty=4,bg=NULL)
}

#set layout, margin sizes, axis labels and plot limits
ats=c(log10(0.0001+0.0001),log10(0.001+0.0001),log10(0.01+0.0001),log10(0.05+0.0001),log10(1+0.0001))
labs=c("0.0001","0.001","0.01","0.05","1")
plims=c(-4,0)

#do the individual plots
#plot significance of correlations epidermal characters against CO2 and temperature
svg(filename = "tests_of_correlations_log scale_new_theta_0.svg", width=7.8, height = 7.5)
par(mfcol=c(5,3),mai=c(0.25,0.45,0.1,0.1),pty="m")

p1=p[which(p$theta==0),]
x=p1[which(p1$climate_name=="Hansen"),]
plotz(x$P_gcl,x$P_gcl_sigma)
plotz(x$P_SD,x$P_SD_sigma)
plotz(x$P_EA,x$P_EA_sigma)
plotz(x$P_GS,x$P_GS_sigma)
plotz(x$P_SI,x$P_SI_sigma)
x=p1[which(p1$climate_name=="Beerling.CO2"),]
plotz(x$P_gcl,x$P_gcl_sigma)
plotz(x$P_SD,x$P_SD_sigma)
plotz(x$P_EA,x$P_EA_sigma)
plotz(x$P_GS,x$P_GS_sigma)
plotz(x$P_SI,x$P_SI_sigma)
x=p1[which(p1$climate_name=="Anag.CO2"),]
plotz(x$P_gcl,x$P_gcl_sigma)
plotz(x$P_SD,x$P_SD_sigma)
plotz(x$P_EA,x$P_EA_sigma)
plotz(x$P_GS,x$P_GS_sigma)
plotz(x$P_SI,x$P_SI_sigma)
dev.off()

svg(filename = "tests_of_correlations_log scale_new_theta_0.1.svg", width=7.8, height = 7.5)
par(mfcol=c(5,3),mai=c(0.25,0.45,0.1,0.1),pty="m")
p1=p[which(p$theta==0.1),]
x=p1[which(p1$climate_name=="Hansen"),]
plotz(x$P_gcl,x$P_gcl_sigma)
plotz(x$P_SD,x$P_SD_sigma)
plotz(x$P_EA,x$P_EA_sigma)
plotz(x$P_GS,x$P_GS_sigma)
plotz(x$P_SI,x$P_SI_sigma)
x=p1[which(p1$climate_name=="Beerling.CO2"),]
plotz(x$P_gcl,x$P_gcl_sigma)
plotz(x$P_SD,x$P_SD_sigma)
plotz(x$P_EA,x$P_EA_sigma)
plotz(x$P_GS,x$P_GS_sigma)
plotz(x$P_SI,x$P_SI_sigma)
x=p1[which(p1$climate_name=="Anag.CO2"),]
plotz(x$P_gcl,x$P_gcl_sigma)
plotz(x$P_SD,x$P_SD_sigma)
plotz(x$P_EA,x$P_EA_sigma)
plotz(x$P_GS,x$P_GS_sigma)
plotz(x$P_SI,x$P_SI_sigma)
dev.off()

svg(filename = "tests_of_correlations_log scale_new_theta_0.2.svg", width=7.8, height = 7.5)
par(mfcol=c(5,3),mai=c(0.25,0.45,0.1,0.1),pty="m")
p1=p[which(p$theta==0.2),]
x=p1[which(p1$climate_name=="Hansen"),]
plotz(x$P_gcl,x$P_gcl_sigma)
plotz(x$P_SD,x$P_SD_sigma)
plotz(x$P_EA,x$P_EA_sigma)
plotz(x$P_GS,x$P_GS_sigma)
plotz(x$P_SI,x$P_SI_sigma)
x=p1[which(p1$climate_name=="Beerling.CO2"),]
plotz(x$P_gcl,x$P_gcl_sigma)
plotz(x$P_SD,x$P_SD_sigma)
plotz(x$P_EA,x$P_EA_sigma)
plotz(x$P_GS,x$P_GS_sigma)
plotz(x$P_SI,x$P_SI_sigma)
x=p1[which(p1$climate_name=="Anag.CO2"),]
plotz(x$P_gcl,x$P_gcl_sigma)
plotz(x$P_SD,x$P_SD_sigma)
plotz(x$P_EA,x$P_EA_sigma)
plotz(x$P_GS,x$P_GS_sigma)
plotz(x$P_SI,x$P_SI_sigma)
dev.off()

svg(filename = "tests_of_correlations_log scale_new_theta_0.5.svg", width=7.8, height = 7.5)
par(mfcol=c(5,3),mai=c(0.25,0.45,0.1,0.1),pty="m")
p1=p[which(p$theta==0.5),]
x=p1[which(p1$climate_name=="Hansen"),]
plotz(x$P_gcl,x$P_gcl_sigma)
plotz(x$P_SD,x$P_SD_sigma)
plotz(x$P_EA,x$P_EA_sigma)
plotz(x$P_GS,x$P_GS_sigma)
plotz(x$P_SI,x$P_SI_sigma)
x=p1[which(p1$climate_name=="Beerling.CO2"),]
plotz(x$P_gcl,x$P_gcl_sigma)
plotz(x$P_SD,x$P_SD_sigma)
plotz(x$P_EA,x$P_EA_sigma)
plotz(x$P_GS,x$P_GS_sigma)
plotz(x$P_SI,x$P_SI_sigma)
x=p1[which(p1$climate_name=="Anag.CO2"),]
plotz(x$P_gcl,x$P_gcl_sigma)
plotz(x$P_SD,x$P_SD_sigma)
plotz(x$P_EA,x$P_EA_sigma)
plotz(x$P_GS,x$P_GS_sigma)
plotz(x$P_SI,x$P_SI_sigma)
dev.off()

svg(filename = "tests_of_correlations_log scale_new_theta_1.svg", width=7.8, height = 7.5)
par(mfcol=c(5,3),mai=c(0.25,0.45,0.1,0.1),pty="m")
p1=p[which(p$theta==1),]
x=p1[which(p1$climate_name=="Hansen"),]
plotz(x$P_gcl,x$P_gcl_sigma)
plotz(x$P_SD,x$P_SD_sigma)
plotz(x$P_EA,x$P_EA_sigma)
plotz(x$P_GS,x$P_GS_sigma)
plotz(x$P_SI,x$P_SI_sigma)
x=p1[which(p1$climate_name=="Beerling.CO2"),]
plotz(x$P_gcl,x$P_gcl_sigma)
plotz(x$P_SD,x$P_SD_sigma)
plotz(x$P_EA,x$P_EA_sigma)
plotz(x$P_GS,x$P_GS_sigma)
plotz(x$P_SI,x$P_SI_sigma)
x=p1[which(p1$climate_name=="Anag.CO2"),]
plotz(x$P_gcl,x$P_gcl_sigma)
plotz(x$P_SD,x$P_SD_sigma)
plotz(x$P_EA,x$P_EA_sigma)
plotz(x$P_GS,x$P_GS_sigma)
plotz(x$P_SI,x$P_SI_sigma)
dev.off()

