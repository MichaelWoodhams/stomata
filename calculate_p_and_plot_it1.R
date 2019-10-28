library(plyr)
library(doBy)
library(data.table)
#setwd("/home/mdw2/Documents/stomata/code")

#SI then loggcl, then Log_epi, then logSD, then logGmax

correls=data.frame(c(rep(1,5),rep(2,5),rep(4,5)),
                   c(rep("Beerling.CO2",5),rep("Anag.CO2",5),rep("Hansen",5)),
                   (rep(c("SI","GCL","EA","SD","Gsmax"),3)),
                   (c(-0.143083,-0.201605,-0.039522,0.021195,-0.072471,0.196907,-0.424963,
                      -0.06725,0.215795,0.078833,0.22372,-0.49242,-0.1312,0.3369,0.1726)),
                   (rep(c(2.4308543, 1.18148, 1.14046, 0.69834, 1.71374),3)))
names(correls)=c("climate","climate_name","measure","r","real.sigma2")

p <- read.table("combined.txt",header=TRUE)

#create function to calculate probabilities for r

probs=function(measure1,climate1){
  p1=p[which(p$measure==measure1 & p$climate==climate1),]
  c1=correls[which(correls$measure==measure1 &correls$climate==climate1),]
  px1=summaryBy(rand.corr1~measure+bdratio+climate,data=p1,FUN=length)
  pv=data.frame(measure=character(),bdratio=double(),climate=double(),p=double())
  
  for(var in 6:10){
    p2=p1[which(p1[var]>c1$r),]
    p2=p2[c(1,2,3,var)]
    names(p2)=c("measure","bdratio","climate","rand.corr")
    py1=summaryBy(rand.corr~measure+bdratio+climate,data=p2,FUN=length)
    py1$p=1-(py1[[4]]/px1[[4]])
    py1=py1[-4]
    pv=rbind(pv,py1)
  }
return(pv)  
}

#calculate probabilities for each scenario, could be looped but I couldn't be bothered
x=probs("EA",1)
x1=rbind(x,probs("EA",2))
x1=rbind(x1,probs("EA",4))
x1=rbind(x1,probs("GCL",1))
x1=rbind(x1,probs("GCL",2))
x1=rbind(x1,probs("GCL",4))
x1=rbind(x1,probs("Gsmax",1))
x1=rbind(x1,probs("Gsmax",2))
x1=rbind(x1,probs("Gsmax",4))
x1=rbind(x1,probs("SD",1))
x1=rbind(x1,probs("SD",2))
x1=rbind(x1,probs("SD",4))
x1=rbind(x1,probs("SI",1))
x1=rbind(x1,probs("SI",2))
x1=rbind(x1,probs("SI",4))


pvalues=summaryBy(p~measure+bdratio+climate,data=x1,FUN=mean)
pvalues$climate_name=ifelse(pvalues$climate==1,"Beerling.CO2",pvalues$climate)
pvalues$climate_name=ifelse(pvalues$climate_name==2,"Anag.CO2",pvalues$climate_name)
pvalues$climate_name=ifelse(pvalues$climate_name==4,"Hansen",pvalues$climate_name)

pvalues$mu_on_lambda=1/pvalues$bdratio

#write.csv(p,"new.probs.csv")



#create function to plot the things  
plotz=function(yvar){
  plot(x$mu_on_lambda,yvar,ylim=plims,pch=16,ylab="",xlab="",cex=0.8,col="red",yaxt="n",xaxt="n",bg=NULL)
  title(ylab="probability",line=2.2,cex=2)
  title(xlab=expression(paste(mu,"/",lambda)),line=0.8,cex=2)
  axis(side=1,at=c(0.0,0.2,0.4,0.6,0.8),labels=c("0.0","0.2","0.4","0.6","0.8"),
       las=1,mgp=c(2,0.1,0),tck=-.015,cex=2,bg=NULL)
  axis(side=2,at=ats,labels=labs,las=2,mgp=c(2,0.4,0),tck=-.015,cex=2,bg=NULL)
  #text(x=0.15,y=0.2,"log of guard cell length")
  fitP=lm(yvar~x$mu_on_lambda)
  abline(fitP,col="red",bg=NULL)
  abline(h=log10(0.05+0.0001),lty=2,bg=NULL)
  abline(h=log10(0.01+0.0001),lty=3,bg=NULL)
  abline(h=log10(0.001+0.0001),lty=4,bg=NULL)
}

#set layout, margin sizes, axis labels and plot limits

#little fudge to deal with zero values; this is corrected in the axis labelling
pvalues$prob=log10(pvalues$p.mean+0.00005)

#define plotting parameters
ats=c(log10(0.0001+0.0001),log10(0.001+0.0001),log10(0.01+0.0001),log10(0.05+0.0001),log10(1+0.0001))
labs=c("0.0001","0.001","0.01","0.05","1")
plims=c(-4,0)

#do the individual plots
#plot significance of correlations epidermal characters against CO2 and temperature
svg(filename = "tests_of_correlations_log scale_newest.svg", width=7.8, height = 7.5)
par(mfcol=c(5,3),mai=c(0.25,0.45,0.1,0.1),pty="m")

x=pvalues[which(pvalues$climate_name=="Hansen" & pvalues$measure=="GCL"),]
plotz(x$prob)
x=pvalues[which(pvalues$climate_name=="Hansen" & pvalues$measure=="SD"),]
plotz(x$prob)
x=pvalues[which(pvalues$climate_name=="Hansen" & pvalues$measure=="EA"),]
plotz(x$prob)
x=pvalues[which(pvalues$climate_name=="Hansen" & pvalues$measure=="Gsmax"),]
plotz(x$prob)
x=pvalues[which(pvalues$climate_name=="Hansen" & pvalues$measure=="SI"),]
plotz(x$prob)
x=pvalues[which(pvalues$climate_name=="Beerling.CO2" & pvalues$measure=="GCL"),]
plotz(x$prob)
x=pvalues[which(pvalues$climate_name=="Beerling.CO2" & pvalues$measure=="SD"),]
plotz(x$prob)
x=pvalues[which(pvalues$climate_name=="Beerling.CO2" & pvalues$measure=="EA"),]
plotz(x$prob)
x=pvalues[which(pvalues$climate_name=="Beerling.CO2" & pvalues$measure=="Gsmax"),]
plotz(x$prob)
x=pvalues[which(pvalues$climate_name=="Beerling.CO2" & pvalues$measure=="SI"),]
plotz(x$prob)
x=pvalues[which(pvalues$climate_name=="Anag.CO2" & pvalues$measure=="GCL"),]
plotz(x$prob)
x=pvalues[which(pvalues$climate_name=="Anag.CO2" & pvalues$measure=="SD"),]
plotz(x$prob)
x=pvalues[which(pvalues$climate_name=="Anag.CO2" & pvalues$measure=="EA"),]
plotz(x$prob)
x=pvalues[which(pvalues$climate_name=="Anag.CO2" & pvalues$measure=="Gsmax"),]
plotz(x$prob)
x=pvalues[which(pvalues$climate_name=="Anag.CO2" & pvalues$measure=="SI"),]
plotz(x$prob)
dev.off()

