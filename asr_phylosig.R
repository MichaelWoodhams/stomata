#install.packages(c('caper', 'ape', 'picante', 'nlme', 'lmodel2', 'geiger', 'phytools', 'doBy'))
library("caper")
library("ape")
library("picante")
library("nlme")
library("lmodel2")#, lib.loc="~/R/win-library/3.1")
library("geiger")
library("phytools")#, lib.loc="~/R/win-library/3.1")
library(doBy)


#read data and transform variables
#note, the data has multiple representatives for several tips 
#i.e. the tips of the tree are [mostly] genera, and we have multiple species
#for some of these tips

size = read.csv("size.csv")
size$GCL=log(size$GCL)
size$SI=(size$stom_index)
size$EA=log(size$epid_areas)
size$SD=log(size$stom_dens)
size$GS=log(size$gsmax)
tree = read.nexus("tree.nex")

#select one data set
size1=size
include=size[c("terminal","order1")]
rand1=runif(nrow(include))
include$rand=rand1
include=include[order(include$terminal,include$rand),]
size1$include=include$order1

#calculate phylogenetic signal for each of the variables
size1=size[which(size1$include==1),]
phylosig(tree,size1$GCL,method="K",test=TRUE,nsim=10000)
phylosig(tree,size1$SI,method="K",test=TRUE,nsim=10000)
phylosig(tree,size1$EA,method="K",test=TRUE,nsim=10000)
phylosig(tree,size1$SD,method="K",test=TRUE,nsim=10000)
phylosig(tree,size1$GS,method="K",test=TRUE,nsim=10000)

# calculate phylogenetic signal for each of the variables, 
# but do it for each of the alternative sets of data 
k=size[1,]
k$k_GCL=phylosig(tree,size1$GCL,method="K",test=FALSE)
k$k_SI=phylosig(tree,size1$SI,method="K",test=FALSE)
k$k_EA=phylosig(tree,size1$EA,method="K",test=FALSE)
k$k_SD=phylosig(tree,size1$SD,method="K",test=FALSE)
k$k_GS=phylosig(tree,size1$GS,method="K",test=FALSE)
k=k[c("k_GCL","k_SI","k_EA","k_SD","k_GS")]
k1=k[NULL]

for(i in 1:100){
  size1=size
  include=size[c("terminal","order1")]
  rand1=runif(nrow(include))
  include$rand=rand1
  include=include[order(include$terminal,include$rand),]
  size1$include=include$order1
  size1=size1[which(size1$include==1),]

  k1$k_GCL=phylosig(tree,size1$GCL,method="K",test=FALSE)
  k1$k_SI=phylosig(tree,size1$SI,method="K",test=FALSE)
  k1$k_EA=phylosig(tree,size1$EA,method="K",test=FALSE)
  k1$k_SD=phylosig(tree,size1$SD,method="K",test=FALSE)
  k1$k_GS=phylosig(tree,size1$GS,method="K",test=FALSE)
  k=rbind(k,k1)
}

k$not=1
summaryBy(k_GCL+k_SI+k_SD+k_EA+k_GS~not,data=k,FUN=c(mean,sd))

#Ancestral state analysis, including the estimation of sigma 
#(I think the function is wrong in that it calls it sigma2, but we should check that)
data=size[ which(size$order1==1), ]
gcl=data$GCL
names(gcl)=data$terminal
ace_gcl=ace(gcl,tree,type="continuous")
ace_gcl$sigma2
node_ages=branching.times(tree)
node_ages1=t(node_ages)
anc_gcl=ace_gcl$ace
anc_gcl1=rbind(node_ages1,anc_gcl)

#plot the ASR onto the phylogeny (n.b. this actually recalculates the ASR)
contMap(tree, gcl, res=100, fsize=NULL, ftype=NULL, lwd=4, legend=NULL, lims=NULL, outline=TRUE, sig=3)

rownames(size1) <-  paste(size1$terminal, 1:nrow(size1), sep='_') # size1$terminal

#Do it all again for SI
data=size[ which(size$order1==1), ]
SI=data$SI
names(SI)=data$terminal
ace_SI=ace(SI,tree,type="continuous")
ace_SI$sigma2
node_ages=branching.times(tree)
node_ages1=t(node_ages)
anc_SI=ace_SI$ace
anc_SI1=rbind(node_ages1,anc_SI)

contMap(tree, SI, res=100, fsize=NULL, ftype=NULL, lwd=4, legend=NULL, lims=NULL, outline=TRUE, sig=3)
rownames(data) <- paste(data$terminal, 1:nrow(data), sep='_') # data$terminal
phylosig(tree,data$SI,method="K",test=T)


#Do it all again for epidermal cell size
data=size[ which(size$order1==1), ]
EA=data$EA
names(EA)=data$terminal
ace_EA=ace(EA,tree,type="continuous")
ace_EA$sigma2
node_ages=branching.times(tree)
node_ages1=t(node_ages)
anc_EA=ace_EA$ace
anc_EA1=rbind(node_ages1,anc_EA)

contMap(tree, EA, res=100, fsize=NULL, ftype=NULL, lwd=4, legend=NULL, lims=NULL, outline=TRUE, sig=3)
rownames(data) <- paste(data$terminal, 1:nrow(data), sep='_') # data$terminal
phylosig(tree,data$EA,method="K",test=T)

#Do it all again for stomatal density
data=size[ which(size$order1==1), ]
SD=data$SD
names(SD)=data$terminal
ace_SD=ace(SD,tree,type="continuous")
ace_SD$sigma2
node_ages=branching.times(tree)
node_ages1=t(node_ages)
anc_SD=ace_SD$ace
anc_SD1=rbind(node_ages1,anc_SD)

contMap(tree, SD, res=100, fsize=NULL, ftype=NULL, lwd=4, legend=NULL, lims=NULL, outline=TRUE, sig=3)
rownames(data) <- paste(data$terminal, 1:nrow(data), sep='_') # data$terminal
phylosig(tree,data$SD,method="K",test=T)

#Do it all again for gsmax
data=size[ which(size$order1==1), ]
GS=data$GS
names(GS)=data$terminal
ace_GS=ace(GS,tree,type="continuous")
ace_GS$sigma2
node_ages=branching.times(tree)
node_ages1=t(node_ages)
anc_GS=ace_GS$ace
anc_GS1=rbind(node_ages1,anc_GS)

contMap(tree, GS, res=100, fsize=NULL, ftype=NULL, lwd=4, legend=NULL, lims=NULL, outline=TRUE, sig=3)
rownames(data) <- paste(data$terminal, 1:nrow(data), sep='_') # data$terminal
phylosig(tree,data$GS,method="K",test=T)

#calculate sigma for standardised tree
tree1=tree
tree1$edge.length=tree1$edge.length/93.2

data=size[ which(size$order1==1), ]
sd_gcl=sd(data$GCL)
gcl=data$GCL/sd_gcl
names(gcl)=data$terminal
ace_gcl=ace(gcl,tree1,type="continuous")
ace_gcl$sigma2

sd_SI=sd(data$SI)
SI=data$SI/sd_SI
names(SI)=data$terminal
ace_SI=ace(SI,tree1,type="continuous")
ace_SI$sigma2

sd_EA=sd(data$EA)
EA=data$EA/sd_EA
names(EA)=data$terminal
ace_EA=ace(EA,tree1,type="continuous")
ace_EA$sigma2

sd_SD=sd(data$SD)
SD=data$SD/sd_SD
names(SD)=data$terminal
ace_SD=ace(SD,tree1,type="continuous")
ace_SD$sigma2

sd_GS=sd(data$GS)
GS=data$GS/sd_GS
names(GS)=data$terminal
ace_GS=ace(GS,tree1,type="continuous")
ace_GS$sigma2

#check how sigma2 is affected by scaling tree depth
tree1=tree
tree1$edge.length=tree1$edge.length/9.32
sd_gcl=sd(data$GCL)
gcl=data$GCL/sd_gcl
names(gcl)=data$terminal
ace_gcl=ace(gcl,tree1,type="continuous")
ace_gcl$sigma2

tree1=tree
tree1$edge.length=tree1$edge.length/0.932
sd_gcl=sd(data$GCL)
gcl=data$GCL/sd_gcl
names(gcl)=data$terminal
ace_gcl=ace(gcl,tree1,type="continuous")
ace_gcl$sigma2



#output results
anc_gcl2=t(anc_gcl1)
write.table(anc_gcl2,"gcl_anc.txt", append=FALSE, sep="\t")

anc_SI2=t(anc_SI1)
write.table(anc_SI2,"SI_anc.txt", append=FALSE, sep="\t")

anc_EA2=t(anc_EA1)
write.table(anc_EA2,"EA_anc.txt", append=FALSE, sep="\t")

anc_SD2=t(anc_SD1)
write.table(anc_SD2,"SD_anc.txt", append=FALSE, sep="\t")

anc_GS2=t(anc_GS1)
write.table(anc_GS2,"GS_anc.txt", append=FALSE, sep="\t")
