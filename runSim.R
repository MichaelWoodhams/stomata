setwd("/home/mdw2/Documents/botany/stomata/code")
# Derived from file 'bmMultipleSampleTimes nick version.R'
# Kept simulation code, removed analysis code.
# Modified to take bdratio and theta on command line.
# (This will allow simulations to be run in parallel.)
# removed 'proc' and 'procs' variables.
library(ape)
library(TreeSim)
source("fast.tree.r")
dyn.load('tree_climb.so')

#let's see if we can write a version of what BM option of rTraitCont does
#code here is adapted from rTraitCont and the C function that it calls
#set.seed(123)

climate = 6

# data rows: SI, log(gcl), log(Epi Area), log(Stom Dens), log(gmax)
# SI data gotten from modified Greg's ASR code (as now not logged and very slightly different)
# otherwise all data from sigma values.xlsx
# only other differences are factor of 100 in GCL and 2 in SD for some reason
real.sigma = c(2.4308543, 1.18148, 1.14046, 0.69834, 1.71374)
real.se = c(0.4080669, 0.19833, 0.19145, 0.11724, 0.28769)

REPS = 5000
root.value = 0 #starting value for the trait arbitrarily set to zero
numtips = 1700 #this is set to match what is known about Proteaceae

#bdratio_values = c(1.1, 1.3, 1.5, 2, 3, 5, 10, 20, 50, 100)
#theta_values = c(0,0.1,0.5,1)
args <- commandArgs(trailingOnly = TRUE)
bdratio <- as.numeric(args[[1]])
theta <- as.numeric(args[[2]])
if (length(args)>=3) {
  REPS=as.numeric(args[[3]])
}
if (length(args)>=4) {
  climate=as.numeric(args[[4]])
  print(climate)
}

sigma = 1
data = read.csv('zachos_correl_and_calc_gsmax1.csv')

clim.vec = data[,c(3:6, 5:6)[climate]]
ages = data[,2]

if (climate <= 4)
{
  keep = which(!is.na(clim.vec) & data$delete == 0)
} else {
  keep = which(!is.na(clim.vec))
}

clim.vec = clim.vec[keep]
ages = ages[keep]

#real sample times (assume sorted)
tab = table(ages)
fossil.ages <- as.numeric(rownames(tab)) #c(15, 18, 20, 24, 30, 32, 33, 34, 40, 43, 45, 47, 49, 50, 51, 56)
fossil.ages <- 93 - sort(fossil.ages, decreasing=TRUE)
numSamplePoints = length(fossil.ages)
#number of fossils at each time point
numfossils = as.integer(tab) #c(2, 2, 2, 13, 1 , 17, 4, 2, 20, 7, 13, 5, 4, 3, 2, 12) 
numfossils  <- numfossils[length(numfossils):1]

#for (theta in theta_values){
#  for (bdratio in bdratio_values){
    print(paste('Theta:',theta, '   BD RATIO', bdratio))
  
    rand.correlations = matrix(0, REPS, 5)
    rand.sds = matrix(0, REPS, 5)
    rand.sigma2 = numeric(REPS)
    rand.ace.sigma2 = numeric(REPS)
    rand.ace.se = numeric(REPS)
  
    for (reps in 1:REPS){
      # RNG is used for initial trait value, and for sampling from fossil trait values
      # at a given fossil age.
      set.seed(reps+100*bdratio+10000*theta)
      if (reps%%100 == 0) print(reps)
      if (theta!=0) {
        # For OU process, equilibrium distribution of trait is Gaussian with mean 0 and variance sigma^2/(2 theta)
        startTrait = rnorm(1,mean=0,sd=sigma/sqrt(2*theta))
      } else {
        # theta=0 is Brownian motion case for trait. There is no equilibrium distribution.
        startTrait = 0
      }
      #print(reps)
      res = fast.tree(n = numtips, lambda = bdratio, mu = 1, theta = theta, frac = 1, 
                      sigma = sigma, sampleTimes = fossil.ages / 93, traits = TRUE, seed=reps)
      # output names: tree, samples, n.leaves, sampleTimes
      
      node.trait = numeric(max(res$tree$edge[,2]))
      node.trait[res$tree$edge[,2]] = res$tree$edge.trait
      
      for (j in 1:5) # different data
      {
        data.vec = data[keep, 10 + j]
        ij = which(!is.na(data.vec))
        #real sample times (assume sorted)
        tab = table(ages[ij])
        
        # assume fossil ages will be the same!!!
        
        #      fossil.ages <- as.numeric(rownames(tab)) #c(15, 18, 20, 24, 30, 32, 33, 34, 40, 43, 45, 47, 49, 50, 51, 56)
        #fossil.ages <- 93 - sort(fossil.ages, decreasing=TRUE)
        #numSamplePoints = length(fossil.ages)
        #number of fossils at each time point
        numfossils = as.integer(tab) #c(2, 2, 2, 13, 1 , 17, 4, 2, 20, 7, 13, 5, 4, 3, 2, 12) 
        numfossils  <- numfossils[length(numfossils):1]
        
        
        #set up a list of lists to store the samples in
        df <- NULL #this data frame will store samples ages and trait values 
        
        for (i in 1:numSamplePoints) {
          theseSamples <- res$samples[1:res$n.leaves[i], i] #all the species at this time slice
          theseFossils <- sample(theseSamples, numfossils[i], replace=TRUE) #replacement or not?
          df <- rbind(df, cbind(theseFossils, res$sampleTimes[i]))
        }
        
        #data frame containing the correlation between random values and zachos
        #each row corresponds to an individual sample
        rand.correlations[reps, j] = cor(df[ ,1], clim.vec[ij])
        rand.sds[reps, j] = sd(df[,1])
        # normalised sigma = tree length / sd@tips
        #    scale = max(node.depth.edgelength(res$tree))/sd(node.trait[1:numtips])^2
        #    rand.sigma2[reps] = scale
        
      }
      
      # TRIM TREE FOR ASR
      res$tree$node.label = paste('n', 1:res$tree$Nnode, sep='')
      drop = (numtips+1):length(res$tree$tip.label)
      old.edge.label = c(res$tree$tip.label, res$tree$node.label)
      old.edge.label = old.edge.label[res$tree$edge[,2]]
      bush = drop.tip(res$tree, drop)
      new.edge.label = c(bush$tip.label, bush$node.label)
      new.edge.label = new.edge.label[bush$edge[,2]]
      all(old.edge.label[old.edge.label %in% new.edge.label] == new.edge.label)
      bush$edge.trait = res$tree$edge.trait[which(old.edge.label %in% new.edge.label)]
      ltt = ltt.plot.coords(bush, backward = FALSE)
      ii = min(which(ltt[,2] == 71))
      #plot(bush, show.tip.label = FALSE); axis(1); abline(v = ltt[ii,1], lty=2, col='grey')
      babies = phytools::treeSlice(bush, slice = mean(ltt[(ii-1):ii,1]), trivial = TRUE)
      tokill = NULL
      for (i in 1:71)
        tokill = c(tokill, setdiff(babies[[i]]$tip.label, sample(babies[[i]]$tip.label, 1)))
      old.edge.label = new.edge.label
      old.edge.trait = bush$edge.trait
      bush = drop.tip(bush, tokill)
      new.edge.label = c(bush$tip.label, bush$node.label)
      new.edge.label = new.edge.label[bush$edge[,2]]
      all(old.edge.label[old.edge.label %in% new.edge.label] == new.edge.label)
      bush$edge.trait = old.edge.trait[which(old.edge.label %in% new.edge.label)]
      node.trait = numeric(max(bush$edge[,2]))
      node.trait[bush$edge[,2]] = bush$edge.trait
      asr = ace(node.trait[1:(bush$Nnode+1)], bush, type="continuous")
      
      # normalised sigma = tree length / sd@tips
      scale = max(node.depth.edgelength(res$tree))/sd(node.trait[1:71])^2
      rand.sigma2[reps] = sigma * scale
      rand.ace.sigma2[reps] = asr$sigma2[1] * scale
      rand.ace.se[reps] = asr$sigma2[2] * scale
    }
    
    results <- data.frame(rand.corr1 = rand.correlations[,1],
                          rand.corr2 = rand.correlations[,2], 
                          rand.corr3 = rand.correlations[,3], 
                          rand.corr4 = rand.correlations[,4], 
                          rand.corr5 = rand.correlations[,5], 
                          rand.sd1 = rand.sds[,1], 
                          rand.sd2 = rand.sds[,2], 
                          rand.sd3 = rand.sds[,3], 
                          rand.sd4 = rand.sds[,4], 
                          rand.sd5 = rand.sds[,5], 
                          rand.sigma2 = rand.sigma2,
                          rand.ace.sigma2 = rand.ace.sigma2,
                          rand.ace.se = rand.ace.se)
    filename = paste("../results/results", theta, bdratio ,REPS, "climate", climate, ".txt", sep="_")
    write.table(results, file = filename, quote=FALSE, row.names=FALSE)
#  }
#}

