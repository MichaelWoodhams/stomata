# Derived from file 'bmMultipleSampleTimes nick version.R'
# Kept simulation code, removed analysis code.
# Modified to take bdratio and alpha on command line.
# (This will allow simulations to be run in parallel.)
# removed 'proc' and 'procs' variables.
library(ape)
library(TreeSim)
library(fasttree2019)

# Args: <measure> <bdratio> <climate> [<reps>] [<first rep>]

measure.stats <- read.csv("OU_parameters_full_data.csv")
measurement.names <- substr(measure.stats[seq(2,10,2),1],7,100) # "GCL"   "SI"    "SD"    "EA"    "Gsmax"
measurement.sigma <- sqrt(measure.stats[seq(1,9,2),"mean"])
measurement.alpha <-      measure.stats[seq(2,10,2),"mean"]

args <- commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  print("No command line args - using debugging defaults")
  args <- c("GCL",1.5,1,10)
}

measure_name <- args[[1]]
measurement <- which(measure_name==measurement.names)
if (length(measurement)!=1) {
  stop(paste(c("First arg must be one of",measurement.names)))
}
bdratio <- as.numeric(args[[2]])
REPS = 10000
REP.START = 1
climate = as.numeric(args[[3]])
if (! climate %in% 1:6) {
  stop(paste(c("Climate must be in 1:6",climate)))
}
if (length(args)>=4) {
  REPS=as.numeric(args[[4]])
}
if (length(args)>=5) {
  REP.START=as.numeric(args[[5]])
}

sigma <- measurement.sigma[measurement]
alpha <- measurement.alpha[measurement]

root.value = 0 #starting value for the trait arbitrarily set to zero
numtips = 1700 #this is set to match what is known about Proteaceae

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

print(sprintf("Simulating measurement %s with sigma=%f, alpha=%f. BD ratio=%f, climate=%d, reps %d..%d",
              measure_name,sigma,alpha,bdratio,climate,REP.START,REPS))
print("Fossil sample times:")
print(fossil.ages)

#for (alpha in alpha_values){
#  for (bdratio in bdratio_values){
    print(sprintf("BD ratio=%f, sigma=%f, alpha=%f",bdratio,sigma,alpha))

    rand.correlations = matrix(0, REPS, 5)
    rand.sds = matrix(0, REPS, 5)
    rand.sigma2 = numeric(REPS)
    rand.ace.sigma2 = numeric(REPS)
    rand.ace.se = numeric(REPS)
    tree.depth = numeric(REPS)

    for (reps in REP.START:REPS){
      # RNG is used for initial trait value, and for sampling from fossil trait values
      # at a given fossil age.
      set.seed(reps+100*bdratio+10000*alpha)
      if (reps%%100 == 0) print(reps)
      if (alpha!=0) {
        # For OU process, equilibrium distribution of trait is Gaussian with mean 0 and variance sigma^2/(2 alpha)
        startTrait = rnorm(1,mean=0,sd=sigma/sqrt(2*alpha))
      } else {
        # alpha=0 is Brownian motion case for trait. There is no equilibrium distribution.
        startTrait = 0
      }
      #print(reps)
      res = fasttree(n = numtips, lambda = bdratio, mu = 1, alpha = alpha, frac = 1,
                     sigma = sigma, sampleTimes = fossil.ages / 93, traits = TRUE, seed=reps)
      # output names: tree, samples, n.leaves, sampleTimes
      tree.depth[reps] <- max(node.depth.edgelength(res$tree))

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

    results <- data.frame(rep=1:reps,
                          depth=tree.depth,
                          rand.corr1 = rand.correlations[,1],
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
    filename = paste("results/results", measure_name, bdratio ,REPS, "climate", climate, ".txt", sep="_")
    write.table(results, file = filename, quote=FALSE, row.names=FALSE)
#  }
#}

