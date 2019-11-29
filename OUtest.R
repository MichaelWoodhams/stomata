# Compare our OU simulations with those of phylolm
# command line usage: measure_name bdratio numsims comparison_per_sim
# Where measure_name is one of  "GCL"   "SI"    "SD"    "EA"    "Gsmax"
#     bdratio = birth/death ratio
#     numsims = number of times to run our sim (defaults to 100)
#     comparisons_per_sim = number of times to run phylolm sim for each of our sims
#         (defaults to 5)

library(ape)
library(phylolm)
library(TreeSim)
library(fasttree2019)

measure.stats <- read.csv("OU_parameters_full_data.csv")
measurement.names <- substr(measure.stats[seq(2,10,2),1],7,100) # "GCL"   "SI"    "SD"    "EA"    "Gsmax"
measurement.sigma <- sqrt(measure.stats[seq(1,9,2),"mean"])
measurement.alpha <-      measure.stats[seq(2,10,2),"mean"]

args <- commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  print("No command line args - using debugging defaults")
  args <- c("GCL",1.5,10,5)
}

measure_name <- args[[1]]
measurement <- which(measure_name==measurement.names)
if (length(measurement)!=1) {
  stop(paste(c("First arg must be one of",measurement.names)))
}
bdratio <- as.numeric(args[[2]])
REPS = 100
REP.START = 1
if (length(args)>=3) {
  REPS=as.numeric(args[[3]])
}
COMPS = 5
if (length(args)>=4) {
  COMPS=as.numeric(args[[4]])
}
if (length(args)>=5) {
  REP.START=as.numeric(args[[5]])
}

sigma <- measurement.sigma[measurement]
alpha <- measurement.alpha[measurement]

print(sprintf("Simulating measurement %s with sigma=%f, alpha=%f",measure_name,sigma,alpha))

n <- 1700; lambda <- bdratio; mu <- 1; sampleTimes <- c(50); treeDepth <- 93; upper.bound <- 1
# Don't care about sample times, just give one time to avoid potential errors.
# TO DO: try with zero times.

results <- NULL
warnings <- data.frame(rep=integer(),comp=integer(),sigma=double(),alpha=double())
for (rep in REP.START:REPS) {
  res = fasttree(n = n, lambda = lambda, mu = mu, alpha = alpha, sigma = sigma,
                 sampleTimes = sampleTimes, treeDepth = treeDepth, traits = TRUE, seed=rep)
  # Trait measurements on the extant tips
  trait <- res$tree$edge.trait[n>=res$tree$edge[,2]]

  #restrict to extant tips only
  tree=keep.tip(res$tree,seq(1:1700))
  #phylolm requires no root edge or root edge of zero
  tree$root.edge=0

  traits <- data.frame(trait,row.names=tree$tip.label)
  OU=phylolm(trait ~ 1, data = traits, phy = tree, model = "OUrandomRoot", upper.bound = 1)
  OU.sigma=OU[[2]]^0.5
  OU.alpha=OU[[3]]

  row <- data.frame(rep,sigma=OU.sigma,alpha=OU.alpha)

  #comp.fits <- c()
  for (comp in 1:COMPS) {
    #simulate new trait on the tree using rTrait
    set.seed(rep*1000000+comp)
    new.trait=rTrait(n=1,tree,model="OU",parameters=list(optimal.value=0, sigma2=sigma^2, 3,alpha=alpha))
    #avoids warning regarding mismatched labels. Be careful that the order of the trait and the tree are the same
    traits$new.trait <- new.trait

    #calculates OU parameters
    #be careful to look for warnings. phylolm often gets lost and can't find the parameter.
    #print(sprintf("rep=%d, comp=%d",rep,comp))
    sigAlpha <- tryCatch(
      expr={
        OU.newtrait<-NULL
        OU.newtrait=phylolm(new.trait ~ 1, data = traits, phy = tree, model = "OUrandomRoot",upper.bound=1)
        c(OU.newtrait[[2]]^0.5,OU.newtrait.alpha=OU.newtrait[[3]])
      },
      warning = function(w){
        print(sprintf("Caught warning rep=%d,comp=%d",rep,comp))
        return(c(NA,NA))
      })

    row[paste("sigma",comp,sep="")] = sigAlpha[1]
    row[paste("alpha",comp,sep="")] = sigAlpha[2]
    # Debugging: if got a warning, what value do we get when we ignore the warning?
    if (is.na(sigAlpha[1])) {
      OU.newtrait=phylolm(new.trait ~ 1, data = traits, phy = tree, model = "OUrandomRoot")
      newwarning <- data.frame(rep,comp,OU.newtrait[[2]]^0.5,OU.newtrait.alpha=OU.newtrait[[3]])
      names(newwarning)<-names(warnings)
      warnings <- rbind(warnings,newwarning)
    }
    #comp.fits <- c(comp.fits,OU.newtrait.sigma,OU.newtrait.alpha)
    #print(comp.fits)
  }
  if (is.null(results)) {
    results <- row
  } else {
    results <- rbind(results,row)
  }
}

# For each sim, count how many rTrait values are less than the fasttree value.
# If all is good, this should be distributed uniformly in the range 0:COMPS.

results$orderFastSigma = -1
results$orderFastAlpha = -1
results$orderTrueSigma = -1
results$orderTrueAlpha = -1

for (rep in 1:REPS) {
  fastSigma <- results[rep,"sigma"]
  fastAlpha <- results[rep,"alpha"]
  sigmas = results[rep,seq(4,2+2*COMPS,2)]
  alphas = results[rep,seq(5,3+2*COMPS,2)]
  results$orderFastSigma[rep]<-sum(sigmas<fastSigma)
  results$orderFastAlpha[rep]<-sum(alphas<fastAlpha)
  results$orderTrueSigma[rep]<-sum(sigmas<sigma)
  results$orderTrueAlpha[rep]<-sum(alphas<alpha)
}

write.csv(results,file=sprintf("results/compare_%s_%f.csv",measure_name,bdratio))

# results <- read.csv(file="../results/compare_GCL_1.500000.csv")
print("Order fast alpha vs rTrait alpha")
print(table(results$orderFastAlpha))

print("Order fast sigma vs rTrait sigma")
print(table(results$orderFastSigma))

print("Order true alpha vs rTrait alpha")
print(table(results$orderTrueAlpha))

print("Order true sigma vs rTrait sigma")
print(table(results$orderTrueSigma))





