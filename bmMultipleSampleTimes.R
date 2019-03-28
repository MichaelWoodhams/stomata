library(ape)
library(TreeSim)
setwd("/home/mdw2/Documents/botany/stomata/code")
dyn.load('tree_climb.so')

#let's see if we can write a version of what BM option of rTraitCont does
#code here is adapted from rTraitCont and the C function that it calls
set.seed(123)

REPS = 1000
#REPS = 10 # debugging
root.value = 0 #starting value for the trait arbitrarily set to zero
numtips = 1700 #this is set to match what is known about Proteaceae
#might want to play with next two
bdratio_values = c(1.1, 1.5, 2, 5, 10, 20, 50, 100)
theta_values = c(0,0.1,0.5,1)
sigma = 1


#the global temperature estimate
zachos <- c(2.744, 2.744, 2.572, 2.572, 2.495, 2.495, 2.763, 2.763, 2.763, 2.763, 2.763, 2.763, 2.763, 2.763, 2.763, 2.763, 2.763, 2.763, 2.763, 1.998, 2.113, 2.113, 2.113, 2.113, 2.113, 2.113, 2.113, 2.113, 2.113, 2.113, 2.113, 2.113, 2.113, 2.113, 2.113, 2.113, 2.113, 1.903, 1.903, 1.903, 1.903, 2.419, 2.419, 3.126, 3.126, 3.126, 3.126, 3.126, 3.126, 3.126, 3.126, 3.126, 3.126, 3.126, 3.126, 3.126, 3.126, 3.126, 3.126, 3.126, 3.126, 3.126, 3.126, 3.384, 3.384, 3.384, 3.384, 3.384, 3.384, 3.384, 3.499, 3.499, 3.499, 3.499, 3.499, 3.499, 3.499, 3.499, 3.499, 3.499, 3.499, 3.499, 3.499, 4.159, 4.159, 4.159, 4.159, 4.159, 4.14, 4.14, 4.14, 4.14, 4.197, 4.197, 4.197, 4.235, 4.235, 3.748, 3.748, 3.748, 3.748, 3.748, 3.748, 3.748, 3.748, 3.748, 3.748, 3.748, 3.748)
gcl <- c(28.2, 28.9, 26.1, 22.7, 31.5, 17.4, 27.5, 24, 33.4, 23.3, 53.7, 74.3, 30.7, 32.3, 30.5, 30.8, 50.4, 19.3, 24.5, 26, 20.2, 37.7, 45.3, 26.9, 28.7, 44.6, 35.2, 66.5, 39, 34.8, 25, 28.3, 32.6, 26.9, 54.3, 57.7, 23.4, 38.86371615, 39.07233467, 35.59714843, 48.04071378, 23.2, 28.2, 23.1, 22.3, 24.9, 18.6, 25, 27.4, 19.9, 26.9, 23.9, 25.2, 27.1, 23.7, 53, 27.9, 27.9, 26.7, 26.6, 27.8, 30.7, 24.6, 20.6, 15.9, 19.8, 24.2, 20.4, 27.7, 30.1, 21.7, 21.6, 21, 25.1, 27.1, 30.5, 39.8, 22.9, 21, 63.2, 26.8, 26.2, 29, 17.6, 20.8, 21.8, 22.8, 23.1, 17.3, 31.1, 27.6, 20.9, 23.1, 21.6, 24, 25.8, 15.6, 17.1, 15.8, 17.7, 23, 23, 19.2, 23.5, 24.8, 24, 24.2, 26.3, 26.2)

#real sample timees
fossil.ages <- c(15, 18, 20, 24, 30, 32, 33, 34, 40, 43, 45, 47, 49, 50, 51, 56)
fossil.ages <- 93 - sort(fossil.ages, decreasing=TRUE)
numSamplePoints = length(fossil.ages)
#number of fossils at each time point
numfossils=c(2, 2, 2, 13, 1 , 17, 4, 2, 20, 7, 13, 5, 4, 3, 2, 12) 
numfossils  <- numfossils[length(numfossils):1]

fast.tree = function(n, lambda=1, mu=1, sigma = 1, theta = 0)
{
  # start off doing things in reverse 
  # (start with n extant species, create species at rate mu, destroy at rate lambda)
  # events = -1 with prob l/(l+mu), 1 with prob mu/(l+mu)
  
  # calculate mean of distribution to predict number of likely events before
  # all n nodes are removed, then run for double that time to make sure we 
  # have enough samples to get there (TODO: make this more robust)
  
  m = (mu - lambda)/(mu + lambda)
  # track events
  changes = sample(c(-1, 1), size = -2*n/m, prob = c(lambda, mu), replace = TRUE)
  # track number of leaves initially (n), then after each event
  leaves = n + c(0, cumsum(changes))
  # find out where we hit zero leaves and crop just before that
  n.events = min(which(leaves == 0)) - 1
  # now turn things the right way round, excluding the first change (always +1)
  changes = -changes[(n.events-1):1]
  leaves = leaves[n.events:1]
  m=-m
  
  # generate times till next event based on current number of leaves
  times = rexp(n = n.events, rate = leaves * (mu + lambda))
  # for each event, randomly allocate a leaf for the event to happen to
  changer = ceiling(runif(n.events) * leaves)
  # allocate edge IDs using protocol for phylo class: first IDs go to tips
  # (i.e. 1:(n+L) below) then remaining IDs to internal nodes
  # note that the first n IDs go to extant tips
  nodes = 0*leaves
  ii = which(changes == -1)
  L = length(ii)
  nodes[ii] = n + (L:1)
  ii = which(changes == 1)
  nodes[ii] = n + L + 1:length(ii)
  # number of edges
  n.edge = n + n.events - 2
  
  curr.time = cumsum(times[-1])
  tree.time = sum(times[-1]) # time from first split to current day
  sampleTimes = fossil.ages * tree.time / 93
  n.samples = length(sampleTimes)
  sample.events = sapply(sampleTimes, function(x) which.min(curr.time < x))
  n.leaves = leaves[sample.events + 1]
  max.leaves = max(n.leaves)
  t.el = sampleTimes - curr.time[sample.events - 1]
  
  # use precompiled Fortran code to convert tree data to
  # phylo form, using a recursive function to calculate edges between nodes
  # also carries along trait information along edges, as well as at sampled points
  
  #  dyn.load('tree_climb.dll')
  
  #recursive subroutine tree_climb(n, n_events, leaves, changes, changer, nodes, times, time, a,&
  #edge, edge_length, edge_trait, n_samples, se, n_leaves, ml, t_el, samples, trait, sigma, ws)
  
  test = .Fortran('tree_climb', 
                  n = as.integer(n), 
                  n_events = as.integer(n.events), 
                  leaves = as.integer(leaves), 
                  changes = as.integer(changes), 
                  changer = as.integer(changer), 
                  nodes = as.integer(nodes), 
                  times = as.double(times), 
                  time = 1L, 
                  a = 1L, 
                  edge = matrix(0L, n + n.events - 2, 2), 
                  edge_length = numeric(n + n.events - 2),
                  edge_trait = numeric(n + n.events - 2),
                  n_samples = as.integer(n.samples),
                  se = as.integer(sample.events),
                  n_leaves = as.integer(n.leaves),
                  ml = as.integer(max.leaves),
                  t_el = as.double(t.el),
                  samples = matrix(as.double(0.0), max.leaves, n.samples),
                  trait = 0.0, #numeric(n + n_events - 2) # change value to change trait start
                  sigma = as.double(sigma),
                  theta = as.double(theta),
                  ws = 1L
  )
  res = list(edge = test$edge, 
             edge.length = test$edge_length, 
             edge.trait = test$edge_trait,
             tip.label = paste('t', 1:((n + n.events)/2), sep=''),
             root.edge = times[1],
             Nnode = ((n + n.events)/2) - 1
  )
  class(res) = 'phylo'
  list(tree = res, samples = test$samples, n.leaves = n.leaves, sampleTimes = sampleTimes)
}

#dyn.load('tree_climb.dll')

for (theta in theta_values) {
  for (bdratio in bdratio_values){
    print(sprintf("Theta: %f, BD Ratio: %f",theta,bdratio))
    print(paste('BD RATIO', bdratio))
    rand.correlations = numeric(REPS)
    rand.sds = numeric(REPS)
    rand.sigma2 = numeric(REPS)
    
    for (reps in 1: REPS){
      if (reps%%100 == 0) print(reps)
      res = fast.tree(n = numtips, lambda = bdratio, mu = 1, sigma = sigma)
      # output names: tree, samples, n.leaves, sampleTimes
      
      node.trait = numeric(max(res$tree$edge[,2]))
      node.trait[res$tree$edge[,2]] = res$tree$edge.trait
      
      #set up a list of lists to store the samples in
      df <- NULL #this data frame will store samples ages and trait values 
      
      for (i in 1:numSamplePoints) {
        theseSamples <- res$samples[1:res$n.leaves[i], i] #all the species at this time slice
        theseFossils <- sample(theseSamples, numfossils[i], replace=TRUE) #replacement or not?
        df <- rbind(df, cbind(theseFossils, res$sampleTimes[i]))
      }
      
      #data frame containing the correlation between random values and zachos
      #each row corresponds to an individual sample
      rand.correlations[reps] = cor(df[ ,1], zachos)
      rand.sds[reps] = sd(df[,1])
      # normalised sigma = tree length / sd@tips
      rand.sigma2[reps] = max(node.depth.edgelength(res$tree))/sd(node.trait[1:numtips])^2
    }
    
    results <- data.frame(rand.corr = rand.correlations, rand.sd = rand.sds, rand.sigma2 = rand.sigma2)
    filename = paste("results", bdratio, theta, REPS, ".txt", sep="_")
    write.table(results, file = filename, quote=FALSE, row.names=FALSE)
  }
}
#dyn.unload('tree_climb.dll')

## plot the last tree coloured by traits
#palette(rainbow(1024)[1:700])
#col = res$tree$edge.trait
#col = (col - min(col))/(max(col) - min(col))*700 + 1
#plot(res$tree, edge.color = col, show.tip.label = FALSE)

## test methods outputs against each other to make sure they look sensible
#plot(c(rand.correlations, rand.correlations.old), col = rep(c('black', 'red'), times = c(1000, 100)))
#plot(c(rand.sds, rand.sds.old), col = rep(c('black', 'red'), times = c(1000, 100)))


teststat <- cor(zachos, gcl)

#process results of random correlation sims
for (theta in theta_values) {
  for (bdratio in bdratio_values){
    filename = paste("results", bdratio, theta, REPS, ".txt", sep="_")
    simres <- read.table(file=filename, header=TRUE)
    print(sprintf('Theta %f, bdratio %f',theta, bdratio))
    print(paste('p = ', sum(abs(simres$rand.corr) > abs(teststat), na.rm = TRUE), '/', sum(!is.na(simres$rand.corr))))
    print(paste('overlap = ', sum(dnorm(simres$rand.sigma2, 1.181, 0.1983), na.rm = TRUE) ))
    print('')
  #  print(summary(simres))
    tmp = hist(simres$rand.sigma2, 100, main = paste('bdratio = ', bdratio)); abline(v = 1.181)
    lines(seq(0.5,3,0.01), dnorm(seq(0.5,3,0.01), 1.181, 0.1983)*max(tmp$counts)/dnorm(0,0,0.1983), col='red')
    Sys.sleep(2)
  }
}
