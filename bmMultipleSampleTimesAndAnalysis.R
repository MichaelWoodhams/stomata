setwd("/home/mdw2/Documents/botany/stomata/code")
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

proc = 1
procs = 2

REPS = 5000
root.value = 0 #starting value for the trait arbitrarily set to zero
numtips = 1700 #this is set to match what is known about Proteaceae
#might want to play with next two
#bdratio_values = c(1.1, 1.5, 2, 5, 10, 20, 50, 100)
bdratio_values = c(1.1, 1.3, 1.5, 2, 3, 5, 10, 20, 50, 100)
theta_values = c(0,0.1,0.5,1)
sigma = 1
debug=TRUE
if(debug) {
  # Set values for small runs
  REPS=10  
}


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

for (theta in theta_values){
  for (bdratio in bdratio_values){
    print(paste('Theta:',theta, '   BD RATIO', bdratio))
  
    rand.correlations = matrix(0, REPS, 5)
    rand.sds = matrix(0, REPS, 5)
    rand.sigma2 = numeric(REPS)
    rand.ace.sigma2 = numeric(REPS)
    rand.ace.se = numeric(REPS)
  
    for (reps in 1:REPS){
      if (reps%%100 == 0) print(reps)
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
    filename = paste("../results/results", theta, bdratio ,REPS, "proc", proc, "climate", climate, ".txt", sep="_")
    write.table(results, file = filename, quote=FALSE, row.names=FALSE)
  }
}


## plot the last tree coloured by traits
palette(rainbow(1024)[1:700])
col = res$tree$edge.trait
col = (col - min(col, na.rm = TRUE))/(max(col, na.rm = TRUE) - min(col, na.rm = TRUE))*700 + 1
l = max(node.depth.edgelength(res$tree))

par(lwd = 2, cex = 1.5)
plot(res$tree, edge.color = col, show.tip.label = FALSE)
#axis(1)
v = c(93, 75, 60, 45, 30, 15, 0)
axis(1, at = (93-v)*l/93, labels = v)
mtext('Millions of years ago', 1, line = 2.5, cex=1.5)


#df[,2] = 93 - (df[,2]*93/l)
for (t in unique(df[,2]))
{
  tr = df[df[,2]==t, 1]
  col = (tr - min(res$tree$edge.trait, na.rm = TRUE))/(max(res$tree$edge.trait, na.rm = TRUE) - min(res$tree$edge.trait, na.rm = TRUE))*700 + 1
  points(rep(t, length(tr)), seq(0, res$tree$Nnode, length.out = length(tr)+2)[-c(1, (length(tr)+2))], bg = col, 
         cex = 1, pch=21, col = 'black')
}

tr= res$tree$edge.trait[sample(1:1700, 71)]
t = l
col = (tr - min(res$tree$edge.trait, na.rm = TRUE))/(max(res$tree$edge.trait, na.rm = TRUE) - min(res$tree$edge.trait, na.rm = TRUE))*700 + 1
points(rep(t, length(tr)), seq(0, res$tree$Nnode, length.out = length(tr)+2)[-c(1, (length(tr)+2))], bg = col, 
       cex = 1, pch=21, col = 'black')


bush$root.edge = max(node.depth.edgelength(res$tree)) - max(node.depth.edgelength(bush))
col = bush$edge.trait
col = (col - min(col, na.rm = TRUE))/(max(col, na.rm = TRUE) - min(col, na.rm = TRUE))*700 + 1
plot(bush, edge.color = col, show.tip.label = FALSE, root.edge = TRUE); axis(1)
nodelabels(text = round(node.trait[71 + 1:bush$Nnode], 3), cex =0.2)
tiplabels(text = round(node.trait[1:71], 3), cex =0.2)

## test methods outputs against each other to make sure they look sensible
#plot(c(rand.correlations, rand.correlations.old), col = rep(c('black', 'red'), times = c(1000, 100)))
#plot(c(rand.sds, rand.sds.old), col = rep(c('black', 'red'), times = c(1000, 100)))


stop('Random sims done')

# Process output

setwd("/home/mdw2/Documents/botany/stomata/code")
library(ape)
library(TreeSim)
real.sigma = c(2.4308543, 1.18148, 1.14046, 0.69834, 1.71374)
real.se = c(0.4080669, 0.19833, 0.19145, 0.11724, 0.28769)
REPS = 5000
root.value = 0 #starting value for the trait arbitrarily set to zero
numtips = 1700 #this is set to match what is known about Proteaceae
#might want to play with next two
bdratio_values = c(1.1, 1.3, 1.5, 2, 3, 5, 10, 20, 50, 100)
theta_values = c(0,0.1,0.5,1)
sigma = 1


poo = array(NA, c(3, 10, 5, 6))

#dataset = 5
#climate = 6

for (dataset in 1:5){
  for (climate in 1:6){
    
    corr = sprintf("rand.corr%d", dataset)
    data = read.csv('zachos correl and calc gsmax1.csv')
    clim.vec = data[,c(3:6, 5:6)[climate]]
    if (climate <= 4)
    {
      keep = which(!is.na(clim.vec) & data$delete == 0)
    } else {
      keep = which(!is.na(clim.vec))
    }
    
    clim.vec = clim.vec[keep]
    
    data.vec = data[keep, 10 + dataset]
    ij = which(!is.na(data.vec))
    print(paste('dataset', dataset, 'climate', climate, 'samples', length(ij), 'ages', min(ages[ij]), max(ages[ij])))
    
    
    teststat <- cor(data.vec[ij], clim.vec[ij]) # cor(zachos, gcl)
    procs = 2
    
    #process results of random correlation sims
    p.values = numeric(length(bdratio_values))
    overlap = numeric(length(bdratio_values))
    overlap2 = numeric(length(bdratio_values))
    i = 0
    for (bdratio in bdratio_values){
      i = i + 1
      simres = NULL
      for (j in 1:procs)
      {
        filename = paste("../results/results", bdratio, REPS, "proc", j, "climate", climate, ".txt", sep="_")
        simres <- rbind(simres, read.table(file=filename, header=TRUE))
      }
      print(paste('bdratio',bdratio))
      p.values[i] = sum(abs(simres[[corr]]) > abs(teststat), na.rm = TRUE)/
        sum(!is.na(simres[[corr]]))
      print(paste('p = ', p.values[i]))
      overlap[i] = sum(
        #dnorm(simres$rand.sigma2, real.sigma[dataset], real.se[dataset]), 
        pnorm(simres$rand.sigma2, real.sigma[dataset], real.se[dataset]),     
        na.rm = TRUE)/(REPS*procs)
      print(paste('overlap = ', overlap[i] ))
      sums2 = simres$rand.ace.se^2 + real.se[dataset]^2
      overlap2[i] = sum(
        #    exp( -(simres$rand.ace.sigma2 - real.sigma[dataset])^2/2/sums2) / 
        #    sqrt(2*pi*sums2)
        #dnorm(real.sigma[dataset], simres$rand.ace.sigma2, sqrt(sums2)), 
        pnorm(real.sigma[dataset],  simres$rand.ace.sigma2, sqrt(sums2)), 
        na.rm = TRUE) / (REPS*procs)
      print(paste('overlap2 = ', overlap2[i] ))
      print('')
      #  print(summary(simres))
      #tmp = hist(simres$rand.sigma2, 100, main = paste('bdratio = ', bdratio)); abline(v = 1.181)
      #  lines(seq(0.5,3,0.01), dnorm(seq(0.5,3,0.01), 1.181, 0.1983)*max(tmp$counts)/dnorm(0,0,0.1983), col='red')
      #  Sys.sleep(2)
    }
    
    overlap = 1 - 2*abs(overlap - 0.5)
    overlap2 = 1 - 2*abs(overlap2 - 0.5)
    
    poo[1, , dataset, climate] = p.values
    poo[2, , dataset, climate] = overlap
    poo[3, , dataset, climate] = overlap2
    
  }
}





par(cex=1.5, lwd=2, mgp = c(3.5, 1, 0), mar = c(5,5,4,2)+0.1); 
plot(c(-3, 6), c(0,4), type='n', log='', xlab = expression(frac(lambda, mu)), 
     ylab = expression(sigma^{2}), axes = FALSE)
axis(1, at = log(bdratio_values - 1), labels = bdratio_values)
axis(2)
box(bty = 'O')
abline(v = log(bdratio_values - 1), col='grey', lty=2)
ds = 0.05
max.sigma = 10
breaks = seq(0, max.sigma, ds)
mids = seq(ds/2, max.sigma - ds/2, ds)
t = seq(0, 5.5, 0.01)
for (i in 1:length(bdratio_values))
{
  bdratio = bdratio_values[i]
  bdratio.next = ifelse(bdratio == 100, 200, bdratio_values[i+1])
  simres = NULL
  for (j in 1:procs)
  {
    filename = paste("../results/results", bdratio, REPS, "proc", j, "climate", climate, ".txt", sep="_")
    simres <- rbind(simres, read.table(file=filename, header=TRUE))
  }
  tmp = hist(simres$rand.sigma2, breaks, plot = FALSE)$counts
  tmp = tmp/max(tmp)*0.8
  for (j in 1:length(mids)) 
    lines(log(bdratio - 1) + 
            c(0, tmp[j]*(log(bdratio.next - 1) - log(bdratio - 1))), 
          rep(mids[j], 2), col = '#000000B0', lwd=4)

  tmp = hist(simres$rand.ace.sigma2, breaks, plot = FALSE)$counts
  tmp = tmp/max(tmp)*0.8
  for (j in 1:length(mids)) 
    lines(log(bdratio - 1) + 
            c(0, tmp[j]*(log(bdratio.next - 1) - log(bdratio - 1))), 
          rep(mids[j], 2), col = '#0000FFB0', lwd=4)
  
  for (dataset in 1:5)
  {
    tmp = dnorm(t, real.sigma[dataset], real.se[dataset])
    tmp = tmp/max(tmp) *0.8
    
    lines(log(bdratio - 1) + tmp*(log(bdratio.next - 1) - log(bdratio - 1)), 
          t, col=c('red', 'orange', 'yellow', 'green', 'purple')[dataset], lwd = 3)
  }
}


# pvt = paste('=', format(p.values, sci = FALSE))
# pvt[pvt=="= 0.0000"] = "< 0.0001"
# text(log(bdratio_values - 1), rep(4, 10), 
#      paste('p', pvt, 
#            '\noverlap =', round(overlap,3), 
#            '\noverlap2 =', round(overlap2,3)), cex=0.5, srt = 90)



# 
par(cex = 1.5, lwd = 2)
palette(rainbow(16))
for (i in 1:length(bdratio_values))
{
  bdratio = bdratio_values[i]
  simres = NULL
  for (j in 1:procs)
  {
    filename = paste("../results/results", bdratio, REPS, "proc", proc, "climate", climate, ".txt", sep="_")
    simres <- rbind(simres, read.table(file=filename, header=TRUE))
  }
  if (i == 1)
    plot(simres$rand.sigma2, simres$rand.ace.sigma2,
         xlim = c(0,5), ylim = c(0,5), col = i, pch = 19,
         xlab = 'Full tree sigma^2', ylab = 'Pruned tree estimated sigma^2')
  else
    points(simres$rand.sigma2, simres$rand.ace.sigma2, col = i, pch = 19)

  print(sum(simres$rand.ace.sigma2 > simres$rand.sigma2))
}
abline(c(0,1), col='red', lty=2)
legend('bottomright', title = 'lambda/mu', legend = bdratio_values, col = 1:length(bdratio_values), pch=19, cex = 0.5)


L = lm(simres$rand.ace.sigma2 ~ simres$rand.sigma2)
conv = simres$rand.sigma2*L$coefficients[2] + L$coefficients[1]
t = seq(0, max(c(conv, simres$rand.ace.sigma2))+0.05, 0.05)
hist(simres$rand.ace.sigma2, t, ylim = c(0,4000), main = bdratio_values[i])
hist(conv, t, col = '#FF000080', add = TRUE)

# 
# 
