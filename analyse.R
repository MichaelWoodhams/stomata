setwd("/home/mdw2/Documents/botany/stomata/code")
library(ape)
library(TreeSim)
real.sigma = c(2.4308543, 1.18148, 1.14046, 0.69834, 1.71374)
real.se = c(0.4080669, 0.19833, 0.19145, 0.11724, 0.28769)
REPS = 10000 # Just used for finding filenames.
root.value = 0 #starting value for the trait arbitrarily set to zero
numtips = 1700 #this is set to match what is known about Proteaceae
#might want to play with next two
bdratio_values = c(1.1, 1.3, 1.5, 2, 3, 5, 10, 20, 50, 100)
#theta_values = c(0,0.1,0.5,1)
theta_values = c(0) # get it working like this first.
#climate_values = 1:6
climate_values = 6 # get it working like this first.

sigma = 1


poo = array(NA, c(3, 10, 5, 6))

#dataset = 5
#climate = 6

for (dataset in 1:5){
  for (climate in climate_values){
    
    corr = sprintf("rand.corr%d", dataset)
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
    
    data.vec = data[keep, 10 + dataset]
    ij = which(!is.na(data.vec))
    print(paste('dataset', dataset, 'climate', climate, 'samples', length(ij), 'ages', min(ages[ij]), max(ages[ij])))
    
    
    teststat <- cor(data.vec[ij], clim.vec[ij]) # cor(zachos, gcl)

    #process results of random correlation sims
    p.values = numeric(length(bdratio_values))
    overlap = numeric(length(bdratio_values))
    overlap2 = numeric(length(bdratio_values))
    for (theta in theta_values) {
      i = 0
      for (bdratio in bdratio_values) { 
        i = i + 1
        filename = paste("../results/results", theta, bdratio, REPS, "climate", climate, ".txt", sep="_")
        simres <- read.table(file=filename, header=TRUE)
        
        print(paste('bdratio',bdratio))
        p.values[i] = sum(abs(simres[[corr]]) > abs(teststat), na.rm = TRUE)/
          sum(!is.na(simres[[corr]]))
        print(paste('p = ', p.values[i]))
        overlap[i] = sum(
          #dnorm(simres$rand.sigma2, real.sigma[dataset], real.se[dataset]), 
          pnorm(simres$rand.sigma2, real.sigma[dataset], real.se[dataset]),     
          na.rm = TRUE)/REPS
        print(paste('overlap = ', overlap[i] ))
        sums2 = simres$rand.ace.se^2 + real.se[dataset]^2
        overlap2[i] = sum(
          #    exp( -(simres$rand.ace.sigma2 - real.sigma[dataset])^2/2/sums2) / 
          #    sqrt(2*pi*sums2)
          #dnorm(real.sigma[dataset], simres$rand.ace.sigma2, sqrt(sums2)), 
          pnorm(real.sigma[dataset],  simres$rand.ace.sigma2, sqrt(sums2)), 
          na.rm = TRUE) / REPS
        print(paste('overlap2 = ', overlap2[i] ))
        print('')
        #  print(summary(simres))
        #tmp = hist(simres$rand.sigma2, 100, main = paste('bdratio = ', bdratio)); abline(v = 1.181)
        #  lines(seq(0.5,3,0.01), dnorm(seq(0.5,3,0.01), 1.181, 0.1983)*max(tmp$counts)/dnorm(0,0,0.1983), col='red')
        #  Sys.sleep(2)
      }
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
for (k in 1:length(theta_values)) {
  for (i in 1:length(bdratio_values))
  {
    bdratio = bdratio_values[i]
    theta = theta_values[k]
    bdratio.next = ifelse(bdratio == 100, 200, bdratio_values[i+1])
    filename = paste("../results/results", theta, bdratio, REPS, "climate", climate, ".txt", sep="_")
    simres <- read.table(file=filename, header=TRUE)
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
for (k in 1:length(theta_values)) {
  for (i in 1:length(bdratio_values))
  {
    bdratio = bdratio_values[i]
    theta = theta_values[k]
    filename = paste("../results/results", theta, bdratio, REPS, "climate", climate, ".txt", sep="_")
    simres <- read.table(file=filename, header=TRUE)
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
}


L = lm(simres$rand.ace.sigma2 ~ simres$rand.sigma2)
conv = simres$rand.sigma2*L$coefficients[2] + L$coefficients[1]
t = seq(0, max(c(conv, simres$rand.ace.sigma2))+0.05, 0.05)
hist(simres$rand.ace.sigma2, t, ylim = c(0,4000), main = bdratio_values[i])
hist(conv, t, col = '#FF000080', add = TRUE)

# 
# 
