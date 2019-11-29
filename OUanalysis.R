setwd("results")

measure.stats <- read.csv("../OU_parameters_full_data.csv")
measurement.names <- substr(measure.stats[seq(2,10,2),1],7,100) # "GCL"   "SI"    "SD"    "EA"    "Gsmax"
measurement.sigma <- sqrt(measure.stats[seq(1,9,2),"mean"])
measurement.alpha <-      measure.stats[seq(2,10,2),"mean"]

files <- list.files(pattern="compare.*csv")
all.data <- NULL
for (file in files) {
  data <- read.csv(file)
  measurement <- strsplit(file,"_")[[1]][2]
  index <- which(measurement==measurement.names)
  data$measurement <- measurement
  data$true.sigma <- measurement.sigma[index]
  data$true.alpha <- measurement.alpha[index]
  bdratio <- strsplit(file,"_")[[1]][3] # looks like "1.100000.csv" currently
  data$bdratio <- as.numeric(substring(bdratio,1,nchar(bdratio)-4)) # looks like 1.1
  print(sprintf("%d %s %s %s %s",index,measurement,measurement.sigma[index],measurement.alpha[index],data$bdratio[1]))
  if (is.null(all.data)) {
    all.data <- data
  } else {
    all.data <- rbind(all.data,data)
  }
}

# Need to name this file so it doesn't get caught by list.files() above if we rerun.
write.csv(all.data,file="CompareAll.csv")
# Check there are no NAs
all.data[is.na(all.data$orderFastAlpha),]
good.data <- all.data[!is.na(all.data$orderFastAlpha),]  # No NAs, so not needed
err.sigma  <- (good.data$true.sigma-good.data$sigma )/good.data$true.sigma
err.alpha  <- (good.data$true.alpha-good.data$alpha )/good.data$true.alpha
err.sigma1 <- (good.data$true.sigma-good.data$sigma1)/good.data$true.sigma
err.alpha1 <- (good.data$true.alpha-good.data$alpha1)/good.data$true.alpha

sderr <- sd(err.sigma )/sqrt(length(err.sigma ))
print(sprintf("%s mean=%f, sd=%f, sderr=%f, mean/sderr=%f","sigma", mean(err.sigma ),sd(err.sigma ),sderr,mean(err.sigma )/sderr))
sderr <- sd(err.sigma1)/sqrt(length(err.sigma1))
print(sprintf("%s mean=%f, sd=%f, sderr=%f, mean/sderr=%f","sigma1",mean(err.sigma1),sd(err.sigma1),sderr,mean(err.sigma1)/sderr))
sderr <- sd(err.alpha )/sqrt(length(err.alpha ))
print(sprintf("%s mean=%f, sd=%f, sderr=%f, mean/sderr=%f","alpha", mean(err.alpha ),sd(err.alpha ),sderr,mean(err.alpha )/sderr))
sderr <- sd(err.alpha1)/sqrt(length(err.alpha1))
print(sprintf("%s mean=%f, sd=%f, sderr=%f, mean/sderr=%f","alpha1",mean(err.alpha1),sd(err.alpha1),sderr,mean(err.alpha1)/sderr))


