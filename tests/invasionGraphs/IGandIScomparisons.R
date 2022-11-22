library(tidyverse)
source("tests/invasionGraphs/invasion_graph_main_functions.R")
source("R/ISgraph.R")
source("R/ISbuild.R")
source("R/ISbuildImproved.R")


kVals <- c(5,7,10,12,14)
nVals <- c(0,0.1,0.3,0.5,0.7)
nMat <- 10

results <- matrix(0,nrow = length(kVals)*length(nVals)*nMat,ncol=8)
colnames(results) <- c("kVal", "nVal", "IGtime", "IStime", "QP","LCP", 
                       "subcTst","subcChk")
nextRec <- 1
for(k in kVals){
  print(paste("Checking for k = ",k," species",sep=""))
  for(n in nVals){
    print(paste("   * n = ",n,sep=""))
    for(i in 1:nMat){
      print(paste("       matrix number = ",i,sep=""))
      # Generate A and b
      m <- 1/k
      A <- -diag(k)+m*matrix(runif(k^2)-0.5,k,k)
      diag(A) <- -1
      b <- matrix(runif(k)-n,k,1)
      # Time of IG
      start_time <- Sys.time()
      # Invasion graph
      IG.function(LV.IS(A,b))
      end_time <- Sys.time()
      timeIG <- difftime(end_time, start_time, units="secs")[[1]]
      # Time of IS
      start_time <- Sys.time()
      ISbuildThird(as.data.frame(t(b)),A)
      end_time <- Sys.time()
      timeIS <- difftime(end_time, start_time, units="secs")[[1]]
      # IS statistics
      IS <- ISbuildThirdWithCounters(as.data.frame(t(b)),A)
      nQP <- IS$nQPsolve
      nLCP <- IS$nLCPsolve
      nTry <- IS$nREUSEmed
      nChk <- IS$nREUSEchecked
      # Record data
      data <- c(k, n, timeIG, timeIS, nQP, nLCP, nTry, nChk)
      results[nextRec,] <- data
      print(paste("         Recorded data: ",toString(data),sep = ""))
      # Print and prepare next one
      nextRec <- nextRec+1
    }
  }
}


#results2 <- as.data.frame(results)
##write.csv(results2,"tests/invasionGraphs/IGandIScomparison.csv", row.names = FALSE)
results2 <- read.csv("tests/invasionGraphs/IGandIScomparison.csv")

#### 
tableTimes <- data.frame(species=integer(), n=double(),seconds=double(),algorithm=character(),stringsAsFactors = FALSE)
for(i in 1:nrow(results2)){
  tableTimes[1+2*(i-1),] <- c(results2[i,c("kVal", "nVal","IGtime")],"IG")
  tableTimes[2+2*(i-1),] <- c(results2[i,c("kVal", "nVal","IStime")],"IS")
}

n_names <- c(
  `0` = "alphas in [0,1]",
  `0.1` = "alphas in [-0.1,0.9]",
  `0.3` = "alphas in [-0.3,0.7]",
  `0.5` = "alphas in [-0.5,0.5]",
  `0.7` = "alphas in [-0.7,0.3]")


ggplot(data = filter(tableTimes,k>1)) +
  geom_point(mapping = aes(x = species, y = seconds, color = algorithm), position = "jitter") + 
  scale_size(trans="log10") +
  facet_wrap(~ n, scales = "free_y", nrow = 2, labeller=as_labeller(n_names)) + 
  scale_y_continuous(trans=scales::log_trans(), breaks=c(1,5,10,50,100,500,1000,5000))



ggplot(filter(tableTimes,k>1),
       aes(x = species, y = seconds, color = algorithm)) +
  geom_boxplot(aes(factor(species), seconds)) +
  scale_size(trans="log10") +
  facet_wrap(~ n, scales = "free_y", nrow = 2, labeller=as_labeller(n_names)) + 
  scale_y_continuous(trans=scales::log_trans(), breaks=c(1,5,10,50,100,500,1000,5000))
