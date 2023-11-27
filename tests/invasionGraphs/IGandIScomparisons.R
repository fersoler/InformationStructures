#################################################################
# Perfornance comparison of IS and IG creation
# Author: Fernando Soler-Toscano - fsoler@us.es
#################################################################

library(tidyverse)
source("tests/invasionGraphs/invasion_graph_main_functions.R")
source("R/ISgraph.R")
source("R/ISbuild.R")

kVals <- c(5,7,10,12,14)
nVals <- c(0,0.1,0.3,0.5,0.7)
nMat <- 10

results <- matrix(0,nrow = length(kVals)*length(nVals)*nMat,ncol=5)
colnames(results) <- c("kVal", "nVal", "IGtime", "IGQuicktime", "IStime")
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
      # Time of InvasionScheme
      start_time <- Sys.time()
      InvSch <- LV.IS(A,b)
      end_time <- Sys.time()
      timeInvSch <- difftime(end_time, start_time, units="secs")[[1]]
      # Time of IG
      start_time <- Sys.time()
      IG.function(InvSch)
      end_time <- Sys.time()
      timeIG <- timeInvSch+difftime(end_time, start_time, units="secs")[[1]]
      # Time of IGquick
      start_time <- Sys.time()
      IG.quick.function(InvSch)
      end_time <- Sys.time()
      timeIGquick <- timeInvSch+difftime(end_time, start_time, units="secs")[[1]]
      # Time of IS
      start_time <- Sys.time()
      ISbuild(as.data.frame(t(b)),A)
      end_time <- Sys.time()
      timeIS <- difftime(end_time, start_time, units="secs")[[1]]
      # Record data
      data <- c(k, n, timeIG, timeIGquick, timeIS)
      results[nextRec,] <- data
      print(paste("         Recorded data: ",toString(data),sep = ""))
      # Print and prepare next one
      nextRec <- nextRec+1
    }
  }
}


#write.csv(results,"tests/invasionGraphs/IGandIScomparison.csv", row.names = FALSE)

results <- read.csv("tests/invasionGraphs/IGandIScomparison.csv")

tableTimes <- data.frame(species=integer(), n=double(),seconds=double(),
                         algorithm=character(),stringsAsFactors = FALSE)
for(i in 1:nrow(results)){
  tableTimes[1+3*(i-1),] <- c(results[i,c("kVal", "nVal","IGtime")],"IG")
  tableTimes[2+3*(i-1),] <- c(results[i,c("kVal", "nVal","IGQuicktime")],"IGq")
  tableTimes[3+3*(i-1),] <- c(results[i,c("kVal", "nVal","IStime")],"IS")
}

n_names <- c(
  `0`   = "b in [0,1]",
  `0.1` = "b in [-0.1,0.9]",
  `0.3` = "b in [-0.3,0.7]",
  `0.5` = "b in [-0.5,0.5]",
  `0.7` = "b in [-0.7,0.3]")

pdf(file="tests/invasionGraphs/pointsIGandIS.pdf", width = 10, height = 5)
ggplot(data = filter(tableTimes, species>1)) +
  geom_point(mapping = aes(x = species, y = seconds,color = algorithm), position = "jitter") + 
  scale_size(trans="log10") +
  facet_wrap(~ n, scales = "free_y", nrow = 2, labeller=as_labeller(n_names)) + 
  scale_y_continuous(trans=scales::log_trans(), breaks=c(1,5,10,100,1000,5000)) +
  ggtitle("Creation of Invasion Graphs (general 'IG' and quick 'IGq' functions) and Info. Structures (IS)")
dev.off()

pdf(file="tests/invasionGraphs/boxesIGandIS.pdf", width = 10, height = 5)
ggplot(filter(tableTimes,species>1),
       aes(x = species, y = seconds, color = algorithm)) +
  geom_boxplot(aes(factor(species), seconds)) +
  scale_size(trans="log10") +
  facet_wrap(~ n, scales = "free_y", nrow = 2, labeller=as_labeller(n_names)) + 
  scale_y_continuous(trans=scales::log_trans(), breaks=c(1,5,10,100,1000,5000))  +
  ggtitle("Creation of Invasion Graphs (general 'IG' and quick 'IGq' functions) and Info. Structures (IS)")
dev.off()
