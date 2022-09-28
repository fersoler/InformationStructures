source("R/ISbuild.R")
source("R/ISgraph.R")
source("R/ISmeasures.R")

library(tidyverse)

############################################################
# Random alpha values
############################################################

getRandomAlphas <- function(n){
  alp <- rnorm(n)
  alp/sqrt(sum(alp^2))
}

nodes   <- 2    # Size of each alpha vector
alpVals <- 100   # How many random alpha vectors

#set.seed(20220923)
randAlphas <- matrix(0,nrow = alpVals, ncol = nodes)
for(r in 1:alpVals){
  randAlphas[r,] <- getRandomAlphas(nodes)
}


plot(c(-1, 1), c(-1, 1), type = "n", axes = FALSE, ann = FALSE, asp = 1)
points(randAlphas,pch=16,cex = 0.3)

