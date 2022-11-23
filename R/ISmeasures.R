source("R/ISbuild.R")
source("R/ISgraph.R")

# IS measures from:
# https://doi.org/10.1371/journal.pcbi.1010412

getISmeasures <- function(IS,ISgr,alphas,gMatrix){
  noelAndPaths <- noel_and_npaths(IS, ISgr)
  coopVals <- coopMeasures(IS, alphas)
  list(
    noel   = noelAndPaths$noel,
    frond  = frondosity(IS),
    crit   = criticality(IS, alphas, gMatrix),
    sync   = synchronicity(IS),
    npaths = noelAndPaths$npaths,
    coopLH = coopVals$highestLev,
    coopA  = coopVals$coopA,
    coopB  = coopVals$coopB,
    coopC  = coopVals$coopC
  )
}

# Criticality value of the GASS
criticality <- function(IS,alphas,gMatrix){
  size<-dim(alphas)[2]            # Size of the system
  m <- -gMatrix
  globalSol <- IS$points[IS$gassInd,]
  v<- (m %*% globalSol)-alphas + globalSol # Sum of w and z vectors
  min(abs(v/sqrt(sum(v^2))))          # Criticality value
}

# Synchronicity
# Input:
# - IS: Information Structure
# Output:
# 1 if all nodes are synchronized and 0 otherwise
synchronicity <- function(IS){
  gass <- IS$points[IS$gassInd,]
  sync <- 0
  if(min(gass) > 0 || max(gass)==0){
    sync <- 1
  }
  sync
}

# Get the frondosity of an IS given its stable points
frondosity <- function(IS){
  rs <- colSums(IS$points)
  (dim(IS$points)[1])/(2^length(rs[rs>0]))
}

# Given an IS and its graph returns the number of
# energy levels (NoEL).
noel <- function(IS, ISgr){
  if(IS$gassInd==1)
    1
  else
    max(sapply(all_simple_paths(ISgr$graph, from=1, to=IS$gassInd),length))
}

# NoEL and number of paths from 0 to the GASS
noel_and_npaths <- function(IS, ISgr){
  if(IS$gassInd==1){
    noel = 1
    npaths = 0
  }
  else {
    allPaths = all_simple_paths(ISgr$graph, from=1, to=IS$gassInd)
    noel = max(sapply(allPaths,length))
    npaths = length(allPaths)
  }
  list(
    noel = noel,
    npaths = npaths
  )
}

# Cooperative points
# Input:
# - IS: output of ISbuild
# Output:
# - Matrix of cooperative points (may be empty)
coopPoints <- function(IS){
  points <- IS$points
  conn <- IS$connectivity
  nP <- nrow(points)
  size <- ncol(points)
  coopPoints <- as.data.frame(matrix(0,0,size))
  for(I in 2:nP){
    goingToI <- (1:nP)[conn[,I]==1]
    if(sum(goingToI)>1){
      prev <- colSums(points[goingToI,,drop=FALSE])>0
      newVals <- (1:size)[prev == 0 & points[I,] > 0]
      if(length(newVals)>0){
        coopPoints[nrow(coopPoints)+1,] <- points[I,] 
      }
    }
  }
  coopPoints
}

# Get cooperative measures (for cooperative IS)
# Input:
# - IS: Information structure (output of ISbuild)
# - alphas: vector of alpha values
# Output:
# - highestLev: highest cooperative level
# - coopA: cooperation value A
# - coopB: cooperation value B
# - coopC: cooperation value C
coopMeasures <- function(IS, alphas){
  cooperativePoints <- coopPoints(IS)
  allLevels <- apply(cooperativePoints, 1, function(x) nodeLevel(x))
  sumAll <- colSums(IS$points)
  ca <- sum(allLevels) - length(allLevels)
  cb <- 0
  cc <- 0
  for(i in 1:length(alphas)){
    if(alphas[i]<=0 && sumAll[i]>0){
      found <- FALSE
      j <- 1
      while(!found && j<=dim(cooperativePoints)[1]){
        if(cooperativePoints[j,i]>0){
          found <- TRUE
        }
        j <- j+1
      }
      ln <- nodeLevel(cooperativePoints[j-1,])-1
      cb <- cb+ln
      cc <- cc+(2^ln)
    }
  }
  list(highestLev = max(c(0,allLevels)),
       coopA = ca,
       coopB = cb,
       coopC = cc)
}

# Auxiliary function to test whether sol1 is a subset of sol2 by
# considering non-zero values
isSubset <- function(sol1, sol2){
  all(sol2[(1:length(sol1))[sol1 > 0]] > 0)
}
# Auxiliary function to get the leven (number of non-zero values) of
# a point in a cooperative IS
nodeLevel <- function(sol){
  sum(rep(1,length(sol))[sol > 0])
}

