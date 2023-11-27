# IG measures

library(vegan)

# Compute the node frondosity given the points
IGnodeFrond <- function(points){
  sp <- sum(colSums(points)>0)
  nrow(points)/(2^sp)
}

# Get the indexes of the GASS-like points
IGgassIndex <- function(out, points){
  IS <- out$IS
  npoints <- nrow(IS)
  nspec <- ncol(IS)
  # Number of possitives
  npos <- nspec+1  # -1 indicates none found
  # Number of zeros
  nzer <- 0
  # GASSes
  gasses <- NULL

  for(n in 1:npoints){
    if(out$permanent[[n]] & min(points[n,])>=0){
      g0 <- sum(IS[n,] > 0) # Number of values greater than 0
      nz <- sum(IS[n,] == 0) # Number of 0 values
      if(g0==npos){
        if(nz == nzer){
          gasses <- c(gasses,n)
        }
        if(nz>nzer){
          nzer <- nz
          gasses <- c(n)
        }
      }
      if(g0<npos){
        npos <- g0
        nzer <- nz
        gasses <- c(n)
      }
    }
  }
  if(is.null(gasses)){
    gasses <- c(1)
  }
  gasses
}

# Number of species in GASS-like points
IGspeciesGASS <- function(out,points){
  gasses <- IGgassIndex(out,points)
  nsp <- 0
  for(g in gasses){
    nsp <- nsp + out$number.species[g]
  }
  if(length(gasses)>0){
    nsp <- nsp/length(gasses)
  }
  nsp
}

# Abundance in the GASS, or mean abundance in case of
# several GASS-like points
IGmeanAbundGASS <- function(abund,gss){
  gassBund <- abund[gss,]
  if(length(gss)>1){
    gassBund <- colMeans(gassBund)
  }
  gassBund
}

# Evenness of the GASS, or mean evenness in case of
# several GASS-like points
IGmeanEvenness <- function(abund,gss){
  gassBund <- abund[gss,]
  eVals <- diversity(gassBund)/log(specnumber(gassBund))
  mean(eVals,na.rm = TRUE)
}

# This function provides, for each specie, the sum of the
# normalized vertex betweenness of all communities
# containing that specie.
IGspeciesBetweenness <- function(out){
  species <- ncol(out$IS)
  # Build the graph
  gr <- graph_from_adjacency_matrix(out$IG)
  bet <- betweenness(gr, normalized = TRUE)
  speciesBetw <- rep(0,species)

  for(node in 1:length(bet)){
    speciesBetw[out$composition[[node]]] <-
      speciesBetw[out$composition[[node]]]+bet[node]
  }
  speciesBetw
}
