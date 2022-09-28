# Building Information Structures (IS)
source("R/lemkelcp.R")


# function to get the GASS with the package LemkeLCP
getGASS_LCP_Lemke <- function(alphas, gammas, maxIter = 100){
  sol <- unlist(lemkelcp(matrix(-gammas,nrow = length(alphas)), t(-alphas), maxIter)[1])
  if(is.na(sol[1])){
    stop("LCP has no solution")
  }
  round(sol, digits=12)
}

# Building the Information Structure
# Input:
# - alphas: data frame with n alpha values
# - gammas: data matrix of n*n values
# Output: List with the following elements
# - $points: matrix with all stationary points, one point per row
# - $subsetGASS: data frame indicating the index of the GASS (2n column)
# of each subsystem (1st column)
# - $connectivity: connectivity matrix. Rows and columns are the indexes of
# the points in $points
# - $gassInd: index of the point which is the GASS of the whole system
ISbuild <- function(alphas, gammas){
  # STEP 1: Finding the GASSes for all subsystems
  size            <- length(alphas)
  # matrix to store all stationary points
  allStatPoints   <- matrix(rep(0.0, size),nrow=1)
  # dataframe indicating the GASS (ind) of each subsystem (subset)
  subsetGASS <- data.frame(subset="0", ind=1, stringsAsFactors=FALSE)
  # loop to get the GASSes for all subsystems
  for(s in 1:size){
    allSubS  <-combn(1:size,s)         # All subsystems of size s
    for(i in 1:ncol(allSubS)){
      ss     <- allSubS[,i]            # 'ss': current subsystem to look for the GASS
      # Obtaining the GASS
      gass   <- getGASS_LCP_Lemke(alphas[ss], gammas[ss,ss])  # Python Lemke
      v <- rep(0.0, size)
      v[ss] <- gass                    # v contains the GASS with 0s
      found <-  FALSE                  # was the GASS found before?
      lookAt <- nrow(allStatPoints)
      while(lookAt >= 1 & !found){     # while not found
        check <- allStatPoints[lookAt,] == v
        if(check[1]  & length(unique(check)) == 1){
          found <- TRUE
          # if the GASS was previously found, we only assign in to 'ss'
          subsetGASS[nrow(subsetGASS)+1,]<-list(subset=toString(ss), ind=lookAt)
        }
        lookAt <- lookAt -1
      }
      if(!found){
        # if the GASS was not found we add it to the list of stationary points
        allStatPoints <- rbind(allStatPoints, v)
        # and assign its number to 'ss'
        subsetGASS[nrow(subsetGASS)+1,]<-list(subset=toString(ss), ind=nrow(allStatPoints))
      }
    }
  }
  nPoints <- nrow(allStatPoints)     # number of stationary points
  rownames(allStatPoints) <- 1:nPoints

  # STEP 2: Connecting stationary points
  # Connectivity matrix
  connectivity <- matrix(0, nrow = nPoints, ncol = nPoints)
  # loop to connect all points
  for(p in 1:nPoints){  # p is the index of the point to look for its connections
    point <- allStatPoints[p,]    # stationary point
    I  <- (1:size)[point > 0]     # positions with a value  > 0
    NI <- (1:size)[point == 0]    # positions with a value == 0
    if(length(NI) != 0){          # if there are 0 values
      if(length(I) == 0){   # when all values are 0, we only look at alphas
        R <- alphas[NI]
      } else {              # otherwise, compute r_i values
        R  <- alphas[NI] + gammas[NI,I] %*% as.matrix(point[I])
      }
      J <- NI[(1:length(R))[R>0]] # J is the i's with r_i > 0
      if(length(J) > 0){          # if there are values in J look for connections
        for(t in 1:length(J)){    # subsets of J of size t
          if(length(J) == 1){
            combJ <- as.matrix(J,ncol=1)
          } else {
            combJ <- combn(J,t)
          }
          for(cn in 1:ncol(combJ)){
            cj <- combJ[,cn]      # take each of them (cj)
            K <- sort(c(cj,I))    # K = cj U I
            # indK is the index of the point GASS(K)
            indK <- subsetGASS[subsetGASS$subset==toString(K),][1,2]
            # set a connection from p to indK
            connectivity[p, indK] <- 1
          }
        }
      }
    }
  }
  gassInd <- subsetGASS[subsetGASS$subset==toString(1:size),][1,2]
  list(
    # The matrix with all stationary points: one point per row
    points       = allStatPoints,
    # The dataframe indicating the index of the GASS (2n column)
    # of each subsystem (1st column)
    subsetGASS   = subsetGASS,
    # The connectivity matrix. Rows and columns are the indexes of
    # the points in $points
    connectivity = connectivity,
    # Index of the point which is the GASS of the whole system
    gassInd      = gassInd
  )
}

