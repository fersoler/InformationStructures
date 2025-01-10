#####################################################################
# Rewiring Structural Stability                                     #
# Functions to change interactions according to several conditions  #
#####################################################################

# Function to test if a matrix is diagonally dominant by rows
isDiagDomRows <- function(A){
  A1 <- abs(diag(A))          # Diagonal
  A2 <- rowSums(abs(A)) - A1  # Sum of rows out of the diagonal
  all(A1 > A2)
}

# 1. Two species u_i and u_j reduce their competition when a common cooperator
# appears. This function receives an interaction matrix and reduces competition
# among species having a common cooperator. 
rewiringReduceComp <- function(a){
  # The matrix that will be returned is A
  A <- a
  size <- dim(a)[1]
  # i is the first species
  for (i in 1:(size-1)) {
    # j is the second one
    for(j in (i+1):size){
      aij <- a[i,j]
      aji <- a[j,i]
      # in case that both species compete
      if(aij < 0 & aji < 0){
        # Competition will be divided by sumCoop
        sumCoop <- 1
        for(r in 1:size){
          # cooperation of species r with i and j
          sumCoop <- sumCoop + max(0, min(a[i,r], a[j,r]))
        }
        A[i,j] <- a[i,j]/sumCoop
        A[j,i] <- a[j,i]/sumCoop
      }
    }
  }
  A
}

# 2. Adverse abiotic conditions increase cooperation
# - 'a' is a matrix of interations
# - 'r' is the vector with intrinsic growth rates
rewiringIncreaseCoop <- function(a, r){
  # The matrix to be returned
  A <- a
  size <- dim(a)[1]
  # i is the first species
  for (i in 1:(size-1)) {
    # j is the second one
    for(j in (i+1):size){
      aij <- a[i,j]
      aji <- a[j,i]
      # in case that both species cooperate, 
      # the cooperation is increased for species
      # having negative intrinsic growth rates
      if(aij > 0 & aji > 0){
        if(r[i] < 0){
          A[i,j] <- a[i,j]/min(1,2^r[i])  
        }
        if(r[j] < 0){
          A[j,i] <- a[j,i]/min(1,2^r[j])   
        }
      }
    }
  }
  A
}

# This function is used to modify a matrix 'newA' in order 
# oldA is a diagonally dominant (rows) matrix, newA is a modification that
# could not be diagonally dominant. The modification is supposed not to 
# affect the diagonal. This funciton modulates the changes in order to keep 
# diagonal dominance. 

# Previous functions modify a matrix of interactions 'oldA' to get a new
# one, 'newA', but sometimes the diagonal dominance (by rows) may be lost. In 
# case that we want to preserve the diagonal dominance of 'oldA', this function 
# modifies 'newA' to get a diagonally dominant version. 
keepDiagDomRows <- function(oldA, newA){
  A <- newA
  if(!isDiagDomRows(newA)){
    size <- dim(newA)[1]  # size
    for(i in 1:size){
     oldRow <- oldA[i,]
     newRow <- newA[i,]
     difRow <- newRow-oldRow
     if(max(abs(difRow)) > 0){
      maxSum <- 2*abs(oldA[i,i]) - sum(abs(oldRow))
      A[i,] <- oldRow + difRow*0.99999*(maxSum/sum(abs(difRow)))
     }
    }
  }
  A
}
