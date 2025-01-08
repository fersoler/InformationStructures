# Test if the matrix is diagonally dominant considering rows
isDiagDomRows <- function(A){
  A1 <- abs(diag(A))          # Diagonal
  A2 <- rowSums(abs(A)) - A1  # Sum of rows out of the diagonal
  all(A1 > A2)
}


rewiringReduceComp <- function(a){
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
        sumCoop <- 1
        for(r in 1:size){
          sumCoop <- sumCoop + max(0, min(a[i,r], a[j,r]))
        }
        A[i,j] <- a[i,j]/sumCoop
        A[j,i] <- a[j,i]/sumCoop
      }
    }
  }
  A
}

rewiringIncreaseCoop <- function(a, r){
  A <- a
  size <- dim(a)[1]
  # i is the first species
  for (i in 1:(size-1)) {
    # j is the second one
    for(j in (i+1):size){
      aij <- a[i,j]
      aji <- a[j,i]
      # in case that both species cooperate
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

# oldA is a diagonally dominant (rows) matrix, newA is a modification that
# could not be diagonally dominant. The modification is supposed not to 
# affect the diagonal. This funciton modulates the changes in order to keep 
# diagonal dominance. 
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


#(a <- matrix(c(-1,.5,-.2,.3,-1,.7,-.4,.4,-1),3,3))
#(a2 <- matrix(c(-1,.5,-.6,.3,-1,.7,-.4,.4,-1),3,3))
#keepDiagDomRows(a,a2)


#(b <- c(1,-2,1.3))
#rewiringReduceComp(a)
#keepDiagDomRows(a, rewiringIncreaseCoop(a,b))

