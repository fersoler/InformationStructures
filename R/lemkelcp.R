library(R6)

lemketableau <- R6::R6Class(

  classname = "lemketableau",

  public = list(
    Tmat = NULL,
    n    = NULL,
    wPos = NULL,
    zPos = NULL,
    W = NULL,
    Z = NULL,
    Y = NULL,
    Q = NULL,
    Tind    = NULL,
    maxIter = NULL,

    initialize = function(M, q, maxIter = 100) {
      n <- length(q)
      self$Tmat <- cbind(diag(n), -M, -rep(1,n), q)
      self$n    <- n
      self$wPos <- 1:n
      self$zPos <- n+(1:n)
      self$W <- 0
      self$Z <- 1
      self$Y <- 2
      self$Q <- 3
      TbInd  <- rbind(self$W * rep(1,n), 1:n)
      TnbInd <- rbind(self$Z * rep(1,n), 1:n)
      DriveInd     <- c(self$Y, 1)
      QInd         <- c(self$Q, 1)
      self$Tind    <- cbind(TbInd, TnbInd, DriveInd, QInd)
      self$maxIter <- maxIter
    },

    lemkeAlgorithm = function(){
      initVal <- self$init()
      if(!initVal){
        return(list(z = rep(0,self$n), code = 0,
                    str = 'Solution found after 0 iterations'))
      } else {
        for(k in 1:self$maxIter){
          stepVal <- self$step()
          if(self$Tind[1,ncol(self$Tind)-1] == self$Y){
            return(list(z = self$extractSolution(), code = 0,
                        str = paste('Solution found after',k,'iterations')))
          } else {
            if(!stepVal){
              return(list(z = NA, code = 1, str = 'Secondary ray found'))
            }
          }
        }
        return(list(z = NA, code = 2, str = 'Max iterations exceeded'))
      }
    },

    init = function(){
      q <- self$Tmat[,ncol(self$Tmat)]
      minQ <- min(q)
      if(minQ < 0){
        ind <- match(minQ, q)
        self$clearDriverColumn(ind)
        self$pivot(ind)
        TRUE
      } else {
        FALSE
      }
    },

    step = function(){
      q <- self$Tmat[,ncol(self$Tmat)]
      a <- self$Tmat[,ncol(self$Tmat)-1]
      ind <- NA
      minRatio <- Inf
      for(i in 1:self$n){
        if(a[i]>0){
          newRatio = q[i]/a[i]
          if(newRatio < minRatio){
            ind <- i
            minRatio <- newRatio
          }
        }
      }
      if(minRatio < Inf){
        self$clearDriverColumn(ind)
        self$pivot(ind)
        TRUE
      } else {
        FALSE
      }
    },

    extractSolution = function(){
      z <- rep(0,self$n)
      q <- self$Tmat[,ncol(self$Tmat)]
      for(i in 1:self$n){
        if(self$Tind[1,i] == self$Z){
          z[self$Tind[2,i]] <- q[i]
        }
      }
      z
    },

    partnerPos = function(pos){
      v   <- self$Tind[1,pos]
      ind <- self$Tind[2,pos]
      if(v == self$W){
        ppos <- self$zPos[ind]
      } else if(v == self$Z){
        ppos <- self$wPos[ind]
      } else {
        ppos <- NA
      }
      ppos
    },

    pivot = function(pos){
      ppos <- self$partnerPos(pos)
      if(!is.na(ppos)){
        self$swapColumns(pos,ppos)
        self$swapColumns(pos,ncol(self$Tind)-1)
        TRUE
      } else {
        self$swapColumns(pos,ncol(self$Tind)-1)
        FALSE
      }
    },

    swapMatColumns = function(M, i, j){
      Mi    <- M[,i]
      Mj    <- M[,j]
      M[,i] <- Mj
      M[,j] <- Mi
      M
    },

    swapPos = function(v,ind,newPos){
      if(v == self$W){
        self$wPos[ind] <- (newPos %% (2*self$n+2))
      } else if(v == self$Z){
        self$zPos[ind] <- (newPos %% (2*self$n+2))
      }
    },

    swapColumns = function(i, j){
      iInd <- self$Tind[,i]
      jInd <- self$Tind[,j]
      v    <- iInd[1]
      ind  <- iInd[2]
      self$swapPos(v,ind,j)
      v   <- jInd[1]
      ind <- jInd[2]
      self$swapPos(v,ind,i)
      self$Tind <- self$swapMatColumns(self$Tind,i,j)
      self$Tmat <- self$swapMatColumns(self$Tmat,i,j)
    },

    clearDriverColumn = function(ind){
      a <- self$Tmat[ind,ncol(self$Tmat)-1]
      self$Tmat[ind,] <- self$Tmat[ind,]/a
      for(i in 1:self$n){
        if(i != ind){
          b <- self$Tmat[i,ncol(self$Tmat)-1]
          self$Tmat[i,] <- self$Tmat[i,] - (b * self$Tmat[ind,])
        }
      }
    }
  )
)

# Main function to solve the LCP problem given M and q using Lemke's algorithm
lemkelcp <- function(M,q,maxIter=100){
  tableau <- lemketableau$new(M,q,maxIter)
  tableau$lemkeAlgorithm()
}
