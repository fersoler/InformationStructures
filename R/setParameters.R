library(tidyverse)
library(deSolve)


######################################################################
# Functions to set LV parameters
######################################################################

# The following function sets the parameters of a gLV system using the 
# 'Log Integral method' in: 
# P.H. Kloppers, J.C. Greeff,
# "Lotka–Volterra model parameter estimation using experiential data",
# Applied Mathematics and Computation, Volume 224, 2013, Pages 817-825
# https://doi.org/10.1016/j.amc.2013.08.093.

# The method has been modified to allow manual setting of some parameters

# vals: timeseries with abundances. Each species data is in a column. Abundances
# are measured in regular intervals (rows). The number of species is equal to
# the number of columns of this matrix. 
# intra: this argument allows for manual setting of some or all the intrinsic 
# growth rates. It is NULL by default, so that the function sets them all. 
# Otherwise it will be a vector with numeric values in the positions 
# corresponding to parameters that are manually adjusted and 'NA' otherwise. 
# inter: similar argument for interaction rates between species and of a 
# species with itself. When it is not NULL, it is a square matrix with 
# numerical values in the parameters that are manually adjusted and NA for 
# those adjusted by the function. The values on the diagonal correspond to the 
# interaction of each species with itself. The value in [i,j] is the effect 
# that species i gets from j. 
# limits: this argument is used to take only certain rows of the X and L 
# matrices to adjust the parameters. That is, only a certain time period 
# is taken into account. It allows to adjust parameters by only looking
# a window around a certain instant. 
gLV_params_log <- function(vals, intra = NULL, inter = NULL, limits = NULL){
  size <- ncol(vals)
  b <- rep(NA,size)            # For intra parameters
  A <- matrix(0, size, size)   # For inter parameters
  L <- diff(log(vals))
  X <- as.matrix(cbind(rep(1,nrow(L)), 
             (head(vals,-1)+tail(vals,-1))/2), ncol=1+size)
  for(sp in 1:size){
    Lsp <- L[,sp]
    Xsp <- X
    # In case some a_sp_j is specified
    if(!is.null(inter)){
      colsSel <- Xsp[,which(!is.na(inter[sp,]))+1, drop = FALSE]
      paramsSel <- inter[sp,!is.na(inter[sp,])]
      Lsp <- Lsp-rowSums(scale(colsSel, center=FALSE, 
                               scale=1/paramsSel))
      Xsp <- Xsp[,c(1,which(is.na(inter[sp,]))+1), drop = FALSE]
      A[sp,] <- inter[sp,]
    }
    # In case the b_sp is specified
    if(!is.null(intra) && !is.na(intra[sp])){
      Lsp <- Lsp-intra[sp]
      Xsp <- Xsp[,-1, drop = FALSE]
      b[sp] <- intra[sp]
    }
    if(!is.null(limits)){
      Xsp  <- Xsp[limits,]
      Lsp  <- Lsp[limits]
    }
    Xspt <- t(Xsp)
    solut  <- solve(Xspt %*% Xsp) %*% Xspt %*% Lsp
    if(is.na(b[sp])){
      b[sp] <- solut[1]
    }
    if(!is.null(inter)){
      A[sp,which(is.na(A[sp,]))] <- tail(solut,-1)
    } else {
      A[sp, ] <- tail(solut,-1)
    }
  }
  list(b = b, A = A)
}

# Function to set the parameters of a gLV model such that:
# - The 'inter' matrix is Volterra-Lyapunov stable
# - Absolute value of the matrix is diagonally dominant.
# - Intrinsic grow rates are different at each instant. They are calculated
# by looking at each time point and 'windW' points at the left and the right 
# of it. The greater value of 'windW' the more stable is the evolution of 
# intrinsic growth rates but reconstruction of the time series is worse. 
getLVparamsNonAuto <- function(vals, windW){
  size <- ncol(vals)
  #gammas <- -diag(size)     # Values of diagonal are -1
  #gammas[gammas == 0] <- NA # Other are not specified
  alphas <- matrix(0, ncol=size, nrow=nrow(vals)-1)
  # Set A parameters (matrix)
  #gammas <- gLV_params_log(vals, inter = gammas, intra = NULL)$A
  gammas <- gLV_params_log(vals, inter = NULL, intra = NULL)$A
  # Make the matrix diagonal dominant 
  # for(sp in 1:size){
  #   sumGammas <- sum(abs(gammas[sp,]))-1
  #   if(sumGammas>=1){
  #     gammas[sp,] <- 0.995*(gammas[sp,]/sumGammas)
  #     gammas[sp,sp] <- -1
  #   }
  # }
  # Make the matrix Volterra-Lyapunov stable:
  S <- (gammas+t(gammas))/2
  routh <- max(eigen(S)$values)
  if(routh > 0){
    gammas <- gammas - (routh+0.01)*diag(size)
  }
  # Calculate intrinsic growth rates
  for(a in 1:nrow(alphas)){
    # Set limits
    Linit <- if(a<windW) 1 else a-windW
    Lend <- if(a+windW>nrow(alphas)) nrow(alphas) else a+windW
    # Call the function to set parameters by setting the matrix and limits: 
    alphas[a,] <- gLV_params_log(vals, inter = gammas, intra = NULL, limits=Linit:Lend)$b
  }
  list(alphas = alphas, gammas = gammas)
}

######################################################################
# Solve a 3 species gLV system
######################################################################

# Equations of the gLV system for 3 species
LotkaVolterra3 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dy1 <-  y1 * (a1 + g11*y1 + g12*y2 + g13*y3)
    dy2 <-  y2 * (a2 + g21*y1 + g22*y2 + g23*y3)
    dy3 <-  y3 * (a3 + g31*y1 + g32*y2 + g33*y3)
    list(c(dy1, dy2, dy3))
  })
}

# Solve the evolution of a 3 species gLV system given a starting point
# alphas: intrinsic growth rates
# gammas: matrix of interactions
# initStat: (x0,y0,z0) starting point
# lengTime: parameter to adjust the number of points in the solution. The 
# higher the number, the more accurate the simulation.
# timeI: time at beginning
# timeE: time at end
getEvolLV3 <- function(alphas, gammas, initStat, lengTime, timeI, timeE){
  pars <- c(a1 = alphas[1], a2 = alphas[2], a3 = alphas[3],
            g11 = gammas[1,1], g12 = gammas[1,2],
            g13 = gammas[1,3], g21 = gammas[2,1],
            g22 = gammas[2,2], g23 = gammas[2,3],
            g31 = gammas[3,1], g32 = gammas[3,2],
            g33 = gammas[3,3])
  yI <- as.numeric(initStat)
  y0 <- c(y1 = yI[1], y2 = yI[2], y3 = yI[3])
  sol <- ode(y = y0, times = seq(timeI,timeE,length.out = lengTime),parms = pars, 
             func = LotkaVolterra3, atol = 1e-12)
  sol
}

# For 2 species
LotkaVolterra2 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dy1 <-  y1 * (a1 + g11*y1 + g12*y2)
    dy2 <-  y2 * (a2 + g21*y1 + g22*y2)
    list(c(dy1, dy2))
  })
}

getEvolLV2 <- function(alphas, gammas, initStat, lengTime, timeI, timeE){
  pars <- c(a1 = alphas[1], a2 = alphas[2],
            g11 = gammas[1,1], g12 = gammas[1,2],
            g21 = gammas[2,1], g22 = gammas[2,2])
  yI <- as.numeric(initStat)
  y0 <- c(y1 = yI[1], y2 = yI[2])
  sol <- ode(y = y0, times = seq(timeI,timeE,length.out = lengTime),parms = pars, 
             func = LotkaVolterra2, atol = 1e-12)
  sol
}

# New functions (May 4th 2024)

# Set A matrix assuming b values have minimal changes
gLV_set_A_log <- function(vals){
  size <- ncol(vals)
  A <- matrix(0, size, size)   # For inter parameters
  L1 <- diff(log(vals))
  L <- diff(L1)
  
  X1 <- as.matrix((head(vals,-1)+tail(vals,-1))/2, ncol=size)
  X <- diff(X1)
  
  for(sp in 1:size){
    Lsp <- L[,sp]
    Xsp <- X
    Xspt <- t(Xsp)
    A[sp, ]  <- solve(Xspt %*% Xsp) %*% Xspt %*% Lsp
  }
  A
}

# Once A is set, get b parameters given a window. 
get_b_given_A <- function(vals, gammas, windW){
  size <- ncol(vals)
  alphas <- matrix(0, ncol=size, nrow=nrow(vals)-1)
  # Calculate intrinsic growth rates
  for(a in 1:nrow(alphas)){
    # Set limits
    Linit <- if(a<windW) 1 else a-windW
    Lend <- if(a+windW>nrow(alphas)) nrow(alphas) else a+windW
    # Call the function to set parameters by setting the matrix and limits: 
    alphas[a,] <- gLV_params_log(vals, inter = gammas, intra = NULL, limits=Linit:Lend)$b
  }
  alphas
}

# LV_map with A(t) B(t)
LV_map <- function(vals, theta){
  size <- ncol(vals)
  times <- nrow(vals)  
  # Distances
  distances <- as.matrix(dist(vals[1:(times-1),], diag=TRUE, upper = TRUE))
  # Mean distance
  mean_d <- mean(distances)
  
  L <- diff(log(vals))
  X <- cbind(1,as.matrix((head(vals,-1)+tail(vals,-1))/2, ncol=size))

  A_list <- list()
  B_all <- matrix(0, times-1, size)
  
  for(t in 1:(times-1)){
    A <- matrix(0, size, size)   # For inter parameters
  #  B <- matrix(0, size) # For intrinsic growths
    
    w_t <- exp(-theta*(distances[t,])/mean_d)
    Xsp <- X*w_t
    
    for(sp in 1:size){
      
      Lsp <- L[,sp]*w_t
      
      Xspt <- t(Xsp)
      solut <- solve(Xspt %*% Xsp) %*% Xspt %*% Lsp
      
      B_all[t,sp] <- solut[1]
      A[sp,] <- solut[2:(size+1)] 
    }
    
    A_list[[length(A_list)+1]] <- A
    #B_list[[length(B_list)+1]] <- B
  }
  
  list(A = A_list, B = B_all)
}


# LV_map with constant A
LV_map_A_const <- function(vals, theta){
  size <- ncol(vals)
  times <- nrow(vals)  
  # Distances
  distances <- as.matrix(dist(vals[1:(times-1),], diag=TRUE, upper = TRUE))
  # Distances
  # print(dim(distances))
  # Mean distance
  mean_d <- mean(distances)
  
  A <- matrix(0, size, size)   # For inter parameters
  B <- matrix(0, size, times-1) # For intrinsic growths
  L <- diff(log(vals))
  X <- cbind(1,as.matrix((head(vals,-1)+tail(vals,-1))/2, ncol=size))
  
  for(sp in 1:size){
    #print(sp)
    X_exp <- NULL
    L_exp <- NULL
    for(t in 1:(times-1)){
      w_sp <- exp(-theta*(distances[t,])/mean_d)
      Lsp <- L[,sp]*w_sp
      Xsp <- X*w_sp
      Xsp2 <- matrix(0, nrow(Xsp), times+size-1)
      Xsp2[,t] <- Xsp[,1]
      Xsp2[,times:(times+size-1)] <- Xsp[,2:(size+1)]
      L_exp <- c(L_exp, Lsp)
      X_exp <- rbind(X_exp,Xsp2)
    }
    
    Xspt <- t(X_exp)
    solut <- solve(Xspt %*% X_exp) %*% Xspt %*% L_exp
    
    B[sp,] <- solut[1:(times-1)]
    A[sp,] <- solut[times:(times+size-1)]
    
  }
  
  list(A = A, B = B)
}


####################################################################
# Idea 31 may.
# Ponderar cada punto según las desviaciones estándar respecto a la
# media, bien de valores o de incrementos. 
####################################################################


gLV_params_log_weighted_A <- function(vals, intra = NULL, inter = NULL, limits = NULL){
  size <- ncol(vals)
  b <- rep(NA,size)            # For intra parameters
  A <- matrix(0, size, size)   # For inter parameters
  L <- diff(log(vals))
  X <- as.matrix(cbind(rep(1,nrow(L)), 
                       (head(vals,-1)+tail(vals,-1))/2), ncol=1+size)
  
  # Weights
  #L_abs_dif <- abs(diff(vals))
  L_abs_dif <- abs(L)
  L_log_mean <- log(apply(L_abs_dif,2,mean))
  XW <- abs(sweep(log(L_abs_dif), 2, L_log_mean, "-"))
  
#  sdL <- apply(abs(L),2,mean)
#  XW <- sweep(abs(L), 2, sdL, "/")
  
  for(sp in 1:size){
    Lsp <- L[,sp]/(XW[,sp])                 # Weighted
    Xsp <- sweep(X, 1, (XW[,sp]), "/")      # Weighted 
    # In case some a_sp_j is specified
    if(!is.null(inter)){
      colsSel <- Xsp[,which(!is.na(inter[sp,]))+1, drop = FALSE]
      paramsSel <- inter[sp,!is.na(inter[sp,])]
      Lsp <- Lsp-rowSums(scale(colsSel, center=FALSE, 
                               scale=1/paramsSel))
      Xsp <- Xsp[,c(1,which(is.na(inter[sp,]))+1), drop = FALSE]
      A[sp,] <- inter[sp,]
    }
    # In case the b_sp is specified
    if(!is.null(intra) && !is.na(intra[sp])){
      Lsp <- Lsp-intra[sp]
      Xsp <- Xsp[,-1, drop = FALSE]
      b[sp] <- intra[sp]
    }
    if(!is.null(limits)){
      Xsp  <- Xsp[limits,]
      Lsp  <- Lsp[limits]
    }
    Xspt <- t(Xsp)
    solut  <- solve(Xspt %*% Xsp) %*% Xspt %*% Lsp
    if(is.na(b[sp])){
      b[sp] <- solut[1]
    }
    if(!is.null(inter)){
      A[sp,which(is.na(A[sp,]))] <- tail(solut,-1)
    } else {
      A[sp, ] <- tail(solut,-1)
    }
  }
  list(b = b, A = A, XW = XW)
}

getLVparamsNonAuto_weighted_A <- function(vals, windW){
  size <- ncol(vals)
  #gammas <- -diag(size)     # Values of diagonal are -1
  #gammas[gammas == 0] <- NA # Other are not specified
  alphas <- matrix(0, ncol=size, nrow=nrow(vals)-1)
  # Set A parameters (matrix)
  #gammas <- gLV_params_log(vals, inter = gammas, intra = NULL)$A
  params <- gLV_params_log_weighted_A(vals, inter = NULL, intra = NULL)
  gammas <- params$A
  # Make the matrix diagonal dominant 
  # for(sp in 1:size){
  #   sumGammas <- sum(abs(gammas[sp,]))-1
  #   if(sumGammas>=1){
  #     gammas[sp,] <- 0.995*(gammas[sp,]/sumGammas)
  #     gammas[sp,sp] <- -1
  #   }
  # }
  # Make the matrix Volterra-Lyapunov stable:
  S <- (gammas+t(gammas))/2
  routh <- max(eigen(S)$values)
  if(routh > 0){
    gammas <- gammas - (routh+0.01)*diag(size)
  }
  # Calculate intrinsic growth rates
  for(a in 1:nrow(alphas)){
    # Set limits
    Linit <- if(a<windW) 1 else a-windW
    Lend <- if(a+windW>nrow(alphas)) nrow(alphas) else a+windW
    # Call the function to set parameters by setting the matrix and limits: 
    alphas[a,] <- gLV_params_log(vals, inter = gammas, intra = NULL, limits=Linit:Lend)$b
  }
  list(alphas = alphas, gammas = gammas, XW = params$XW)
}
