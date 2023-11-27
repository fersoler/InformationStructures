#################################################################
# Functions to draw Information Fields
# Author: Fernando Soler-Toscano - fsoler@us.es
#################################################################

library(plotly)
library(pracma)
library(tidyverse)
library(shape)
library(randomcoloR)
library(RColorBrewer)
library(deSolve)
source("R/ISbuild.R")
source("R/figs/sphereCones.R")

# Function to set the distances of a given point to all points and
# edges (segments joining two vertices) of an IS.
# The output is a matrix. Values in [i,i] indicate the distance of
# 'point' to vertex number i of the IS. Values in other [i,j] are
# the distance of 'point' to the segments joining IS nodes i,j.
distancesPointToIS <- function(IS, point){
  npoints <-nrow(IS$points)     # Number of IS states
  distances <- IS$connectivity  # Matrix to store the distances
  for(a in 1:npoints){
    for(b in 1:npoints){
      # There is a link from a to b (or a=b, so it's a point)
      if(distances[a,b] == 1 || a == b){
        pa <- IS$points[a,]
        pb <- IS$points[b,]
        # Distance from 'point' to the segment pa--pb. In case
        # pa==pb, it's the distance from 'point' to IS node 'pa'.
        distances[a,b] <- segm_distance(pa,pb,point)$d
      } else {
        # 0s in the connectivity matrix are replaced by NA
        distances[a,b] <- NA
      }
    }
  }
  # Return the matrix of distances
  distances
}

# Lyapunov function. Valid for 2 species.
# Input: IS for 2 species and a point.
# Output: Value of the Lyapunov Function at the given point.
lyapunovFunc <- function(IS, point){
  # Matrix with distances from all nodes (diagonal) and edges of the IS
  # to the given point
  distances <- distancesPointToIS(IS, point)
  # GASS index
  gi <- IS$gassInd
  # Index of point 0
  i0 <- 1 # default with current IS build functions
  # Number of points
  npoints <- nrow(IS$points)
  # There are different cases to compute the Lyapunov function
  # 1. Trivial IS (only 0 point)
  if(npoints == 1){
    return(distances[i0,i0])
  }
  # 2. IS with two points: (00 -> 10) or (00 -> 01)
  if(npoints == 2){
    d1 <- distances[gi,gi]
    d2 <- distances[i0,i0]
    return(d1/(d1+d2))
  }
  # 3. Three points: (00 -> 10 -> 01) or (00 -> 01 -> 10)
  if(npoints == 3){
    d1 <- distances[gi,gi]
    d2 <- min(distances[c(-gi),c(-gi)], na.rm = TRUE)
    z1 <- d1 / (d1 + d2)

    d3 <-  min(distances[c(-i0),c(-i0)], na.rm = TRUE)
    d4 <- distances[i0,i0]
    z2 <- d3 / (d3 + d4)

    return((z1+z2)/2)
  }

  # 4. Four points (00 -> 01,10 -> 11)
  d1 <- distances[gi,gi]
  d2 <- min(distances[c(-gi),c(-gi)], na.rm = TRUE)
  z1 <- d1 / (d1 + d2)

  intermPoints <- c(1:npoints)[c(-i0,-gi)]
  d3a <- distances[intermPoints[1],gi]
  d4a <- distances[i0,intermPoints[2]]
  d3b <- distances[intermPoints[2],gi]
  d4b <- distances[i0,intermPoints[1]]
  z2a <- d3a / (d3a + d4a)
  z2b <- d3b / (d3b + d4b)
  z2  <- min(z2a,z2b)

  d5 <-  min(distances[c(-i0),c(-i0)], na.rm = TRUE)
  d6 <- distances[i0,i0]
  z3 <- d5 / (d5 + d6)

  return((z1+z2+z3)/3)
}


# Function to plot the Information Field of a given
# IS with 2 species. The optional argument 'np' sets
# the divisions of the grid. The higher value, the
# best resolution but taking more time.
plotLyapFunc <- function(IS, np = 100){
  # Set x and y points:
  # Values are from 0 to 1.3 times the maximum in the IS points. The
  # coordinates of IS points are added for better visualization.
  x <- unique(sort(c(seq(0, max(c(1,IS$points[,1]))*1.3, length.out = np),
                     array(IS$points[,1]))))
  y <- unique(sort(c(seq(0, max(c(1,IS$points[,2]))*1.3, length.out = np),
                     array(IS$points[,2]))))

  # Value of the Lyapunov function for all (x,y):
  zetas <- matrix(0,length(x), length(y))
  for(i in 1:length(x)){
    for(j in 1:length(y)){
      zetas[i,j] <- lyapunovFunc(IS,c(x[i],y[j]))
    }
  }
  # Transpose for plotting:
  z <- t(zetas)

  # IS points:
  zetas <- rep(0, nrow(IS$points))
  for(r in 1:nrow(IS$points)){
    zetas[r] <- lyapunovFunc(IS,c(IS$points[r,1], IS$points[r,2]))
  }

  # Figure with the plot of the Lyapunov function:
  kd <- with(MASS::geyser, MASS::kde2d(duration, waiting, n = 50))
  fig <- plot_ly(x = x, y = y, z = z,
                 showscale = FALSE,
                 showticklabels = FALSE,
                 colorscale= 'Viridis',
                 reversescale = TRUE,
                 showgrid = TRUE
  ) %>% add_surface()
  # Add marks for the IS points:
  fig <- fig %>%
    add_trace(
      type = "scatter3d",
      mode = "markers",
      x = array(IS$points[,1]),
      y = array(IS$points[,2]),
      z = zetas,
      showlegend = FALSE,
      marker = list(color = 'black', size = 3)
    )
  fig <- fig %>%
    layout(
      scene = list(
        xaxis = list(title = "u1", showticklabels = FALSE),
        yaxis = list(title = "u2", showticklabels = FALSE),
        zaxis = list(title = "", showticklabels = FALSE)
      )
    )
  fig
}

# Function to display the Lyapunov function with the trajectories
# starting at some given points. Arguments:
# - IS. An Information Structure
# - g: connectivity matrix
# - a: intrinsic growt rates
# - initPoints: dataframe with initial points (cols 1 and 2) and 
# colors (col 3)
# - time: time for ode solutions. 
plotLyapFuncWithInitPoints <- function(IS, g, a, initPoints, time, np = 100){
  fig <- plotLyapFunc(IS, np)
  for (nInPo in 1:nrow(initPoints)) {
    ## Add the lines representing the evolution
    u10 <- initPoints[nInPo,1]
    u20 <- initPoints[nInPo,2]
    
    pars <- c(a1 = as.numeric(a[1]), a2 = as.numeric(a[2]), g12 = g[1,2], g21 = g[2,1])
    times <- seq(0,time*3,length.out = 2)
    
    y0 <- c(y1 = u10, y2 = u20)
    
    sol2 <- ode(y = y0, times = seq(0,time,length.out = 100), parms = pars, func = LotkaVolterra, atol = 1e-12)
    
    zetas <- mapply(function(x, y) lyapunovFunc(IS, c(x, y)), sol2[, "y1"], sol2[, "y2"])
    
    fig <- fig %>%
      add_trace(
        type = "scatter3d",
        mode = "lines",
        x = sol2[, "y1"],
        y = sol2[, "y2"],
        z = rep(0, length(sol2[, "y1"])),
        line = list(color = initPoints[nInPo,3], width = 4),
        showlegend = FALSE
      )
    
    fig <- fig %>%
      add_trace(
        type = "scatter3d",
        mode = "lines",
        x = sol2[, "y1"],
        y = sol2[, "y2"],
        z = zetas,
        line = list(color = initPoints[nInPo,3], width = 6),
        showlegend = FALSE
      )
    
    for(i in seq(1, length(sol2[, "y1"]), 20)){
      fig <- fig %>% add_trace(
        x = rep(sol2[i, "y1"], 2), y = rep(sol2[i, "y2"], 2),
        z = c(zetas[i],0),
        type = "scatter3d",
        mode = "lines",
        color = I("grey"),
        line = list(color = 'grey'), showlegend = F)
    }
  }
  fig
}

# Given 'alpha' and 'gamma' values, this function returns the
# endpoint of a vector starting at 'point' with the given
# 'size' and the direction given by gLV equations.
getDirection <- function(alphas, gammas, point, size){
  u1 <- point[1]
  u2 <- point[2]

  al <- unlist(alphas)

  e1 = u1 + u1*(al[1] - u1 + gammas[1,2]*u2)
  e2 = u2 + u2*(al[2] - u2 + gammas[2,1]*u1)
  end <- c(e1, e2)

  vec <- c(end[1]-u1, end[2]-u2)
  length <- sqrt(sum(vec^2))
  unit_vec <- vec / length
  new_vec <- size * unit_vec

  c(u1 + new_vec[1], u2 + new_vec[2])
}

# Equations of the gLV system
LotkaVolterra <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dy1 <-  y1 * (a1 - y1 + g12*y2)
    dy2 <-  y2 * (a2 - y2 + g21*y1)
    list(c(dy1, dy2))
  })
}

# Function to draw a solution of gLV with a line of the given 'color'
# and by drawing arrows at intervals separated by 'disArrow'
drawLineArrow <- function(sol, color, disArrow){
  # Accumulated distance
  acc <- 0
  for(i in 1:(nrow(sol)-1)){
    # Distance of the next segment
    d <- segm_distance(sol[i,2:3],sol[i,2:3],sol[i+1,2:3])$d
    acc <- acc+d
    if(acc>=disArrow){
      # Paint an arrow
      Arrows(
        x0 = sol[i,2],
        y0 = sol[i,3],
        x1 = sol[i+1,2],
        y1 = sol[i+1,3],
        arr.length = 0.4,
        arr.col = color,
        arr.type = 'triangle'
      )
      acc <- 0
    } else {
      # Draw the segment
      lines(sol[i:(i+1),'y1'], sol[i:(i+1),'y2'], type = "l", lty = 1,
            col = color, lwd = 2)
    }
  }
}

# This function returns a string which identifies the scheme of
# an IS. The first character 's' is the number of species. The second
# 'c' is the number of communities. Then there are 'c*s' characters
# indicating the presence of each specie in the communities and
# c*c with the connectivity. In order to have a 1-to-1 correspondence
# of IS and strings, the same algorithm to build IS has to be
# used for all strings.
ISschemeToString <- function(IS){
  theP <- c(t(IS$points))
  theP[theP>0]<-1
  theP <- c(ncol(IS$points),nrow(IS$points),theP,t(IS$connectivity))
  paste(theP,collapse="")
}

# Function to get the intermediate point between 'init' and 'end' which
# is closer to 'end' but the target color is 'targetCol'.
getBorderPoint <- function(init, end, targetCol, gammas, minD){
  d <- euclidean(init,end)
  if(d<=minD){
    return(init)
  } else {
    mp <- midPoint(init,end)
    ISmp <- ISbuild(
      as.data.frame(
        matrix(round(mp,5),1,3)),gammas)
    # Get IS type (color)
    ISstr <- ISschemeToString(ISmp)
    if(ISstr == targetCol){
      return(getBorderPoint(mp, end, targetCol, gammas,minD))
    } else {
      return(getBorderPoint(init, mp, targetCol, gammas,minD))
    }
  }
}



# Auxiliary functions to get the euclidean distance between two points
# and the midPoint vetween two vectors
euclidean <- function(p1,p2) sqrt(sum((p1-p2)^2))
midPoint <- function(p1,p2) (p1+p2)/2
# Normalize a vector
normVec <- function(vec){
  vec/sqrt(sum(vec^2))
}