library(plotly)
library(pracma)
library(tidyverse)
library(shape)
library(randomcoloR)
library(RColorBrewer)
library(deSolve)
source("R/ISbuildImproved.R")
source("figures/sphereCones.R")

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
  # # Condition to draw an arrow at the end
  # if(d>disArrow/2){
  #   i <- nrow(sol)-1
  #   Arrows(
  #     x0 = sol[i,2],
  #     y0 = sol[i,3],
  #     x1 = sol[i+1,2],
  #     y1 = sol[i+1,3],
  #     arr.length = 0.4,
  #     arr.col = color,
  #     arr.type = 'triangle'
  #   )
  # }
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
    ISmp <- ISbuildThird(
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


# Function to build the 3D mesh (points and triangles) that can be used to
# display the sub-cones of a gLV system with 3 specied.
# Input:
# - gammas, interactions between species. The diagonal is -1, negative
# values are competition and possitive values are cooperation.
# - cone: array indicating the cone to be displayed, c(1,1,0) for example
# is the cone where only the first two species appear.
# - resol: initial number of division of the cone into triangles. Good values
# go from 80 to 150.
# - minD: distance to stop dividing the triangles. The lower value, the more
# times the function takes but the better resolution of the output.
sphereSubCones <- function(gammas, cone, resol, minD){
  # Initial mesh of triangles
  rawMesh <- BuildMeshPattern(resol)
  # Points to instantiate the mesh to the actual cone
  thePointsRaw  <- setPointsCone(-gammas,cone,rawMesh$vb)
  # Triangles of the mesh
  theTriang     <- as.data.frame(rawMesh$it)
  # A column is added with the kind of IS that will be converted into the color
  theTriang$col <- "no"  # not color yet
  # Dataframe with the points:
  thePoints <- data.frame(a1  = thePointsRaw[1,],
                          a2  = thePointsRaw[2,],
                          a3  = thePointsRaw[3,],
                          col = "0"
  )

  # The color of each point is set:
  for(p in 1:nrow(thePoints)){
    # build the IS of the current points
    ISp <- ISbuildThird(as.data.frame(
      matrix(round(as.numeric(thePoints[p,1:3]),5),1,3)),gammas)
    # Get the GASS
    gass <- ISp$points[ISp$gassInd,]
    # Check is the GASS is in the intented cone.
    # It may happen that some points at the border of the cone
    # produce IS out of the cone by rounding issues.
    nonZerosCone <- c(1:3)[gass > 0]
    nonZerosGass <- c(1:3)[cone > 0]
    # Is the IS is in the cone
    if(identical(nonZerosCone,nonZerosGass)){
      # Register in the dataframe of points the kind of IS
      thePoints[p,4] <- ISschemeToString(ISp)
    } else {
      # Otherwise mark as 0
      thePoints[p,4] <- "0"
    }
  }

  # Now we go through all triangles setting the color if possible or
  # creating new triangles. The general procedure is:
  # - If all vertices of the triangle have the same color, it's the color
  # of the triangle.
  # - Otherwise:
  #   * If the highest distance between vertices is lower than 'minD', the
  # color of the triangle is set to the color of the center when the three
  # vertices have different colors. When two vertices have the same color and
  # the third one is different, the triangle is divided into 3 sub-triangles
  # and colors are set, trying to find the division of the colors.
  #   * Else, the triangle is divided into 4 new triangles.

  # Total number of triangles. We're going to go across all triangles, so the
  # function stops when a counter 'curntTriangle' indicating the triangle we're
  # exploring reaches totalTriangles. But totalTriangles increases as new
  # triangles are added to the mesh because of the division of previous
  # triangles.
  totalTriangles <- nrow(theTriang)
  curntTriangle  <- 1
  # And this is the lastPoint of the mesh (total number of points)
  lastPoint <- nrow(thePoints)

  # Loop that goes across all triangles
  while (curntTriangle <= totalTriangles) {

    # If the color of the triangle is set, then
    # go to the following triangle
    if(theTriang[curntTriangle,4] != "no"){
      curntTriangle <- curntTriangle+1
    }

    # Vertex numbers of the triangle (in the mesh)
    triVert <- unlist(theTriang[curntTriangle,1:3])

    # Points of the triangle
    triPoin <- thePoints[triVert,]

    # Colors of the points in the triangle
    theCols <- triPoin[,4]
    # Different colors
    difCols <- unique(theCols)

    # If there are different colors
    if(length(difCols) > 1){
      # Distances between triangle vertices
      d1 <- pdist2(unlist(triPoin[1,1:3]), unlist(triPoin[2,1:3]))
      d2 <- pdist2(unlist(triPoin[1,1:3]), unlist(triPoin[3,1:3]))
      d3 <- pdist2(unlist(triPoin[2,1:3]), unlist(triPoin[3,1:3]))

      # In case the triangle is too small:
      # - The color of the middle point is selected in case that
      # the three colors are different
      # - When there are two different colors, the division line
      # is searched.
      if(max(d1,d2,d3) <= minD){
        # If there are two different colors, we divide the
        # triangle into three parts. One of the different color (trying to
        # find the 'line' dividing the colors) and the other two of the
        # color shared by two vertices.
        if(length(difCols) == 2){
          # Get the point with the different color and the two
          # other points with equal color
          if(identical(triPoin[1,4], triPoin[2,4])){
            pDf <- 3 # Index of the point with different color
            pE1 <- 1 # Points with equal color
            pE2 <- 2
          } else { # Other cases
            if(identical(triPoin[1,4], triPoin[3,4])){
              pDf <- 2
              pE1 <- 1
              pE2 <- 3
            } else {
              pDf <- 1
              pE1 <- 2
              pE2 <- 3
            }
          }
          pointDif <- array(unlist(triPoin[pDf,1:3]))
          pointEq1 <- array(unlist(triPoin[pE1,1:3]))
          pointEq2 <- array(unlist(triPoin[pE2,1:3]))

          # Get the intermediate points between the different point and the other
          # two points. These points are taken as the line dividing the triangle
          # into both colors.
          difColor <- triPoin[pDf,4]
          med1 <- getBorderPoint(pointDif,pointEq1,difColor,gammas,minD/4)
          med2 <- getBorderPoint(pointDif,pointEq2,difColor,gammas,minD/4)

          # The different vertex is too close to a point with the same color
          # of the two other vertices. We finish by setting the color
          # of the equal vertices to the triangle
          if(identical(pointDif,med1) | identical(pointDif,med2)){
            theTriang[curntTriangle,4] <- triPoin[pE1,4]
          } else {
            # In other case the triangle is divided

            # The current triangle is updated and the color is set
            theTriang[curntTriangle,1:3] <- c(triVert[pDf],lastPoint+1,lastPoint+2)
            theTriang[curntTriangle,4] <- difColor

            # Two new triangles are created and two new points
            # First one:
            thePoints[lastPoint+1,1:3] <- med1 # Coordinates
            thePoints[lastPoint+1,  4] <- difColor
            # Second
            thePoints[lastPoint+2,1:3] <- med2
            thePoints[lastPoint+2,  4] <- difColor

            theTriang[totalTriangles+1,1:3] <- c(lastPoint+1, triVert[pE1], triVert[pE2])
            theTriang[totalTriangles+2,1:3] <- c(lastPoint+1, lastPoint+2, triVert[pE2])
            # Set the color:
            theTriang[(totalTriangles+1):(totalTriangles+2),4] <- triPoin[pE1,4]

            # There are 2 more points and 2 more triangles
            lastPoint <- lastPoint+2
            totalTriangles <- totalTriangles+2
          }
        } else {
          # Three different colors, in this case the color of the center
          # is given to the triangle.
          # Middle point
          mp <- (colSums(triPoin[1:3,1:3]))/3
          # IS of the middle point
          ISmp <- ISbuildThird(
            as.data.frame(
              matrix(round(mp,5),1,3)),gammas)
          # Get IS type (color)
          ISstr <- ISschemeToString(ISmp)
          # Set the IS type
          theTriang[curntTriangle,4] <- ISstr
        }
      } else {
        # In other case, divide the triangle into 4 new ones. This is when the
        # size of the triangle is greater than the given limit.

        # Intermediate points between triangle vertices
        n12 <- (unlist(triPoin[1,1:3])+unlist(triPoin[2,1:3]))/2
        n13 <- (unlist(triPoin[1,1:3])+unlist(triPoin[3,1:3]))/2
        n23 <- (unlist(triPoin[2,1:3])+unlist(triPoin[3,1:3]))/2

        # Add the three new points
        # First one:
        thePoints[lastPoint+1,1:3] <- n12 # Coordinates
        thePoints[lastPoint+1,4] <-  # IS type (color)
          ISschemeToString(
            ISbuildThird(round(thePoints[lastPoint+1,1:3],5),gammas))
        # Second
        thePoints[lastPoint+2,1:3] <- n13
        thePoints[lastPoint+2,4] <-
          ISschemeToString(
            ISbuildThird(round(thePoints[lastPoint+2,1:3],5),gammas))
        # Third
        thePoints[lastPoint+3,1:3] <- n23
        thePoints[lastPoint+3,4] <-
          ISschemeToString(
            ISbuildThird(round(thePoints[lastPoint+3,1:3],5),gammas))

        # New triangles
        # Set triangle vertices
        theTriang[totalTriangles+1,1:3] <- c(triVert[1], lastPoint+1, lastPoint+2)
        theTriang[totalTriangles+2,1:3] <- c(triVert[2], lastPoint+1, lastPoint+3)
        theTriang[totalTriangles+3,1:3] <- c(triVert[3], lastPoint+2, lastPoint+3)
        theTriang[totalTriangles+4,1:3] <- c(lastPoint+1, lastPoint+2, lastPoint+3)
        # There is no color yet
        theTriang[(totalTriangles+1):(totalTriangles+4),4] <- "no"

        # There are 3 more points and 4 more triangles
        lastPoint <- lastPoint+3
        totalTriangles <- totalTriangles+4
      }
    }

    # All points are of the same color (the nice easy case)
    if(length(difCols) == 1){
      # It's the color of the triangle
      theTriang[curntTriangle,4] <- difCols[1]
    }

    # Go to the following triangle
    curntTriangle <- curntTriangle+1
  }
  # Return points and triangles
  list(points = thePoints, triang = theTriang)
}

