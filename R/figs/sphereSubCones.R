#################################################################
# Functions to draw the biodiversity sub-cones 
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
    ISp <- ISbuild(as.matrix(round(as.numeric(thePoints[p,1:3]),5),nrow=1),gammas)
    # Get the GASS
    gass <- ISp$points[ISp$gassInd,]
    # Check is the GASS is in the intended cone.
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
    while (theTriang[curntTriangle,4] != "no"){
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
        # Three different colors, in this case the color of the center
        # is given to the triangle.
        # Middle point
        mp <- (colSums(triPoin[1:3,1:3]))/3
        # IS of the middle point
        ISmp <- ISbuild(as.matrix(round(mp,5),nrow=1),gammas)
        # Get IS type (color)
        ISstr <- ISschemeToString(ISmp)
        # Set the IS type
        theTriang[curntTriangle,4] <- ISstr
        #}
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
            ISbuild(as.matrix(round(thePoints[lastPoint+1,1:3],5),nrow=1),gammas))
        # Second
        thePoints[lastPoint+2,1:3] <- n13
        thePoints[lastPoint+2,4] <-
          ISschemeToString(
            ISbuild(as.matrix(round(thePoints[lastPoint+2,1:3],5),nrow=1),gammas))
        # Third
        thePoints[lastPoint+3,1:3] <- n23
        thePoints[lastPoint+3,4] <-
          ISschemeToString(
            ISbuild(as.matrix(round(thePoints[lastPoint+3,1:3],5),nrow=1),gammas))
        
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

