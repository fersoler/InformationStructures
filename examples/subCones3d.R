#####################################################################
# Example of conversion of sub-cones surfaces into 3D volumes and   #
# st files for 3d print                                             #
#####################################################################

library(plotly)
library(pracma)
library(tidyverse)
library(shape)
library(dplyr)
library(randomcoloR)
library(RColorBrewer)
library(deSolve)
library(rgl)
source("R/ISbuild.R")
source("R/figs/InformationField.R")
source("R/figs/sphereCones.R")
source("R/figs/sphereSubCones.R")
source("examples/auxFiles/subCones3d_aux.R")


# This file begins by building the subcones representing the zones in the
# feasibility domain of a 3-species LV system. 
# 
# Subcones are represented by a mesh of points and triangles. So the next step
# will be to find the points and lines separating the subcones. Points will be set
# by finding the intersecion of geodesics fitted from the points in the mesh 
# belonging to triangles of different colors. 

# Matrix of interactions: 
(g3 <- t(matrix(c(-1,.25,.3,.3,-1,.25,.34,.32,-1),3,3)))

# Relevant points of the resulting cone.
# These points will be used to correct the intersections of geodesics. 
# Relevant points are the limits of the cone and limits of the
# positive octant. 
relevantPoints <- rbind(t(apply(-g3, 2, normalize)), diag(3))

# Now be build the mesh of points and triangles. 

# Initial division of the cone into triangles.
divT      <- 60       # Quick with poorer quality
#divT <- 120             # Slower but better quality
#divT <- 150
# Triangles will stop dividing when the distance among vertices is smaller than
# 'subDivTol': 
#subDivTol <- 0.01     # Quick with poorer quality
#subDivTol <- 0.003      # Slower but better quality
subDivTol <- 0.001
# Build the mesh: 
sph <- sphereSubCones(g3, c(1,1,1),divT,subDivTol)

# Get the points (vertices of the triangles)
dataPoints <- as.matrix(t(sph$points[,1:3]))
# Triangles to be plotted are those with a color
nonZeroTriang <- as.matrix(sph$triang[sph$triang[,4] != "no",1:3])

## Remove duplicate points
## This is a necessary step to avoid duplication of points with
## different indices
points_t <- t(dataPoints)
unique_points <- unique(points_t)
unique_indices <- match(data.frame(t(points_t)), data.frame(t(unique_points)))
dataPoints_unique <- t(unique_points)
Triangles_updated <- matrix(unique_indices[nonZeroTriang], ncol = 3)
dataPoints <- dataPoints_unique
nonZeroTriang <- Triangles_updated-1

## We can display the initial mesh of points and triangles. 
## This is what will be later improved. 
colorCodes <- sph$triang[sph$triang[,4] != "no",4]
allCodes <- unique(colorCodes)
(n <- length(allCodes))
palette <- distinctColorPalette(n)
for(c in 1:n){
  colorCodes[colorCodes == allCodes[c]] <- palette[c]
}
# Plot
fig <- plot_ly(
  x = array(dataPoints[1,]),
  y = array(dataPoints[2,]),
  z = array(dataPoints[3,]),
  i = array(nonZeroTriang[,1]),
  j = array(nonZeroTriang[,2]),
  k = array(nonZeroTriang[,3]),
  facecolor = toRGB(colorCodes,alpha=1),
  type = "mesh3d",
  flatshading = TRUE,
  lighting = list(ambient = 0.6,
                  diffuse = 0.8,
                  fresnel = 0.1,
                  specular = 0.5,
                  roughness = 0.5,
                  facenormalsepsilon = 0,
                  vertexnormalsepsilon = 0)) %>%
  layout(scene = list(
    xaxis = list(visible=FALSE),
    yaxis = list(visible=FALSE),
    zaxis = list(visible=FALSE)
  ))
fig

##################################################################
# Set the intersections of subcones
##################################################################

# Colos codes of each colored triangle
colorCodes2 <- sph$triang[sph$triang[,4] != "no",4]

# Remove triangles with less than 10 occurrences
# These can be errors due to parameters in the limits of
# subcones producing incorrect IS. 
cleanData <- cleanNonSignificantZones(nonZeroTriang, colorCodes2, 10)
colorCodes2 <- cleanData$colorCodes
nonZeroTriang <- cleanData$triangles

# Different colors
allCodes2 <- unique(colorCodes2)
# Number of colors
(nCols <- length(allCodes2))
# Associate each color code with a number
for(c in 1:nCols){
  colorCodes2[colorCodes2 == allCodes2[c]] <- c
}
colorCodes2 <- as.numeric(colorCodes2)

# Number of points and triangles
nPoints <- dim(dataPoints)[2]
nTriang <- dim(nonZeroTriang)[1]

# Matrix to record to which color is associated each point. A point can be the
# vertex of more than one triangle, with different colors. 
# Each row correspond to a points, columns are:
# - Column 1: Number of the point
# - Column 2: 1 is the coordinates have been changed, 0 otherwise
# - Cols.  3 to 5: Coordinates of the point (x, y, z)
# - Column 6: Indicates the number of triangles in which the point appears
# - Cols.  7-end: Each column is one color, the value is 1 if the point appears in
# at least a triangle of the corresponding color and 0 otherwise. 
PointsColors <- matrix(0, ncol = 6+nCols, nrow = nPoints)
for(nt in 1:nTriang){
  if(colorCodes2[nt] > 0){
    PointsColors[nonZeroTriang[nt,]+1, 6+colorCodes2[nt]] <- 1+PointsColors[nonZeroTriang[nt,]+1, 6+colorCodes2[nt]]
    PointsColors[nonZeroTriang[nt,]+1, 6] <- 1+PointsColors[nonZeroTriang[nt,]+1, 6]
  }
}
PointsColors[,1] <- c(1:nPoints)
PointsColors[,2] <- rep(0, nPoints)
PointsColors[,3] <- dataPoints[1,]
PointsColors[,4] <- dataPoints[2,]
PointsColors[,5] <- dataPoints[3,]

# Border points, points in the border of the great cone
borderPoints <- PointsColors[which(PointsColors[,6] < 6 ),]
borderPoints <- borderPoints[1:(3*divT),]

# Table to record found limits between regions
# Column 1: inner color
# Column 2: outer color (or 0 for 'border' of external regions)
# Cols 3-5: coordinates of the starting points
# Cols 6-8: coordinates of the end point
# Cols 9-11: normal vector
foundLimitsGeneral <- NULL

for(subCone in 1:nCols){
  
  print(subCone)  
  # Setting the border of a subcone
  # For the single subcone
  foundLimits <- NULL
  
  # Select all points in triangles with color 'subCone' and another one. These 
  # points are the borders of the cone 
  selectedPoints <- PointsColors[which(PointsColors[,6+subCone] > 0 & 
                                         rowSums(PointsColors[,7:(6+nCols)]) >= 2),]
  
  # Adjacent cones
  adjacentColors <- which(colSums(selectedPoints[,7:(6+nCols)]) > 0)
  adjacentColors <- adjacentColors[adjacentColors != subCone]
  # Points of the border of the cone belonging to the sub-cone
  selectedBorder <- borderPoints[which(borderPoints[,6+subCone] > 0),]
  # Look at the adjacent cones and draw the borders
  for (adjCol in adjacentColors){
    
    # Points of the corresponding border
    pointsToInterpolate <- selectedPoints[which(selectedPoints[,c(adjCol+6)] > 0),]
    
    if(adjCol < subCone & length(pointsToInterpolate) >= 10*(nCols+6)){
      findInGeneral <- which(foundLimitsGeneral[, 1] == adjCol & foundLimitsGeneral[, 2] == subCone)
      foundLimits <- rbind(foundLimits, c(subCone,adjCol, 
                                          foundLimitsGeneral[findInGeneral, 3:11]))
    } 
    
    # Require at least 10 points (x3 coordinates) to avoid corner borders
    if(adjCol > subCone & length(pointsToInterpolate) >= 10*(nCols+6)){
      
      # Adjust the points
      adjusted_points_normal <- adjustPoints2(pointsToInterpolate[,3:5])
      adjusted_points <- adjusted_points_normal$points
      adjusted_normal <- adjusted_points_normal$normal
      
      foundLimits <- rbind(foundLimits, c(subCone,adjCol, 
                                          adjusted_points[1,], adjusted_points[2,],adjusted_normal))
    }
  }
  if(length(selectedBorder) > 4*(nCols+6)){
    pointsToInterpolate <- selectedBorder[,3:5]
    adjusted_points_normal <- adjustPoints2(pointsToInterpolate)
    adjusted_points <- adjusted_points_normal$points
    adjusted_normal <- adjusted_points_normal$normal
    foundLimits <- rbind(foundLimits, c(subCone,0, 
                                        adjusted_points[1,], adjusted_points[2,],adjusted_normal))
  }
  
  # Order and adjust the points with the found limits
  foundLimits_order <- sort_edges(foundLimits)
  foundLimits_adjust <- adjustIntersections(foundLimits_order)
  
  # Store the found limits
  foundLimitsGeneral <- rbind(foundLimitsGeneral, foundLimits_adjust)
  
}

foundLimitsGeneral <- adjustWithRelevantPoints(foundLimitsGeneral, relevantPoints, 0.01)

## Display limits of cones
fig2 <- fig
fig2 <- fig2 %>% add_trace(x = foundLimitsGeneral[,3],
                           y = foundLimitsGeneral[,4],
                           z = foundLimitsGeneral[,5],
                           type = 'scatter3d', mode = 'markers',
                           showlegend = FALSE,
                           name = 'IS points',
                           hovertemplate = paste('<b>sp1</b>: %{x:.2f}',
                                                 '<br><b>sp2</b>: %{y:.2f}',
                                                 '<br><b>sp3</b>: %{z:.2f}'),
                           marker = list(size = 3, color = "yellow", opacity = 1, symbol = 'circle'))

fig2 <- fig2 %>% add_trace(x = foundLimitsGeneral[,6],
                           y = foundLimitsGeneral[,7],
                           z = foundLimitsGeneral[,8],
                           type = 'scatter3d', mode = 'markers',
                           showlegend = FALSE,
                           name = 'IS points',
                           hovertemplate = paste('<b>sp1</b>: %{x:.2f}',
                                                 '<br><b>sp2</b>: %{y:.2f}',
                                                 '<br><b>sp3</b>: %{z:.2f}'),
                           marker = list(size = 3, color = "yellow", opacity = 1, symbol = 'circle'))

fig2

# Read/write
write.csv(foundLimitsGeneral, "examples/auxFiles/foundLimitsGeneral.csv", row.names = FALSE)
#foundLimitsGeneral <- as.matrix(read.csv("examples/auxFiles/foundLimitsGeneral.csv"))


#########################################################################
# EXTRACTION OF INDIVIDUAL CONES AND SPHERE
#########################################################################

# Figure to control the quality of the output
figC <- plot_ly() %>%
  layout(scene = list(
    xaxis = list(visible=FALSE),
    yaxis = list(visible=FALSE),
    zaxis = list(visible=FALSE)
  ))

### Individual subcones

for (nC in 1:nCols){
  print(nC)
  # Select limits
  selRows <- which(foundLimitsGeneral[,1] == nC)
  # Create subcone
  newMesh <- meshSubConeFromPoints(foundLimitsGeneral[selRows,3:5], 1, .5, 120, top = TRUE)
  # Extract points and triangles
  dataPoints <- newMesh$points
  Triang <- newMesh$triangles
  
  # Remove duplicate points
  points_t <- dataPoints
  unique_points <- unique(points_t)
  unique_indices <- match(data.frame(t(points_t)), 
                          data.frame(t(unique_points)))
  Triangles_updated <- matrix(unique_indices[Triang], ncol = 3)
  dataPoints <- unique_points
  Triang <- Triangles_updated
  
  # Remove points not appearing in triangles
  usedPoints <- unique(as.vector(Triang))
  Triang <- matrix(match(Triang,usedPoints), ncol = 3)
  dataPoints <- dataPoints[usedPoints,]
  
  # Create 3d mesh
  mesh <- tmesh3d(
    vertices = t(dataPoints),
    indices = t(Triang)
  )
  
  figC <- figC %>% add_trace(
    x = array(dataPoints[,1]),
    y = array(dataPoints[,2]),
    z = array(dataPoints[,3]),
    i = array(Triang[,1]-1),
    j = array(Triang[,2]-1),
    k = array(Triang[,3]-1),
    facecolor = rep(palette[nC],nrow(Triang)),
    type = "mesh3d",
    flatshading = TRUE,
    lighting = list(ambient = 0.6,
                    diffuse = 0.8,
                    fresnel = 0.1,
                    specular = 0.5,
                    roughness = 0.5,
                    facenormalsepsilon = 0,
                    vertexnormalsepsilon = 0))
  
  # Save stl file
  open3d()
  shade3d(mesh)
  writeSTL(paste("examples/stl_files/subcone",nC,".stl", sep=""))
  
  
}
figC

#########
# Rest of the sphere
##########

# Center part
newMesh <- meshSubConeFromPoints(normalize(relevantPoints[1:3,1:3]), 1, .5, 150, top = FALSE)
dataPoints <- newMesh$points
Triang <- newMesh$triangles 

# Sides
n <- 150
# Raw mesh that will be instantiated in each cone
rawMesh <- BuildMeshPattern(n)
# Information of triangles:
it <- rawMesh$it 

otherSides <- rbind(diag(3), 1-diag(3), rep(0,3))
for(sn in 1:nrow(otherSides)){
  newPoints <- t(setPointsCone(-g3,otherSides[sn,],rawMesh$vb))
  Triang <- rbind(Triang, it+nrow(dataPoints))
  dataPoints <- rbind(dataPoints, newPoints)
}

## Remove duplicate points
points_t <- dataPoints
unique_points <- unique(points_t)
unique_indices <- match(data.frame(t(points_t)), 
                        data.frame(t(unique_points)))
Triangles_updated <- matrix(unique_indices[Triang], ncol = 3)
dataPoints <- unique_points
Triang <- Triangles_updated

## Remove points not appearing in triangles
usedPoints <- unique(as.vector(Triang))
Triang <- matrix(match(Triang,usedPoints), ncol = 3)
dataPoints <- dataPoints[usedPoints,]

# Make a base
convertedPoints <- dataPoints
radMod <- 0.05  # modified sphere radius
angMax <- acos(radMod)
radCut <- sin(angMax)
centerCirc <- normalize(c(-1,-1,-1))*radMod
for(np in 1:nrow(dataPoints)){
  ang <- angle_between_vectors(dataPoints[np,], c(-1,-1,-1))
  disTop <- distance(dataPoints[np,], normalize(c(1,1,1)))
  disBot <- distance(dataPoints[np,], normalize(c(-1,-1,-1)))
  if((abs(ang) < angMax) & (disTop > disBot)){
    convertedPoints[np,] <- dataPoints[np,]*(radMod/cos(ang))
  }
}

for(nt in 1:nrow(Triang)){
  pointNs <- Triang[nt,]
  modSome <- 0
  modAll <- 1
  for(np in pointNs){
    if(distance(convertedPoints[np,], 
                dataPoints[np,]) > 0){
      modSome <- 1
    } else {
      modAll <- 0
    }
  }
  if(modSome == 1 & modAll == 0){
    for(np in pointNs){
      if(distance(convertedPoints[np,], dataPoints[np,]) > 0){
        convertedPoints[np, ] <- 
          centerCirc + normalize(convertedPoints[np, ] - centerCirc)*radCut
      }
    } 
  }
}


# View
figC <- figC %>% add_trace(
  x = array(convertedPoints[,1]),
  y = array(convertedPoints[,2]),
  z = array(convertedPoints[,3]),
  i = array(Triang[,1]-1),
  j = array(Triang[,2]-1),
  k = array(Triang[,3]-1),
  facecolor = rep("lightgrey",nrow(Triang)),
  type = "mesh3d",
  flatshading = TRUE,
  lighting = list(ambient = 0.6,
                  diffuse = 0.8,
                  fresnel = 0.1,
                  specular = 0.5,
                  roughness = 0.5,
                  facenormalsepsilon = 0,
                  vertexnormalsepsilon = 0)) %>%
  layout(scene = list(
    xaxis = list(visible=FALSE),
    yaxis = list(visible=FALSE),
    zaxis = list(visible=FALSE)
  ))
figC

## 3D print
mesh <- tmesh3d(
  vertices = t(convertedPoints),
  indices = t(Triang)
)

open3d()
shade3d(mesh)
writeSTL("examples/stl_files/baseAl.stl")
