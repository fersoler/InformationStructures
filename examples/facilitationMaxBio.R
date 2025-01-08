source("R/ISbuild.R")
source("R/ISgraph.R")
source("R/figs/circleCones.R")
source("R/figs/sphereCones.R")
source("R/figs/sphereSubCones.R")
source("R/figs/InformationField.R")




######################################################
#  BIODIVERSITY SUB-CONES FOR 3 SPECIES WITH DIFFERENT
#  INTERACTION MATRICES
######################################################

# Visualization
planeColor <- "#606060"  # Color of the planes
alphaPlanes <- 0.6       # Transparency of the planes
pPlane <- t(matrix(c(1,1,0,1,-1,0,-1,1,0,-1,-1,0), 3, 4))

pPlane <- t(matrix(c(1,1,0,1,0,0,0,1,0,0,0,0), 3, 4))
tPlane <- t(matrix(c(1,2,3,2,3,4), 3, 2))-1


# # Cooperation index
# cI <- 0.1
# (g3 <- t(matrix(c(-1,cI,cI,cI,-1,cI,cI,cI,-1),3,3)))
# sph1 <- sphereSubCones(g3,      # Gamma matrix
#                        c(1,1,1), # Cone
#                        120,       # Initial division on triangles
#                        0.003)     # Stop dividing triangles with vertices closer
# cI <- 0.2
# (g3 <- t(matrix(c(-1,cI,cI,cI,-1,cI,cI,cI,-1),3,3)))
# sph2 <- sphereSubCones(g3,      # Gamma matrix
#                        c(1,1,1), # Cone
#                        120,       # Initial division on triangles
#                        0.003)     # Stop dividing triangles with vertices closer
# cI <- 0.3
# (g3 <- t(matrix(c(-1,cI,cI,cI,-1,cI,cI,cI,-1),3,3)))
# sph3 <- sphereSubCones(g3,      # Gamma matrix
#                        c(1,1,1), # Cone
#                        120,       # Initial division on triangles
#                        0.003)     # Stop dividing triangles with vertices closer
# cI <- 0.4
# (g3 <- t(matrix(c(-1,cI,cI,cI,-1,cI,cI,cI,-1),3,3)))
# sph4 <- sphereSubCones(g3,      # Gamma matrix
#                        c(1,1,1), # Cone
#                        120,       # Initial division on triangles
#                        0.003)     # Stop dividing triangles with vertices closer

#save.image(file = "examples/conosCoopFig.RData")
load("examples/conosCoopFig.RData")

dP1 <- as.matrix(t(sph1$points[,1:3]))
dP2 <- as.matrix(t(sph2$points[,1:3]))*0.9
dP3 <- as.matrix(t(sph3$points[,1:3]))*0.8
dP4 <- as.matrix(t(sph4$points[,1:3]))*0.8

dataPoints <- cbind(dP1,dP2,dP4)#,dP4)

nZT1 <- as.matrix(sph1$triang[sph1$triang[,4] != "no",1:3])-1
nZT2 <- as.matrix(sph2$triang[sph2$triang[,4] != "no",1:3])-1+ncol(dP1)
#nZT3 <- as.matrix(sph3$triang[sph3$triang[,4] != "no",1:3])-1+ncol(dP1)+ncol(dP2)
nZT4 <- as.matrix(sph4$triang[sph4$triang[,4] != "no",1:3])-1+ncol(dP1)+ncol(dP2)#+ncol(dP3)

nonZeroTriang <- rbind(nZT1, nZT2, nZT4)#, nZT4)

colorCodes <- c(sph1$triang[sph1$triang[,4] != "no",4],
                sph2$triang[sph2$triang[,4] != "no",4],
                sph4$triang[sph4$triang[,4] != "no",4])#,
#                    sph4$triang[sph4$triang[,4] != "no",4])
allCodes <- unique(colorCodes)
(n <- length(allCodes))
palette <- distinctColorPalette(n)
for(c in 1:n){
  colorCodes[colorCodes == allCodes[c]] <- palette[c]
}
# Plot
plot_ly(
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
  ))  %>% 
  layout(
    scene = list(
      xaxis = list(range = c(-1, 1), 
                   title = list(text = "r\u2081", 
                                font = list(size = 18, color = "black", 
                                            family = "Arial, sans-serif", 
                                            weight = "bold"))),
      yaxis = list(range = c(-1, 1), 
                   title = list(text = "r\u2082", 
                                font = list(size = 18, color = "black", 
                                            family = "Arial, sans-serif", 
                                            weight = "bold"))),
      zaxis = list(range = c(-1, 1), 
                   title = list(text = "r\u2083", 
                                font = list(size = 18, color = "black", 
                                            family = "Arial, sans-serif", 
                                            weight = "bold")))
    )) %>% add_trace(
      type = "mesh3d",
      x = pPlane[,1], y = pPlane[,2], z = pPlane[,3],
      i = tPlane[,1], j = tPlane[,2], k = tPlane[,3],
      facecolor = toRGB(rep(planeColor, 2), alpha=alphaPlanes)
    ) %>% add_trace(
      type = "mesh3d",
      x = pPlane[,1], y = pPlane[,3], z = pPlane[,2],
      i = tPlane[,1], j = tPlane[,2], k = tPlane[,3],
      facecolor = toRGB(rep(planeColor, 2), alpha=alphaPlanes)
    ) %>% add_trace(
      type = "mesh3d",
      x = pPlane[,3], y = pPlane[,2], z = pPlane[,1],
      i = tPlane[,1], j = tPlane[,2], k = tPlane[,3],
      facecolor = toRGB(rep(planeColor, 2), alpha=alphaPlanes)
    )
# end fig --------------------------------------------

#######
#####   VARIATION
################################


dP1 <- as.matrix(t(sph1$points[,1:3]))
dP2 <- as.matrix(t(sph2$points[,1:3]))*0.99
dP3 <- as.matrix(t(sph3$points[,1:3]))*0.98
dP4 <- as.matrix(t(sph4$points[,1:3]))*0.97

dataPoints <- cbind(dP1,dP2,dP3,dP4)

nZT1 <- as.matrix(sph1$triang[sph1$triang[,4] != "no",1:3])-1
nZT2 <- as.matrix(sph2$triang[sph2$triang[,4] != "no",1:3])-1+ncol(dP1)
nZT3 <- as.matrix(sph3$triang[sph3$triang[,4] != "no",1:3])-1+ncol(dP1)+ncol(dP2)
nZT4 <- as.matrix(sph4$triang[sph4$triang[,4] != "no",1:3])-1+ncol(dP1)+ncol(dP2)+ncol(dP3)

nonZeroTriang <- rbind(nZT1, nZT2, nZT3, nZT4)

colorCodes <- rep(c("blue", "green", "red", "orange"), c(nrow(nZT1),nrow(nZT2),nrow(nZT3),nrow(nZT4)))


# Plot
plot_ly(
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
  ))  %>% 
  layout(
    scene = list(
      xaxis = list(range = c(-1, 1), 
                   title = list(text = "r\u2081", 
                                font = list(size = 18, color = "black", 
                                            family = "Arial, sans-serif", 
                                            weight = "bold"))),
      yaxis = list(range = c(-1, 1), 
                   title = list(text = "r\u2082", 
                                font = list(size = 18, color = "black", 
                                            family = "Arial, sans-serif", 
                                            weight = "bold"))),
      zaxis = list(range = c(-1, 1), 
                   title = list(text = "r\u2083", 
                                font = list(size = 18, color = "black", 
                                            family = "Arial, sans-serif", 
                                            weight = "bold")))
    )) %>% add_trace(
      type = "mesh3d",
      x = pPlane[,1], y = pPlane[,2], z = pPlane[,3],
      i = tPlane[,1], j = tPlane[,2], k = tPlane[,3],
      facecolor = toRGB(rep(planeColor, 2), alpha=alphaPlanes)
    ) %>% add_trace(
      type = "mesh3d",
      x = pPlane[,1], y = pPlane[,3], z = pPlane[,2],
      i = tPlane[,1], j = tPlane[,2], k = tPlane[,3],
      facecolor = toRGB(rep(planeColor, 2), alpha=alphaPlanes)
    ) %>% add_trace(
      type = "mesh3d",
      x = pPlane[,3], y = pPlane[,2], z = pPlane[,1],
      i = tPlane[,1], j = tPlane[,2], k = tPlane[,3],
      facecolor = toRGB(rep(planeColor, 2), alpha=alphaPlanes)
    )


#######################################
# MEASURING THE MAX BIODIVERSITY CONE
########################################

library(geosphere)
library(pracma)

# TEST
# Octante positivo:
(p <- rbind(c(0,0), c(0, 90), c(90,0), c(0,0)))
# Proporción de la superficie: 
areaPolygon(p,1,0)/(4*pi)

# Get the longitude and latitude of a point given cartesian coordinates
getLongLat <- function(coord){
  s <- cart2sph(coord)
  s[1:2]*360/(2*pi)
}

# Function to measure the are of the max. biodiversity cone for a system
# of 3 species with the given interaction matrix
areaMaxBioCone <- function(gMatrix){
  n <- 1
  rawMesh <- BuildMeshPattern(n)
  pointsTriang <- setPointsCone(-gMatrix,c(1,1,1),rawMesh$vb)
  polygon <- rbind(getLongLat(pointsTriang[,1]),
                   getLongLat(pointsTriang[,2]),
                   getLongLat(pointsTriang[,3]),
                   getLongLat(pointsTriang[,1]))
  areaPolygon(polygon,1,0)/(4*pi)
}

# Measure the area of the max. bio. cone for a system with a fixed given
# interaction parameter
areaMaxWithConstant <- function(val){
  iC <- val
  g3 <- t(matrix(iC,3,3))
  diag(g3) <- -1
  areaMaxBioCone(g3)
}

# Random matrices
meanVals <- rep(0,900)
sphPorts <- rep(0,900)
for (i in 1:300) {
  m <- matrix(runif(9,-0.45,0.45),3,3)
  diag(m) <- 0
  meanVals[i] <- sum(m)/6
  diag(m) <- -1
  sphPorts[i] <- areaMaxBioCone(m)
}
for (i in 301:600) {
  m <- matrix(runif(9,-0.45,0),3,3)
  diag(m) <- 0
  meanVals[i] <- sum(m)/6
  diag(m) <- -1
  sphPorts[i] <- areaMaxBioCone(m)
}
for (i in 601:900) {
  m <- matrix(runif(9,0,0.45),3,3)
  diag(m) <- 0
  meanVals[i] <- sum(m)/6
  diag(m) <- -1
  sphPorts[i] <- areaMaxBioCone(m)
}

# Constant matrices

coopIndex <- seq(-0.49,0.49,0.01)
spherePortion <- sapply(coopIndex,areaMaxWithConstant)


plot(coopIndex, spherePortion, type = "l", col = "blue", lwd = 2, 
     xlab = "Interaction parameter (mean)", ylab = "Portion of the sphere", 
     main = "Max. biodiversity cone for 3 species")
points(meanVals,sphPorts, col = "red", pch = 16, cex = 1)
lines(coopIndex, spherePortion, col = "blue", lwd = 2)


