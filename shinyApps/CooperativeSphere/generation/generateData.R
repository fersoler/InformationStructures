source("R/ISbuild.R")
source("R/ISgraph.R")
source("R/figs/circleCones.R")
source("R/figs/sphereCones.R")
source("R/figs/InformationField.R")


# Gamma values
(g3 <- t(matrix(c(-1,.3,.3,.3,-1,.3,.3,.3,-1),3,3)))

## Generate all cones to form a sphere
sph111 <- sphereSubCones(g3,
                         c(1,1,1),
                         80,  
                         0.01)
sph110 <- sphereSubCones(g3,
                         c(1,1,0),
                         80,  
                         0.01)
sph101 <- sphereSubCones(g3,
                         c(1,0,1),
                         80,  
                         0.01)
sph011 <- sphereSubCones(g3,
                         c(0,1,1),
                         80,  
                         0.01)
sph100 <- sphereSubCones(g3,
                         c(1,0,0),
                         80,  
                         0.01)
sph010 <- sphereSubCones(g3,
                         c(0,1,0),
                         80,  
                         0.01)
sph001 <- sphereSubCones(g3,
                         c(0,0,1),
                         80,  
                         0.01)
sph000 <- sphereSubCones(g3,
                         c(0,0,0),
                         80,  
                         0.01)

# All spheres
listSpheres <- list(sph111, sph110, sph101, sph011, sph100, sph010, sph001, sph000)

# To reduce the radius of some cones
fM <- .95
multVals <- c(1,fM,fM,fM,1,1,1,fM)
# Points: 
allPoints <- data.frame(p1=numeric(), p2=numeric(), p3=numeric())
for(n in 1:8){
  esfera <- listSpheres[[n]]
  thesePoints <- round(esfera$points[,1:3]*multVals[n],4)
  allPoints <- rbind(allPoints,thesePoints)
}
dataPoints <- as.matrix(t(allPoints))


# Generate triangles: 
prevPoints <- 0
allTriangles <- data.frame(v1=numeric(), v2=numeric(), v3=numeric(), color=character())
for(n in 1:8){
  esfera <- listSpheres[[n]]
  theseTriangles <- esfera$triang
  theseTriangles[,1:3] <- theseTriangles[,1:3]+prevPoints-1
  allTriangles <- rbind(allTriangles,theseTriangles)
  prevPoints <- prevPoints + nrow(esfera$points)
}

nonZeroTriang <- as.matrix(allTriangles[allTriangles[,4] != "no",1:3])

# The lines below are to convert codes representing ISs into random colors
# cange for other options
colorCodes <- allTriangles[allTriangles[,4] != "no",4]
allCodes <- unique(colorCodes)
n <- length(allCodes)
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
  lighting = list(ambient = 0.5,
                  diffuse = 0.8,
                  fresnel = 0.1,
                  specular = 0.5,
                  roughness = 0.5,
                  facenormalsepsilon = 0,
                  vertexnormalsepsilon = 0))

## Output files: 
write.csv(dataPoints, "pointsCoopSphere.csv", row.names = FALSE)
write.csv(nonZeroTriang, "trianglesCoopSphere.csv", row.names = FALSE)
write.csv(colorCodes, "colorsCoopSphere.csv", row.names = FALSE, quote = TRUE)



# Finally, create and save an envoronment with only these files:
spherePoints <- as.matrix(read.csv("pointsCoopSphere.csv"))
sphereTriang <- as.matrix(read.csv("trianglesCoopSphere.csv"))
sphereColors <- as.matrix(read.csv("colorsCoopSphere.csv", stringsAsFactors = FALSE))


