source("R/ISbuildImproved.R")
source("R/ISgraph.R")
source("figures/circleCones.R")
source("figures/sphereCones.R")
source("figures/InformationField.R")

#####
# GENERAL FUNCTIONS

ISgraphDrawLabels2 <- function(IS, ISgr, ISlay, colorEdges){
  plot(ISgr$graph,
       vertex.label = ISgr$vlab,
       layout = ISlay/2.5,
       edge.width=4,
       vertex.shape = "circle",
       vertex.color = "white",
       vertex.label.cex = 1,
       edge.color = colorEdges,
       mark.groups=list(c(IS$gassInd)),
       mark.col="orange",
       rescale = FALSE,
       mark.expand=4,
       mark.shape=-2)
}


showIScolor <- function(aVals, gMatrix, colorEdges){
  IS3 <- ISbuildThird(as.data.frame(t(aVals)),gMatrix)
  gr3 <- ISgraph(IS3, 1:3)
  ISgraphDrawLabels2(IS3,gr3, ISgraphLayout(IS3, gr3, "3Dcube"), colorEdges)
}

#showIScolor(spherePoints[,3000], g3, "blue")

# Look at the first triangle with a given color and draw the
# corresponding IS
drawIsColorNumber <- function(n){
  # Index of the first triangle of the given color
  triangNumber <- which(colorCodes == palette[n])[1]
  # Vector to the center of the triangle
  point <- rowSums(dataPoints[,nonZeroTriang[triangNumber,]])/3
  a <- as.vector(point)
  showIScolor(a, g3, palette[n])
}

#drawIsColorNumber(1)

#################################################################
# Sphere and triangle surfaces
#################################################################

surface_3d_triangle <- function(v1, v2, v3) {
  l1 <- v2 - v1
  l2 <- v3 - v1
  crossed_prod <- cross(l1, l2)
  sqrt(sum(crossed_prod^2)) / 2 
}

# crossed product
cross <- function(v1, v2) {
  if (length(v1) != 3 || length(v2) != 3) {
    stop("incorrect values")
  }
  
  result <- numeric(3)
  
  result[1] <- v1[2] * v2[3] - v1[3] * v2[2]
  result[2] <- v1[3] * v2[1] - v1[1] * v2[3]
  result[3] <- v1[1] * v2[2] - v1[2] * v2[1]
  
  result
}


# Ejemplo de uso de la función
# v1 <- c(1,10,0)
# v2 <- c(-3,-7,-3)
# v3 <- c(-2, -2, 5)
# 
# surface_3d_triangle(v1, v2, v3)


##########################################
#  COOPERATION
##########################################

(g3 <- t(matrix(c(-1,.3,.3,.3,-1,.3,.3,.3,-1),3,3)))

# Slower but very good quality:
sph <- sphereSubCones(g3, c(1,1,1),120,0.003)

# Quick (several minutes) with poorer quality
sph <- sphereSubCones(g3,       # Gamma matrix
                     c(1,1,1), # Cone
                     80,       # Initial division on triangles
                     0.01)     # Stop dividing triangles with vertices closer
# 
# # Get the points (vertices of the triangles)
dataPoints <- as.matrix(t(sph$points[,1:3]))
# # Triangles to be plotted are those with a color
nonZeroTriang <- as.matrix(sph$triang[sph$triang[,4] != "no",1:3])-1
# # The lines below are to convert codes representing ISs into random colors
# # cange for other options
colorCodes <- sph$triang[sph$triang[,4] != "no",4]
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
  ))


n
drawIsColorNumber(2)


##########################################
#  COOPERATION AND COMPETITION
##########################################


(g3 <- t(matrix(c(-1,.3,-.2,.15,-1,.23,-.31,.2,-1),3,3)))

# Slower but very good quality:
sph <- sphereSubCones(g3, c(1,1,1),120,0.005)

# Quick (several minutes) with poorer quality
#sph <- sphereSubCones(g3,       # Gamma matrix
#                     c(1,1,1), # Cone
#                     80,       # Initial division on triangles
#                     0.01)     # Stop dividing triangles with vertices closer
# 
# # Get the points (vertices of the triangles)
dataPoints <- as.matrix(t(sph$points[,1:3]))
# # Triangles to be plotted are those with a color
nonZeroTriang <- as.matrix(sph$triang[sph$triang[,4] != "no",1:3])-1
# # The lines below are to convert codes representing ISs into random colors
# # cange for other options
colorCodes <- sph$triang[sph$triang[,4] != "no",4]
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

n

# 5

drawIsColorNumber(2)

c <- 17
aVals <- c(0.299, 0.396, 0.068)
showIScolor(aVals, g3, palette[c])

# "#D9DFA1" "#D4C669" "#7BE744" "#5EDE75" "#DFB9D2" "#D39175" "#CDE3EA"
# [8] "#EAAA36" "#DD6251" "#DE5ED3" "#74E8E3" "#9AAFE5" "#816F8B" "#DC9BE2"
# [15] "#799E8A" "#6CE6B0" "#E2D1B7" "#DA6293" "#9E3DDF" "#A5DC7B" "#CEE046"
# [22] "#B5E7C5" "#6ABDDB" "#7F75D7"



##########################################
#  COMPETITION
##########################################


(g3 <- t(matrix(c(-1,-0.1,-0.23,-0.3,-1,-0.12,-0.2,-0.17,-1),3,3)))

# Slower but very good quality:
sph <- sphereSubCones(g3, c(1,1,1),120,0.005)

# Quick (several minutes) with poorer quality
#sph <- sphereSubCones(g3,       # Gamma matrix
#                     c(1,1,1), # Cone
#                     80,       # Initial division on triangles
#                     0.01)     # Stop dividing triangles with vertices closer
# 
# # Get the points (vertices of the triangles)
dataPoints <- as.matrix(t(sph$points[,1:3]))
# # Triangles to be plotted are those with a color
nonZeroTriang <- as.matrix(sph$triang[sph$triang[,4] != "no",1:3])-1
# # The lines below are to convert codes representing ISs into random colors
# # change for other options
colorCodes <- sph$triang[sph$triang[,4] != "no",4]
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

n

# 5
palette
drawIsColorNumber(1)

which(colorCodes == palette[3])


