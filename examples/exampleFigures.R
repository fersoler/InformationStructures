source("R/ISbuild.R")
source("R/ISgraph.R")
source("R/figs/circleCones.R")
source("R/figs/sphereCones.R")
source("R/figs/sphereSubCones.R")
source("R/figs/InformationField.R")


##########################################
#  INFORMATION FIELDS (2 species)
##########################################

# Alpha values
a <- c(2,1.8)
# Gamma matrix
g <- matrix(c(-1,0.2,0.3,-1),2,2)

# Build the IS
IS <- ISbuild(a,g)
gr <- ISgraph(IS, 1:2)
# Draw the IS in two different ways:
ISgraphDrawLabels(IS,gr, ISgraphLayout(IS, gr, "tree"))
ISgraphDrawPie(IS, gr, ISgraphLayout(IS, gr, "tree"), c("brown", "purple"))
# Plot the Lyapunov function
plotLyapFunc(IS, 120)



## PLOT THE DIRECTIONS IN THE PHASE SPACE
## The code below uses the same values for alpha (a) and gamma (g) defined
## above. So it's another representation of the IF as a vector space.

# begin fig ------------------------------------------
# Create the canvas with appropriate limits
limX <- max(IS$points[,1])*1.3 # Limit for x
limY <- max(IS$points[,2])*1.3 # Limit for y
plot(0, 0, xlim = c(0, limX), ylim = c(0, limY),
     type = 'n',xlab = "u1", ylab = "u2")
# Plot vectors showing the direction of each point
for(x in seq(0,limX,limX/15)){
  for(y in seq(0,limY,limY/15)){
    dir <- getDirection(a,g,c(x,y),.01)
    Arrows(
      x0 = x, y0 = y,
      x1 = dir[1], y1 = dir[2],
      arr.length = 0.2, arr.type = 'triangle'
    )
  }
}
## Show the trayectories for some initial points
initPoints <- data.frame(y1 = c(0.01, 0.2,0.01,3, 0.3),  # u1 values
                         y2 = c(0.5, 0.01,2.3,0.25, 0.3),   # u2 values
                         cols = c("green", "orange", "blue","brown", "red")) # Colors

# Solve gLV equations and draw trajectories
pars <- c(a1 = a[1], a2 = a[2], g12 = g[1,2], g21 = g[2,1])
times <- seq(0,1000,length.out = 2)
for(p in 1:nrow(initPoints)){
  y0 <- c(y1 = initPoints[p,1], y2 = initPoints[p,2])
  sol2 <- ode(y = y0, times = seq(0,20,length.out = 10000),
              parms = pars, func = LotkaVolterra, atol = 1e-12)
  drawLineArrow(sol2, initPoints[p,3], 1)
  points(sol2[1,'y1'], sol2[1,'y2'], pch = 19, col = initPoints[p,3], cex=1.5)
}
# Finally, plot the IS nodes
points(IS$points[,1], IS$points[,2], cex = 2, pch = 19, col = "black")
# end fig -------------------------------------------

# Look at the Lyapunov function with the trajectories of the 
# given starting points:
plotLyapFuncWithInitPoints(IS, g, a, initPoints, 10)

##########################################
#  BIODIVERSITY CONES FOR 2 SPECIES
##########################################

# With the same values above:
drawCircleCones(g)

# Cooperation:
(gEx = t(matrix(c(-1, 0.5, 0.4, -1), ncol=2)))
drawCircleCones(gEx)
# Other ways of drawing the 2D cones: 
# 1. Draw all cones without labels:  
drawCircleCones(gEx, allCones = TRUE, drawLabels = FALSE)
# 2. Remove cones out of 11:  
drawCircleCones(gEx, allCones = FALSE, drawLabels = FALSE)
# 3. Add labels just for the 11 cone: 
drawCircleCones(gEx, allCones = FALSE, drawLabels = TRUE)

# Competition:
(gEx = t(matrix(c(-1, -0.6, -0.5, -1), ncol=2)))
drawCircleCones(gEx)

# Cooperation (2->1) and competition (1->2):
(gEx = t(matrix(c(-1, 0.45, -0.38, -1), ncol=2)))
drawCircleCones(gEx)


##########################################
#  BIODIVERSITY CONES FOR 3 SPECIES
##########################################

# Gamma values
(g3 <- t(matrix(c(-1,.3,.3,.3,-1,.3,.3,.3,-1),3,3)))
# Alpha values
(a3 <- c(1,-.2,1.3))
# Information Structure:
IS3 <- ISbuild(a3,g3)
gr3 <- ISgraph(IS3, 1:3)
# Draw the IS in different ways:
ISgraphDrawLabels(IS3,gr3, ISgraphLayout(IS3, gr3, "tree"))
ISgraphDrawLabels(IS3,gr3, ISgraphLayout(IS3, gr3, "3Dcube"))
ISgraphDrawPie(IS3,gr3, ISgraphLayout(IS3, gr3, "tree"), c("#806000", "#002080", "#408000"))
ISgraphDrawPie(IS3,gr3, ISgraphLayout(IS3, gr3, "3Dcube"), c("#806000", "#002080", "#408000"))

## Draw the sphere with some cones:

# begin fig --------------------------------------------
# We set n to divide each cone into n^2 triangles
# The higher value, the better visualization
n <- 200
# Raw mesh that will be instantiated in each cone
rawMesh <- BuildMeshPattern(n)
# Information of triangles:
it <- rawMesh$it-1

# Instantiate the points for several cones:
vb111 <- setPointsCone(-g3,c(1,1,1),rawMesh$vb)
vb100 <- setPointsCone(-g3,c(1,0,0),rawMesh$vb)
vb010 <- setPointsCone(-g3,c(0,1,0),rawMesh$vb)
vb001 <- setPointsCone(-g3,c(0,0,1),rawMesh$vb)

# Join all points and triangles.
vbAll <- cbind(vb111, vb100, vb010, vb001)
itAll <- rbind(it, dim(rawMesh$vb)[1]+it,2*dim(rawMesh$vb)[1]+it,3*dim(rawMesh$vb)[1]+it)

# Visualization
plot_ly(
  x = vbAll[1,], y = vbAll[2,], z = vbAll[3,],
  i = itAll[,1], j = itAll[,2], k = itAll[,3],
  facecolor = toRGB(rep(c("#0D9127", "#0D6591", "#0D6591", "#0D6591"), each = n^2), alpha = 0.5),
  type = "mesh3d", flatshading=TRUE
)
# end fig --------------------------------------------

##########################################
#  DRAW JUST ONE CONE
##########################################


# Gamma values
(g3 <- t(matrix(c(-1,.3,-.3,.3,-1,-.3,-.3,.3,-1),3,3)))

# begin fig --------------------------------------------
# We set n to divide each cone into n^2 triangles
# The higher value, the better visualization
n <- 200
# Raw mesh that will be instantiated in each cone
rawMesh <- BuildMeshPattern(n)
# Information of triangles:
it <- rawMesh$it-1

# Instantiate the cone
cone <- c(1,1,1)
vbCone <- setPointsCone(-g3,cone,rawMesh$vb)

# Visualization
coneColor <- "#0D9127"   # Color of the cone
planeColor <- "#606060"  # Color of the planes
alphaCone <- 1           # Transparency of the cone
alphaPlanes <- 0.6       # Transparency of the planes
pPlane <- t(matrix(c(1,1,0,1,-1,0,-1,1,0,-1,-1,0), 3, 4))
tPlane <- t(matrix(c(1,2,3,2,3,4), 3, 2))-1

plot_ly(
  x = vbCone[1,], y = vbCone[2,], z = vbCone[3,],
  i = it[,1], j = it[,2], k = it[,3],
  facecolor = toRGB(rep(coneColor, n^2), alpha = alphaCone),
  type = "mesh3d", 
  flatshading = TRUE,
  lighting = list(ambient = 0.6,
                  diffuse = 0.8,
                  fresnel = 0.1,
                  specular = 0.5,
                  roughness = 0.5,
                  facenormalsepsilon = 0,
                  vertexnormalsepsilon = 0)) %>% 
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


##########################################
#  BIODIVERSITY SUB-CONES FOR 3 SPECIES
##########################################

(g3 <- t(matrix(c(-1,.3,.3,.3,-1,.3,.3,.3,-1),3,3)))

# begin fig --------------------------------------------
# Quick with poorer quality
sph <- sphereSubCones(g3,      # Gamma matrix
                      c(1,1,1), # Cone
                      60,       # Initial division on triangles
                      0.01)     # Stop dividing triangles with vertices closer
# Slower but very good quality:
#sph <- sphereSubCones(g3, c(1,1,1),120,0.003)

# Get the points (vertices of the triangles)
dataPoints <- as.matrix(t(sph$points[,1:3]))
# Triangles to be plotted are those with a color
nonZeroTriang <- as.matrix(sph$triang[sph$triang[,4] != "no",1:3])-1
# The lines below are to convert codes representing ISs into random colors
# change for other options
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
  hovertemplate = paste('x = %{x:.2f}',
                        '<br>y = %{y:.2f}',
                        '<br>z = %{z:.2f}'),
  lighting = list(ambient = 0.6,
                  diffuse = 0.8,
                  fresnel = 0.1,
                  specular = 0.5,
                  roughness = 0.5,
                  facenormalsepsilon = 0,
                  vertexnormalsepsilon = 0)) %>%
  layout(scene = list(
           xaxis = list(visible = FALSE),
           yaxis = list(visible = FALSE),
           zaxis = list(visible = FALSE)#,
           #aspectmode = 'manual'  # Usar proporciones manuales
           #aspectratio = list(x = 1, y = 1, z = 1)  # Proporciones iguales
         ))
# end fig --------------------------------------------
