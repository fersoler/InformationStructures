source("R/ISbuild.R")
source("R/ISgraph.R")
source("R/rewiring.R")
source("R/figs/circleCones.R")
source("R/figs/sphereCones.R")
source("R/figs/sphereSubCones.R")
source("R/figs/InformationField.R")

#### 2D examples with increase of cooperation due to r_i decrease

drawCircleCones_with_R <- function(a, r){
  a2 <- keepDiagDomRows(a, rewiringIncreaseCoop(a,r))
  drawCircleCones(a2)
  point <- r/sqrt(sum(r^2))
  points(point[1],point[2], pch = 21, bg = "yellow", cex = 1)
}


(a <- matrix(c(-1,.2,.2,-1),2,2))
# R both positive
r <- c(1,1)
drawCircleCones_with_R(a, r)
# r_1 = 0
r <- c(0,1)
drawCircleCones_with_R(a, r)
# r_1 = -0.2
r <- c(-0.2,1)
drawCircleCones_with_R(a, r)
# r_1 = -0.4
r <- c(-0.4,1)
drawCircleCones_with_R(a, r)
# R both -0.5
r <- c(-0.5,-0.5)
drawCircleCones_with_R(a, r)
# R both -1
r <- c(-1,-1)
drawCircleCones_with_R(a, r)

## 3D examples with 2 competing species and one cooperating with them

# Gamma values
(a3 <- matrix(c(-1,-.3,.3,-.3,-1,.3,.3,.3,-1),3,3))
(a3b <- rewiringReduceComp(a3))
## Draw the sphere with some cones:

drawCircleCones(a3[1:2,1:2])
drawCircleCones_with_R(a3b[1:2,1:2])

# begin fig --------------------------------------------
# We set n to divide each cone into n^2 triangles
# The higher value, the better visualization
n <- 200
# Raw mesh that will be instantiated in each cone
rawMesh <- BuildMeshPattern(n)
# Information of triangles:
it <- rawMesh$it-1

# Instantiate the points for several cones:
vb111  <- setPointsCone(-a3,c(1,1,1),rawMesh$vb)
vb111b <- setPointsCone(-a3b,c(1,1,1),rawMesh$vb)

# Join all points and triangles.
vbAll <- cbind(vb111, 0.99*vb111b)
itAll <- rbind(it, dim(rawMesh$vb)[1]+it)

# Visualization
plot_ly(
  x = vbAll[1,], y = vbAll[2,], z = vbAll[3,],
  i = itAll[,1], j = itAll[,2], k = itAll[,3],
  facecolor = toRGB(rep(c("#0D6591","orange"), each = n^2), alpha = 1),
  type = "mesh3d", flatshading=TRUE
)
# end fig --------------------------------------------




