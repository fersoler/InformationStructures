#####################################################################
# Rewiring Structural Stability                                     #
# Examples                                                          #
#####################################################################

source("R/ISbuild.R")
source("R/ISgraph.R")
source("R/rewiring.R")
source("R/figs/circleCones.R")
source("R/figs/sphereCones.R")
source("R/figs/sphereSubCones.R")
source("R/figs/InformationField.R")


# 1. Two species reduce competition due to a commmon cooperator ------

# Matrix of interactions. Species 1 and 2 compete, but 3 is a common
# cooperator. 
(a3 <- matrix(c(-1,-.3,.3,-.3,-1,.3,.3,.3,-1),3,3))
# Modified version due to the coperation of species 3: 
(a3b <- rewiringReduceComp(a3))

# Look at the biodiversity cones of species 1 and 2
## Without common cooperator: 
drawCircleCones(a3[1:2,1:2])
## With common cooperator: 
drawCircleCones(a3b[1:2,1:2])

# Biodiversity cone for the 3 species: 
n <- 100
rawMesh <- BuildMeshPattern(n)
it <- rawMesh$it-1
# Without cooperation: 
vb111  <- setPointsCone(-a3,c(1,1,1),rawMesh$vb)
# With cooperation: 
vb111b <- setPointsCone(-a3b,c(1,1,1),rawMesh$vb)

# Join all points and triangles.
vbAll <- cbind(vb111, 0.99*vb111b)
itAll <- rbind(it, dim(rawMesh$vb)[1]+it)

# Visualization (blue for previous and orange with cooperation)
plot_ly(
  x = vbAll[1,], y = vbAll[2,], z = vbAll[3,],
  i = itAll[,1], j = itAll[,2], k = itAll[,3],
  facecolor = toRGB(rep(c("#0D6591","orange"), each = n^2), alpha = 1),
  type = "mesh3d", flatshading=TRUE
)



# 2. Two species increase cooperation when r_i decrease ------

# Function to draw 2d cones with a point representing r_i
drawCircleCones_with_R <- function(a, r){
  a2 <- keepDiagDomRows(a, rewiringIncreaseCoop(a,r))
  drawCircleCones(a2)
  point <- r/sqrt(sum(r^2))
  points(point[1],point[2], pch = 21, bg = "yellow", cex = 1)
}

# Matrix of interactions
(a <- matrix(c(-1,.2,.2,-1),2,2))
# R both positive
r <- c(1,1)
drawCircleCones_with_R(a, r)

# Now, r_1 = 0
r <- c(0,1)
drawCircleCones_with_R(a, r)


# Negative r_1 = -0.2
# The cone increase includes allows biodiversity
r <- c(-0.2,1)
drawCircleCones_with_R(a, r)

# Smaller r_1 = -0.4
# The cone increase is not enough
r <- c(-0.4,1)
drawCircleCones_with_R(a, r)

# Both r_i = -0.5
r <- c(-0.5,-0.5)
drawCircleCones_with_R(a, r)

