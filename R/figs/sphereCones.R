#################################################################
# Functions to draw the biodiversity cones 
# Author: Fernando Soler-Toscano - fsoler@us.es
#################################################################

library(plotly)

# Build a generic mesh pattern to divide a biodiversity cone into n^2 triangles
# The output of this function is a list with two matrices representing a mesh:
# * 'vb' is the list of points coordinates. To be instantiated for a real cone
# and normalized. If the triangle is divided into n^2 subtriangles, coordinate
# values are all list (a,b,c) of integer values such that a+b+c=n.
# Points in 'vb' no not represent real coordinates but proportions of vectors
# (v1,v2,v3) defining a cone. So, (a,b,c) represents the point:
# a*v1 + b*v2 + c*v3. In order to get a point in the surface of a sphere of
# radius 1, the vector will be divided by its norm.
# * 'it' is the list of triangles, each row is (a,b,c) indicating that a
# triangle is built with points a, b, and c in vb. Indexed start from 1, so to
# plot the triangle it-1 has to be used.
BuildMeshPattern <- function(n){

  # Number of points in the mesh
  npoints <- (n+1)*(n+2)/2

  # Matrices to store the points and triangles
  vb <- matrix(0, nrow=npoints, ncol = 3) # Point coordinates
  it <- matrix(0, nrow=(n**2), ncol = 3) # Triangles

  # First point
  vb[1,] <- c(0,n,0)

  nextP <- 2      # Number of the following point to be added
  nextT <- 1      # Number of the following triangle to be added
  prevPs <- c(1)  # Points that were added in the previous iteration

  # Iteration over layers of increasing number of points and triangles
  # The number of points to be added in layer l is equal to l
  for(l in 2:(n+1)){

    currPs <- c(nextP:(nextP+l-1)) # Points numbers to be added

    # Vector for the first point
    value <- c(l-1,n-l+1,0)

    # Adding points and triangles
    # Interation over the new points pn to be added
    for(pn in 1:l){

      # A new point is stored in position 'nextP' with the coordinated in value
      vb[nextP,] <- value

      # Cases in which triangles are added
      if(pn<l-1){
        it[nextT,] <- c(currPs[pn],prevPs[pn],currPs[pn+1])
        nextT <- nextT+1
        it[nextT,] <- c(prevPs[pn],currPs[pn+1],prevPs[pn+1])
        nextT <- nextT+1
      } else if(pn<l){
        it[nextT,] <- c(currPs[pn],prevPs[pn],currPs[pn+1])
        nextT <- nextT+1
      }

      # Preparing next point
      value <- value+c(-1,0,1)  # Value of the following point
      nextP <- nextP+1          # Number of the following point
    }

    # At the end the row of previous points (used to create triangles) is
    # assigned the row of point just added.
    prevPs <- currPs
  }

  # Return the list of point coordinated (to be instantiated and normalized) and
  # the list of triangles
  list(vb = vb, it = it)
}

# Function to normalize a vector. This is used because coordinates in 'vb' do
# not represent real points but coeficients of the vectors defining a cone.
normalize <- function(vec){
  vec/sqrt(sum(vec^2))
}


# A function that receives the three vertices of a triangle
# (one per column of 'verts') and instantiates the
# coordinates of the points forming the mesh (vb part).
# The output of the function is a matrix with the same
# format that 'vb' but instantiated for the indicated triangle.
setPointsTriangle <- function(verts, vb){
  apply(vb,1,function(x) {normalize(verts %*% x)})
}

#######################################
# LV cones
#######################################

# Function to instantiate and normalize the points in a mesh (vb part) for a
# specific cone given a gamma matrix. The output of the function is a matrix
# with the same format that 'vb' but instantiated for the indicated cone.
# The input arguments are:
# * gamma: connectivity matrix
# * cone: a vector of 1s (specie in the cone) and 0s (specie out of the cone).
# * vb: the list of generic points to be instantiated
setPointsCone <- function(gamma, cone, vb){
  # 1s in cone will use columns of gamma
  # 0s in cone will use columns of -diag(3)
  vectors <- -diag(3)
  for(i in 1:3){
    if(cone[i]==1){
      vectors[,i] <- normalize(gamma[,i])
    }
  }
  apply(vb,1,function(x) {normalize(vectors %*% x)})
}
