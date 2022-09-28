# Building and drawing the IS graphs

library(igraph)
library(relations)
library(R6)
source("R/ISbuild.R")

########################################################################
# Obtaining the praph (with igraph)
########################################################################

# Function to build the IS graph
# Input:
# - IS: Information Structure (as obtained by ISbuild)
# - sNames: array with the labels of the species in the system
# (numbers, letters, etc)
# Output: list with the following elements
# - $graph: graph of the IS created by using the igraph library
# - $vlab: labels of the graph vertices
ISgraph <- function(IS, sNames){
  list(
    graph = graph_from_adjacency_matrix(IS$connectivity),
    vlab = nodeStrings(t(IS$points),sNames)
  )
}

# Get the list of vertex labels for all the nodes in the IS
nodeStrings <-function(isPoints, sNames){
  strs <- array("",dim(isPoints)[2])
  for(i in 1:length(strs)){
    p<-paste(sNames[isPoints[,i]>0],sep = "",collapse = "")
    if(p=="")
      p<-0
    strs[i]<-p
  }
  strs
}

########################################################################
# Setting the graph layout
########################################################################


# Function to set the graph layout for drawing it
# Input:
# - IS: Info. Structure (output of ISbuild)
# - ISgr: IS graph (output of ISgraph)
# - lyType: type of graph layout
#   * "tree": use the igraph option layout_as_tree (default)
#   * "3Dcube": for 3 species, nodes in the vertices of a cube
#   * "5Dfix": for 5 species, nodes in fixed positions
# Output: graph layout (matrix with node positions)
ISgraphLayout <- function(IS, ISgr, lyType = "tree"){
  layt <- layout_(ISgr$graph,as_tree())
  if(lyType == "tree"){
    conn <- IS$connectivity
    rownames(conn) <- rownames(IS$points)
    colnames(conn) <- rownames(IS$points)
    reducedConn <- relation_incidence(transitive_reduction(
      as.relation(conn)))[rownames(conn), colnames(conn)]
    grNew <- graph_from_adjacency_matrix(reducedConn)
    layt = layout_(grNew,as_tree())
  }
  if(lyType == "3Dcube" && dim(IS$points)[2] == 3){
    layt <- getCoordIS3D(t(IS$points))
  }
  if(lyType == "5Dfix" && dim(IS$points)[2] == 5){
    layt <- t(apply(IS$points, 1,
                             function(x) as.matrix(fixedNodes5SCoords[fixedNodes5SCoords$subset==toString(c(1:5)[x>0]),][1,2:3])))
  }
  layt
}

# Coordinates of the nodes in 3D graphs
vertexCoords3D<-cbind("0, 0, 0"=c(0,0),
                      "1, 0, 0"=c(2,0),
                      "0, 1, 0"=c(0,2),
                      "1, 1, 0"=c(2,2),
                      "0, 0, 1"=c(-1,-1.5),
                      "1, 0, 1"=c(1,-1.5),
                      "0, 1, 1"=c(-1,0.5),
                      "1, 1, 1"=c(1,0.5))
# Get coordinates for a node in the IS
getCoordNode3D <- function(node){
  vertexCoords3D[,toString(as.integer(node>0))]
}
getCoordIS3D <- function(isP){
  t(apply(isP,2,getCoordNode3D))
}

# fixed positions of the nodes in IS with 5 species
fixedNodes5S <- R6::R6Class(

  classname = "fixedNodes5S",

  public = list(

    pCoords = NULL,

    initialize = function(){
      # First, the fixed positions for all possible points in IS of size 5
      yVals <-  c(4,8,12,16,10)
      xVals <- c(3,5.5,5.5,3)
      fc <- 3
      # Level 0
      pointsCoords <- data.frame(subset="", x=0, y=0, stringsAsFactors=FALSE )
      # Levels 1:4
      for(lev in 1:4){
        allLevN = combn(1:5,lev)
        for(i in 1:ncol(allLevN)){
          pointsCoords[nrow(pointsCoords)+1,] <- list(subset=toString(allLevN[,i]),
                                                    x=fc*(i-xVals[lev]),
                                                    y=yVals[lev])
        }
      }
      # Level 5
      pointsCoords[nrow(pointsCoords)+1,]<-list(subset=toString(1:5), x=0, y=20)
      self$pCoords <- pointsCoords
    }
  )
)
fixedNodes5Spoints <- fixedNodes5S$new()
fixedNodes5Spoints$initialize()
fixedNodes5SCoords <- fixedNodes5Spoints$pCoords

#############################################################
# Drawing IS graphs
#############################################################

# Some example functions, more can be defined
# More arguments can be added to the functions when necessary to
# set node size, colors, etc.

# A function to plot a graph with labeled node
# Input:
# - IS: Info. Structure (output of ISbuild)
# - ISgr: IS graph (output of the ISgraph function)
# - ISlay: layout of the graph (output of ISgraphLayout)
ISgraphDrawLabels <- function(IS, ISgr, ISlay){
  plot(ISgr$graph,
       vertex.label = ISgr$vlab,
       layout = ISlay,
       edge.width=1,
       vertex.shape = "circle",
       vertex.color = "white",
       vertex.label.cex = 0.6,
       mark.groups=list(c(IS$gassInd)),
       mark.col="orange",
       mark.expand=4,
       mark.shape=-2)
}

# A function to plot a graph with pie charts in nodes
# Input:
# - IS: Info. Structure (output of ISbuild)
# - ISgr: IS graph (output of the ISgraph function)
# - ISlay: layout of the graph (output of ISgraphLayout)
# - colors: colors representing the species
# - maxAbund: if 0, then all nodes have the same size, if -1, the size is
# relative to the max sum of abundances in the IS points, otherwise
# the node sizes are relative to the vlaue of maxAbund
ISgraphDrawPie <- function(IS, ISgr, ISlay, colors, maxAbund = 0){
  gr <- ISgr$graph
  isP <- t(IS$points)
  abund <- colSums(isP)
  maxAbundThis <- max(abund)
  for(i in 1:ncol(isP)){
    if(maxAbund == 0){
      V(gr)[i]$size = 20
    }
    else{
      if(maxAbund == -1){
        mm = maxAbundThis
      }else{
        mm = maxAbund
      }
      V(gr)[i]$size = 10+10*(abund[i]/mm)
    }
    if(i==1){
      V(gr)[i]$shape = "circle"
      V(gr)[i]$color = "white"
    }
    else {
      if(sum(as.integer(isP[,i]>0))==1){
        V(gr)[i]$shape = "circle"
        V(gr)[i]$color = colors[isP[,i]>0]
      } else {
        V(gr)[i]$shape = "pie"
        V(gr)[i]$pie = list(10*(isP[,i]))
        V(gr)[i]$pie.color = list(colors)
      }
    }
  }
  plot(gr,
       layout = ISlay,
       edge.width=1,
       mark.groups=list(c(IS$gassInd)),
       mark.col="orange",
       mark.expand=4,
       mark.shape=-2,
       vertex.label=NA)
}





