source("tests/invasionGraphs/invasion_graph_main_functions.R")
source("tests/invasionGraphs/invasion_graph_main_functions-mod.R")
source("R/ISbuild.R")
source("R/ISgraph.R")
source("R/ISmeasures.R")
source("R/IGmeasures.R")


##########################################
#  Informational Structures for 2 species
##########################################

# Intrinsic growth rates:
(b <- c(2,1.8))
# Matrix of interactions:
(a <- matrix(c(-1,0.2,0.3,-1),2,2))

# Build the Informational Structure
(IS <- ISbuild(b,a))
gr <- ISgraph(IS, 1:2)
# Draw the IS in two different ways:
ISgraphDrawLabels(IS,gr, ISgraphLayout(IS, gr, "tree"))
ISgraphDrawPie(IS, gr, ISgraphLayout(IS, gr, "tree"), c("brown", "purple"))

##########################################
#  Informational Structures for 3 species
##########################################

# Matrix of interactions
(a <- matrix(c(-1,.3,-.3,-.3,-1,.3,-.3,.3,-1),3,3))
# Intrinsic growth rates
(b <- c(1,-.2,1.3))
# Informational Structure:
(IS3 <- ISbuild(b,a))
gr3 <- ISgraph(IS3, 1:3)
# Draw the IS in different ways:
ISgraphDrawLabels(IS3,gr3, ISgraphLayout(IS3, gr3, "tree"))
ISgraphDrawLabels(IS3,gr3, ISgraphLayout(IS3, gr3, "3Dcube"))
ISgraphDrawPie(IS3,gr3, ISgraphLayout(IS3, gr3, "tree"), c("#806000", "#002080", "#408000"))
ISgraphDrawPie(IS3,gr3, ISgraphLayout(IS3, gr3, "3Dcube"), c("#806000", "#002080", "#408000"))

#### Example with cycles
# Intrinsic growth rates:
(b <- c(1,1,1))
# Matrix of interactions:
d <- -.5
(a <- matrix(c(d,1,-1,-1,d,1,1,-1,d),3,3))
(IS3 <- ISbuild(b,a))
gr3 <- ISgraph(IS3, 1:3)
# Draw the IS
ISgraphDrawLabels(IS3,gr3,layout_(gr3$graph,nicely()))


##########################################
#  Informational Structures for 4 species
##########################################

# Matrix of interactions
(a <- matrix(c(-1,.2,.2,-.3,.2,-1,-.3,.2,.2,.3,-1,.2,.3,-.2,-.3,-1),4,4))
# Intrinsic growth rates
(b <- c(1,-.2,1.3,1))
# Informational Structure:
(IS4 <- ISbuild(b,a))
gr4 <- ISgraph(IS4, 1:4)
# Draw the IS in different ways:
ISgraphDrawLabels(IS4,gr4, ISgraphLayout(IS4, gr4, "tree"))
ISgraphDrawPie(IS4,gr4, ISgraphLayout(IS4, gr4, "tree"), c("#806000", "#002080", "#408000","#0066FF"))

##########################################
# Comparison with invasion graphs
##########################################

(a <- data.matrix(t(matrix(c(
  -1,     0.16, -0.36, -0.18, -0.04,
   0.27, -1,    -0.13,  0.19,  0.20,
  -0.1,   0.27, -1,     0.06, -0.23,
   0.12,  0.15,  0.16, -1,     0.24,
  -0.08,  0.33, -0.19,  0.20, -1),nrow = 5))))
(b <- c(1,-0.3,1,-1,1.49))

# Information Structure:
(IS5 <- ISbuild(b, a))
gr5 <- ISgraph(IS5,1:5)
ISgraphDrawLabels(IS5,gr5, ISgraphLayout(IS5, gr5, "5Dfix"))
ISgraphDrawLabels(IS5,gr5, ISgraphLayout(IS5, gr5, "tree"))
ISgraphDrawPie(IS5,gr5, ISgraphLayout(IS5, gr5, "tree"), c("#806000", "#002080", "#408000","#0066FF", "#FFDD33"))

# Invasion Graph:
IG=LV.IS(a,b)
# compute the invasion graph
out=IG.function(IG)
# use the plot.IG command to plot the figure with default settings
plot.IG(out)

# Function to sort the communities of and Information Structures
# like a given Invasion Graph. 
sortISlikeIG <- function(IS, IG){
  species <- ncol(IS$points)
  newOrder <- rep(0, nrow(IS$points))
  boolM <- (IS$points > 0)
  p <- 1
  for(point in IG$composition){
    pattern <- rep(FALSE, species)
    pattern[point] <- TRUE
    key <- which(apply(boolM, 1, function(x) all(x == pattern)))
    newOrder[p] <- key
    p <- 1+p
  }
  subsetGASS <- IS$subsetGASS
  subsetGASS[,2] <- newOrder[IS$subsetGASS[,2]]
  list(
    points = IS$points[newOrder,],
    subsetGASS = subsetGASS,
    connectivity = IS$connectivity[newOrder, newOrder],
    gassInd = newOrder[IS$gassInd]
  )
}

## Compare Information Structure and Invasion Graph: 
# Sort the previous IS
IS5sorted <- sortISlikeIG(IS5, out)
# Compare the connectivity
identical(out$IG, IS5sorted$connectivity)

#####################################################
# Measuring the Invasion Graph
#####################################################
# We use the modification of Schreiber's code
# Look at the README.R file in '/tests/invasionGraphs'

# Build and display the IG:
ISComm=LV.ISandComm(a,b)
out=IG.functionComm(ISComm)
plot.IG(out)

# Nodes:
(points <- ISComm$Comm)
# Number of nodes:
(npoints <- nrow(points))
# Number of edges: 
(nedges <- sum(out$IG))
# Number of permanent communities:
sum(out$permanent)
# Node frondosity: 
IGnodeFrond(points)
# Edge frondosity
nedges/(npoints*(npoints-1)/2)
# Index of the GASS(es)
(gass <- IGgassIndex(out,points))    
# (mean) GASS size:
IGspeciesGASS(out,points)
# (mean) GASS abundances:
IGmeanAbundGASS(points,gass)
# (mean) Evenness of the GASS:
IGmeanEvenness(points,gass)

### Cycles with 5 species

d <- -1.0001
(a <- data.matrix(t(matrix(c(
  d, -1, 1, -1, 1,
  1, d, -1, 1, -1,
  -1, 1, d, -1, 1,
  1, -1, 1, d, -1,
  -1, 1, -1, 1, d),nrow = 5))))
(b <- c(1,1,1,1,1))
# Information Structure:
(IS5 <- ISbuild(b, a))
gr5 <- ISgraph(IS5,1:5)
# Draw the IS
ISgraphDrawLabels(IS5,gr5,layout_(gr5$graph,nicely()))

# Invasion graph: 
#IG=LV.ISandComm(a,b)
# compute the invasion graph
#out=IG.functionComm(IG)
# use the plot.IG command to plot the figure with default settings
#plot.IG(out)

### Cycles with 4 species

d <- -1.618033988748
(a <- data.matrix(t(matrix(c(
  d, 0, -1, 1,
  1, d, 0, -1,
  -1, 1, d, 0,
  0, -1, 1, d),nrow = 4))))
(b <- c(1,1,1,1))
# Information Structure:
(IS5 <- ISbuild(b, a))
gr5 <- ISgraph(IS5,1:4)
# Draw the IS
ISgraphDrawLabels(IS5,gr5,layout_(gr5$graph,in_circle()))
#ISgraphDrawLabels(IS5,gr5,layout_(gr5$graph,nicely()))

# Invasion graph: 
#IG=LV.ISandComm(a,b)
# compute the invasion graph
#out=IG.functionComm(IG)
# use the plot.IG command to plot the figure with default settings
#plot.IG(out)

