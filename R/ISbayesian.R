source("R/ISbuild.R")

# Function to build a bayesian IS
# Input:
# - listAlphas, list of alpha values
# - listGammas, list of gamma values
# The output is a list with the following components:
# - points: matrix with columns indicating:
#   * subset: nodes with non-zero values in each point
#   * sumAbund: sum of abundances of the species in all IS in which the point
#   appears.
#   * occur: nummber of IS in which the point accurs
#   * a1, a2, ..., an: sum of abundances of each specie in all occurrences of 
#   the point
# - connect: connectivity matrix of the bayesian IS
getISbayesian <- function(listAlphas, listGammas){
  
  # Set the first IS to get common data
  a1 <- as.data.frame(t(listAlphas[[1]]))
  size <- length(a1)
  g1 <- matrix(data.matrix(listGammas[[1]]),nrow = size)
  IS1 <- ISbuild(a1, g1)
  allSubs <- IS1$subsetGASS
  allSubs[,2] <- 0
  allSubs[,3] <- 0
  allSubs[,4:(4+size-1)] <- 0
  colnames(allSubs) <- c(c("subset", "sumAbund", "occur"), paste("a",1:size,sep = ""))
  newConn <- matrix(0, ncol = 2^size, nrow = 2^size)
  
  for(nis in 1:length(listAlphas)){
    allSubs[1,3] <- allSubs[1,3]+1 # Occurrence of 0
    theIS <- ISbuild(as.data.frame(t(listAlphas[[nis]])), 
                     matrix(data.matrix(listGammas[[nis]]),nrow = size))
    
    # Set keys, abundances and ocurrences
    keys <- rep(0,nrow(theIS$points)) # to store new positions of the IS points
    keys[1] <- 1
    if(nrow(theIS$points) > 1){
      for(node in 2:nrow(theIS$points)){
        found <- FALSE
        toLook <- 2
        while(!found){
          if(toString((1:size)[theIS$points[node,] > 0]) == theIS$subsetGASS[toLook,1]){
            found <- TRUE
          } else {
            toLook <- toLook+1
          }
        }
        keys[node] <- toLook
        allSubs[toLook,2] <- allSubs[toLook,2]+sum(theIS$points[node,])
        allSubs[toLook,3] <- allSubs[toLook,3]+1
        allSubs[toLook,4:(4+size-1)] <- 
          allSubs[toLook,4:(4+size-1)]+theIS$points[node,]
      }
    }
    # Now set the connectivity
    for(n in 1:nrow(theIS$points)){
      for(m in 1:nrow(theIS$points)){
        if(theIS$connectivity[n,m] == 1){
          newConn[keys[n],keys[m]] <- newConn[keys[n],keys[m]]+1
        }
      }
    }
    
  }
  
  # Positions with non-zero occurrences
  nonZeroPos <- (1:nrow(allSubs))[allSubs[,3]>0]
  
  # Get non-zero points and matrix
  points <- allSubs[nonZeroPos,]
  connect <- newConn[nonZeroPos,nonZeroPos] / points[,3] #/ points[1,3]
  
  list(points = points,
       connect = connect
  )
}

##################################################################
# Graphs and layouts
##################################################################

# Create a bayesian IS graph
bayISgraph <- function(bayIS, sNames){
  list(
    graph = graph_from_adjacency_matrix(bayIS$connect, mode="directed", weighted = TRUE),
    vlab = nodeStrings(t(bayIS$points[,4:ncol(bayIS$points)]),sNames)
  )
}

# Bayesian graph layout
bayISgraphLayout <- function(bayIS, bayISgr, lyType = "tree"){c
  
  cMat <- (bayIS$connect)
  
  quant <- apply(cMat,1,function(x) quantile(x,.5))
  cMat[cMat<quant] <- 0
  cMat[lower.tri(cMat)] <- 0
  cMat[cMat>0] <- 1
  transitive_reduction(as.relation(cMat))
  IS <- list(
    points = bayIS$points[,4:ncol(bayIS$points)],
    connectivity = cMat
  )
  ISgraphLayout(IS, bayISgr, lyType)
}

# Bayesian graph drawing (with labels)
bayesianISgraphDrawLabels <- function(bayIS, ISgr, ISlay){
  
  # To get the weights (to display)
  h <- graph.empty() + vertices(V(ISgr$graph))
  h <- h + edges(as.vector(t(get.edgelist(ISgr$graph)[order(E(ISgr$graph)$weight),])))
  E(h)$weight <- E(ISgr$graph)$weight[order(E(ISgr$graph)$weight)]
  
  # Curve edged (display)
  curves <-autocurve.edges2(ISgr$graph)
  
  # Plot:
  plot(ISgr$graph, 
       vertex.label = ISgr$vlab,
       vertex.size = 25,
       edge.width=.5+8*E(h)$weight, 
       edge.arrow.size=.01+10*E(h)$weight,
       vertex.color = paste("#000000",as.hexmode(floor(255*bay$points[,3]/bay$points[1,3])),sep=""),
       edge.curved=curves,
       vertex.label.cex = 1,
       vertex.label.color = "red",
       layout = ISlay
  )
 
}

# Drawing bayesian graphs (pie)
bayesianISgraphDrawPie <- function(IS, ISgr, ISlay, colors){
  
  # To get the weights (to display)
  h <- graph.empty() + vertices(V(ISgr$graph))
  h <- h + edges(as.vector(t(get.edgelist(ISgr$graph)[order(E(ISgr$graph)$weight),])))
  E(h)$weight <- E(ISgr$graph)$weight[order(E(ISgr$graph)$weight)]
  
  gr <- ISgr$graph
  isP <- t(IS$points[,4:ncol(IS$points)])
#  abund <- colSums(isP)
#  maxAbundThis <- max(abund)
  for(i in 1:ncol(isP)){
    V(gr)[i]$size = 10+20*(IS$points[i,3]/IS$points[1,3])
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
  
  # Curve edged (display)
  curves <-autocurve.edges2(ISgr$graph)
  
  plot(gr,
       layout = ISlay,
       edge.width=.3+3*E(h)$weight, 
       edge.arrow.size=.01+10*E(h)$weight,
       edge.width=1,
       edge.curved=curves,
       vertex.label=NA)
}

# Function to curve the edges that are duplicated
# From: https://stackoverflow.com/questions/16875547/
# using-igraph-how-to-force-curvature-when-arrows-point-in-opposite-directions
autocurve.edges2 <-function (graph, start = 0.3)
{
  cm <- count.multiple(graph)
  mut <-is.mutual(graph)  #are connections mutual?
  el <- apply(get.edgelist(graph, names = FALSE), 1, paste,
              collapse = ":")
  ord <- order(el)
  res <- numeric(length(ord))
  p <- 1
  while (p <= length(res)) {
    m <- cm[ord[p]]
    mut.obs <-mut[ord[p]] #are the connections mutual for this point?
    idx <- p:(p + m - 1)
    if (m == 1 & mut.obs==FALSE) { #no mutual conn = no curve
      r <- 0
    }
    else {
      r <- seq(-start, start, length = m)
    }
    res[ord[idx]] <- r
    p <- p + m
  }
  res
}
