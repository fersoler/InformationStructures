source("R/ISgraph.R")
source("R/ISmeasures.R")

library(tidyverse)

#####################################################
# 3 species examples
#####################################################

# Load files
alphas3 <- as.data.frame(read_csv("data/alphas.csv"))
gammas3 <- as.matrix(read_csv("data/gMatrix.csv",col_names = FALSE))

# Set time (change only this)
t <- 23

# Get the IS
IS3 <- ISbuild(alphas3[t,],gammas3)
# And the graph
gr3 <- ISgraph(IS3,1:3)

# IS measures:
getISmeasures(IS3,gr3, alphas3[t,],gammas3)

# Different ways to show the IS
ISgraphDrawLabels(IS3,gr3, ISgraphLayout(IS3, gr3, "tree"))
ISgraphDrawLabels(IS3,gr3, ISgraphLayout(IS3, gr3, "3Dcube"))
ISgraphDrawPie(IS3,gr3, ISgraphLayout(IS3, gr3, "tree"), c("#806000", "#002080", "#408000"))
ISgraphDrawPie(IS3,gr3, ISgraphLayout(IS3, gr3, "3Dcube"), c("#806000", "#002080", "#408000"))


# Join different IS to get a bayesian IS

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
  nonZeroPos <- (1:8)[allSubs[,3]>0]
  
  # Get non-zero points and matrix
  points <- allSubs[nonZeroPos,]
  connect <- newConn[nonZeroPos,nonZeroPos] / points[,3] #/ points[,3]
  
  list(points = points,
       connect = connect
       )
}

aL <- list()
gL <- list()
init <- 1
end <- 250
for(i in init:end){
  aL[i+1-init] <- as.data.frame(t(alphas3[i,]))
  gL[i+1-init] <- as.data.frame(matrix(gammas3))
}

# Test
(bay <- getISbayesian(aL, gL))



(gr <- graph_from_adjacency_matrix(bay$connect, mode="directed", weighted = TRUE))

h <- graph.empty() + vertices(V(gr))
h <- h + edges(as.vector(t(get.edgelist(gr)[order(E(gr)$weight),])))
E(h)$weight <- E(gr)$weight[order(E(gr)$weight)]

plot(gr, 
     vertex.label = str_replace_all(bay$points[,1], ", ", ""),
     vertex.size = 25,
     #vertex.size = 10+(10*bay$points[,3])/bay$points[1,3],
     edge.width=5*E(h)$weight, 
     edge.arrow.size=E(gr)$weight,
     vertex.color = paste("#00FF00",as.hexmode(floor(255*bay$points[,3]/bay$points[1,3])),sep=""),
     curved = curve_multiple(gr, start = 0.5),
     vertex.label.cex = 1,
     vertex.label.color = "black",
     layout = 2*t(apply(bay$points[,4:6],1,getCoordNode3D))
     )


#####################################################
# 5 species examples
#####################################################

## Cooperative/competitive   ####################################
(g <- data.matrix(t(matrix(c(
  -1, -0.16, -0.36, 0.18, -0.04,
  0.27, -1, -1.13, -0.19, -0.20,
  .25, -0.27, -1, 0.06, -0.23,
  0.12, 0.15, -0.16, -1, .35,
  -0.08, -0.33, .19, 0.20, -1),nrow = 5))))
(a <- as.data.frame(t(c(1,-.2,1.3,-.1,1))))


IS5 <- ISbuild(a,g)
gr5 <- ISgraph(IS5,1:5)


# IS measures
getISmeasures(IS5, gr5, a, g)


ISgraphDrawLabels(IS5,gr5, ISgraphLayout(IS5, gr5, "tree"))
ISgraphDrawLabels(IS5,gr5, ISgraphLayout(IS5, gr5, "5Dfix"))
ISgraphDrawPie(IS5,gr5, ISgraphLayout(IS5, gr5, "tree"),
               c("#806000", "#002080", "#408000", "#800060", "#008080"))
ISgraphDrawPie(IS5,gr5, ISgraphLayout(IS5, gr5, "5Dfix"),
               c("#806000", "#002080", "#408000", "#800060", "#008080"))



## Purely cooperative  ###############################################
(g <- data.matrix(matrix(0.1, nrow = 5,ncol = 5) - 1.1 * diag(5)))
(a <- as.data.frame(t(c(1,-.3,1,-.2,1))))

IS5 <- ISbuild(a,g)
gr5 <- ISgraph(IS5,1:5)

# IS measures
getISmeasures(IS5, gr5, a, g)



ISgraphDrawLabels(IS5,gr5, ISgraphLayout(IS5, gr5, "tree"))
ISgraphDrawLabels(IS5,gr5, ISgraphLayout(IS5, gr5, "5Dfix"))
ISgraphDrawPie(IS5,gr5, ISgraphLayout(IS5, gr5, "tree"),
               c("#806000", "#002080", "#408000", "#800060", "#008080"))
ISgraphDrawPie(IS5,gr5, ISgraphLayout(IS5, gr5, "5Dfix"),
               c("#806000", "#002080", "#408000", "#800060", "#008080"))

#########################################################################
# Python example (code on Plos Comp Biol paper)
#########################################################################

## Fix criticality to get the same value

# Connectivity matrix
(gammasEx = t(matrix(c(0,0.1,0.2,0.3,0.3,0,
                                   0.4,0.1,0.15,0.25,
                                   0,0.1,0.1,0.2,0.3,0),ncol = 4))-diag(4))


# Alpha values
(alphasEx = as.data.frame(t(c(-.1,-.1,3,-.2))))


(IS4 <- ISbuild(alphasEx,gammasEx))
(gr4 <- ISgraph(IS4,1:4))

# IS measures
getISmeasures(IS4, gr4, alphasEx, gammasEx)

ISgraphDrawLabels(IS4,gr4, ISgraphLayout(IS4, gr4, "tree"))


