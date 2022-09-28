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


