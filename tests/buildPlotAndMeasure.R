source("R/ISgraph.R")
source("R/ISmeasures.R")
source("R/ISbayesian.R")
source("R/ISbuildImproved.R")

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

# Get the IS with improved algorithm
IS3b <- ISbuildThird(alphas3[t,],gammas3)
gr3b <- ISgraph(IS3b,1:3)
# And the measures
getISmeasures(IS3b,gr3b, alphas3[t,],gammas3)


# Different ways to show the IS
ISgraphDrawLabels(IS3,gr3, ISgraphLayout(IS3, gr3, "tree"))
ISgraphDrawLabels(IS3,gr3, ISgraphLayout(IS3, gr3, "3Dcube"))
ISgraphDrawPie(IS3,gr3, ISgraphLayout(IS3, gr3, "tree"), c("#806000", "#002080", "#408000"))
ISgraphDrawPie(IS3b,gr3b, ISgraphLayout(IS3b, gr3b, "tree"), c("#806000", "#002080", "#408000"))

ISgraphDrawPie(IS3,gr3, ISgraphLayout(IS3, gr3, "3Dcube"), c("#806000", "#002080", "#408000"))
ISgraphDrawPie(IS3b,gr3b, ISgraphLayout(IS3b, gr3b, "3Dcube"), c("#806000", "#002080", "#408000"))

###########################################################################
# Bayesian IS
###########################################################################


getRandomAlphas <- function(n){
  alp <- rnorm(n)
  alp/sqrt(sum(alp^2))
}

# 3 nodes
#set.seed("20220929")
set.seed("20220927")
gFix <- matrix(rnorm(9)/20,3,3)*(1-diag(3))-diag(3)
aFix <- getRandomAlphas(3)/3

aL <- list()
gL <- list()
init <- 1
end <- 20
for(i in init:end){
  aL[i+1-init] <- as.data.frame(aFix + getRandomAlphas(3)/2)
  gL[i+1-init] <- as.data.frame(matrix(gFix + (matrix(rnorm(9)/3,3,3)*(1-diag(3)))))
}


# Build bayesian IS
bay <- getISbayesian(aL, gL)
bay
# Build the graph
grafo <- bayISgraph(bay, 1:3)

bayesianISgraphDrawLabels(bay, grafo, bayISgraphLayout(bay, grafo, lyType = "3Dcube"))
bayesianISgraphDrawLabels(bay, grafo, bayISgraphLayout(bay, grafo, lyType = "tree"))
bayesianISgraphDrawPie(bay, grafo, bayISgraphLayout(bay, grafo, lyType = "3Dcube"),c("#806000", "#002080", "#408000"))
bayesianISgraphDrawPie(bay, grafo, bayISgraphLayout(bay, grafo, lyType = "tree"),c("#806000", "#002080", "#408000"))


# 5 nodes
set.seed("20220929")
gFix <- matrix((rnorm(25)-0.3)/25,5,5)*(1-diag(5))-diag(5)
aFix <- getRandomAlphas(5)/5

aL <- list()
gL <- list()
init <- 1
end <- 10
for(i in init:end){
  aL[i+1-init] <- as.data.frame(aFix + getRandomAlphas(5)/2)
  gL[i+1-init] <- as.data.frame(matrix(gFix + (matrix(rnorm(25)/5,5,5)*(1-diag(5)))))
}


# Build bayesian IS
bay <- getISbayesian(aL, gL)
# Build the graph
grafo <- bayISgraph(bay, 1:5)



bayesianISgraphDrawLabels(bay, grafo, bayISgraphLayout(bay, grafo, lyType = "5Dfix"))
bayesianISgraphDrawLabels(bay, grafo, bayISgraphLayout(bay, grafo, lyType = "tree"))
bayesianISgraphDrawPie(bay, grafo, 2*bayISgraphLayout(bay, grafo, lyType = "5Dfix"),
                       c("#806000", "#002080", "#408000", "#800060", "#008080"))
bayesianISgraphDrawPie(bay, grafo, 2*bayISgraphLayout(bay, grafo, lyType = "tree"),
                       c("#806000", "#002080", "#408000", "#800060", "#008080"))


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

### Improved version
IS5b <- ISbuildThird(a,g)
gr5b <- ISgraph(IS5b,1:5)
getISmeasures(IS5b, gr5b, a, g)



ISgraphDrawLabels(IS5,gr5, ISgraphLayout(IS5, gr5, "tree"))
ISgraphDrawLabels(IS5,gr5, ISgraphLayout(IS5, gr5, "5Dfix"))
ISgraphDrawPie(IS5,gr5, ISgraphLayout(IS5, gr5, "tree"),
               c("#806000", "#002080", "#408000", "#800060", "#008080"))
ISgraphDrawPie(IS5,gr5, ISgraphLayout(IS5, gr5, "5Dfix"),
               c("#806000", "#002080", "#408000", "#800060", "#008080"))
ISgraphDrawPie(IS5b,gr5b, ISgraphLayout(IS5b, gr5b, "5Dfix"),
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


