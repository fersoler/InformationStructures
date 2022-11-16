# R Code for the article "Permanence via invasion graphs: Incorporating community assembly into Modern Coexistence Theory" by Josef Hofbauer and Sebastian J. Schreiber in the Journal of Mathematical Biology.

# R Code Author: Sebastian J. Schreiber

# This file illustrates the use of the main functions in main_functions.R. The example corresponds to the LV system used in the top row of Figure 1 in the article.

# load the invasion_graph_main_function.R commands.

source("tests/invasionGraphs/invasion_graph_main_functions.R")

# set seed for choosing the matrices and vectors for the Lotka Volterra model
seed=26
set.seed(seed)
# define the number of species and the matrices for the Lotka-Volttera model dx/dt=x*(Ax+b)
k=5
A=-diag(k)-1.5*matrix(runif(k^2),k,k)
b=matrix(1,k,1)
# compute the invasion scheme
IS=LV.IS(A,b)
# compute the invasion graph
out=IG.function(IS)
# use the plot.IG command to plot the figure with default settings
plot.IG(out)

## Fernando:
### Compare with our IG and IS
source("R/ISgraph.R")
source("R/ISmeasures.R")
library(tidyverse)

(g <- data.matrix(t(matrix(c(
  -1, -0.16, -0.36, 0.18, -0.04,
  0.27, -1, -1.13, -0.19, -0.20,
  .25, -0.27, -1, 0.06, -0.23,
  0.12, 0.15, -0.16, -1, .35,
  -0.08, -0.33, .19, 0.20, -1),nrow = 5))))
(a <- as.data.frame(t(c(1,2,-1.3,-.1,1))))
#(a <- as.data.frame(t(c(1,-1,1,-1,1))))


#(g <- data.matrix(matrix(0.1, nrow = 5,ncol = 5) - 1.1 * diag(5)))
#(a <- as.data.frame(t(c(1,-.3,1,-.2,1))))

# Out Information Structure:
IS5 <- ISbuild(a,g)
gr5 <- ISgraph(IS5,1:5)
ISgraphDrawLabels(IS5,gr5, ISgraphLayout(IS5, gr5, "5Dfix"))

# Invasion Graph:
IS=LV.IS(g,as.matrix(t(a)))
# compute the invasion graph
out=IG.function(IS)
# use the plot.IG command to plot the figure with default settings
plot.IG(out)

