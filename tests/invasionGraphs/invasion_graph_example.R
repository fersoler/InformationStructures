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