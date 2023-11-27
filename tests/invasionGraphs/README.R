#################################################################
# Modification of Schreiber's code
# Author: Fernando Soler-Toscano - fsoler@us.es
#################################################################

# The following files:
# - 'invasion_graph_main_functions.R'
# - 'invasion_graph_example.R'
# come from:
# Schreiber, S. (2022). R Code for the article "Permanence via invasion graphs: 
# Incorporating community assembly into Modern Coexistence Theory" by 
# Josef Hofbauer and Sebastian J. Schreiber in the Journal of Mathematical 
# Biology. Zenodo. <https://doi.org/10.5281/zenodo.7111753> 
# Published under the Creative Commons Attribution 4.0 International license.

# File `invasion_graph_main_functions-mod.R` contains some modifications of the original file
# 'invasion_graph_main_functions.R' by S. Schreiber. 

# We present an example illustrating the problem solved with the modification.

source("tests/invasionGraphs/invasion_graph_main_functions.R")
source("tests/invasionGraphs/invasion_graph_main_functions-mod.R")
source("R/ISbuild.R")
source("R/ISgraph.R")

## Setting the parameters

## Set the parameters
A <- as.matrix(matrix(c(-1,1.5,1,-1,-1,0.5,-1,1,-1), ncol=3))
b <- c(1,1,1)

# IG with original code
IS=LV.IS(A,b)
out=IG.function(IS)
plot.IG(out)

# IG with modified code
ISComm=LV.ISandComm(A,b)
out=IG.functionComm(ISComm)
plot.IG(out)

# Observations

#  Schreiber code introduces species 1 in communities 12 and 13. 
# This is because the IG construction builds the communities from the zeros 
# in the IS (Invasion Scheme). However, in some cases these zeros do not 
# correspond to positive abundances in the system solutions.

# The case with community 12
I <- c(1,2) 
# Solution, only species 2 have positive abundance: 
(xtemp <- solve(A[I, I], -b[I]))
# Build the row in the IS
xtemp2 <- rep(0,ncol(A))
xtemp2[I] <- xtemp
# Resulting row, there is a 0 also for species 1. 
(rowIS <- A%*%xtemp2+b) 

# Look at the Information Structure created following the algorithm in
# https://arxiv.org/abs/2209.09802

IS3 <- ISbuild(as.data.frame(t(b)),A)
gr3 <- ISgraph(IS3, 1:3)
# Draw the IS in different ways:
ISgraphDrawLabels(IS3,gr3, ISgraphLayout(IS3, gr3, "3Dcube"))


