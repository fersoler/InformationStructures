source("tests/invasionGraphs/invasion_graph_main_functions.R")
source("R/ISgraph.R")
source("R/ISbuildImproved.R")
source("R/ISmeasures.R")
library(tidyverse)

(g <- structure(c(-2.24677441129461, -0.231986807892099, -2.3462732497137,
                 -1.14973693038337, -2.62612960184924, -1.14922749949619, -1.36730381241068,
                 -1.59987767413259, -1.30858003953472, -1.09168730839156, -1.16794390231371,
                 -2.79106080019847, -1.87226756103337, -1.44121627137065, -0.464399928459898,
                 -2.66434512590058), .Dim = c(4L, 4L)))


(a <- as.data.frame(t(structure(c(0.264781965175644, 0.584047143114731, 0.267519016517326,
                 0.171760221244767), .Dim = c(4L, 1L)))))

eigen(g,only.values = TRUE)

# Out Information Structure:
(IS5 <- ISbuildThird(a, g))
gr5 <- ISgraph(IS5,1:4)
ISgraphDrawLabels(IS5,gr5, ISgraphLayout(IS5, gr5, with_kk()))

# Invasion Graph:
IS=LV.IS(g,as.matrix(t(a)))
# compute the invasion graph
(out=IG.function(IS))
# use the plot.IG command to plot the figure with default settings
plot.IG(out)

# non-permanent point 1,3,4
p <- c(1,3,4)
# There is a solution
solve(g[p,p],-a[p])
# But it's not the GASS
getGASS_LCP_Lemke(a[p], g[p,p])

# A permanent exmaple
p <- c(1,4)
# There is a solution
solve(g[p,p],-a[p])
# it's also the GASS
getGASS_LCP_Lemke(a[p], g[p,p])

