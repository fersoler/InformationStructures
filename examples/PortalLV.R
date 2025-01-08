library(tidyverse)
library(plotly)
source("R/ISgraph.R")
source("R/figs/sphereCones.R")
source("R/figs/sphereSubCones.R")
source("R/setParameters.R")

#########################################################################
# Loading example data
#########################################################################
#
# Data from:
# https://zenodo.org/records/10892031
# 
# @Article{christensen2019portalr,
#   title = {portalr: an R package for summarizing and using the Portal
#     Project Data},
#   author = {Erica M. Christensen and Glenda M. Yenni and Hao Ye and
#     Juniper L. Simonis and Ellen K. Bledsoe and Renata M. Diaz and
#     Shawn D. Taylor and Ethan P. White and S. K. Morgan Ernest},
#   year = {2019},
#   journal = {Journal of Open Source Software},
#   volume = {4},
#   number = {33},
#   pages = {1098},
#   doi = {10.21105/joss.01098},
# }
#
# The data table was obtained from the 'portalr' package:
# https://cran.r-project.org/web/packages/portalr/index.html 
# Using: 
# abundance(".", level = "site", shape = "crosstab", time = "period")
# Three species are selected:
# DM: Dipodomys merriami
# DO: Dipodomys ordii
# DS: Dipodomys spectabilis
# The sole purpose of selecting these species is to illustrate how 
# to adjust LV parameters and visualise IS.

PortalR_abundances <- read.csv("data/PortalData_abundances.csv")


# Data table
data <- PortalR_abundances[,c(1,c(4,5,9))]
# Abbreviations of the three species
spNs <- colnames(data)[2:4]
# Avoid 0 values (because of the use of log)
data[,2:4] <-data[,2:4]+1
# Divide by 50 to scale resulting parameters 
#data[,2:4] <- (data[,2:4]/50)


#########################################################################
# Setting LV parameters
#########################################################################

# Obtain alpha and gamma parameters: 
winW <- 10 # Margins to set alpha parameters
params <- getLVparamsNonAuto(as.matrix(data[,2:4]),winW)
# Matrix of interspecific interactions: 
(gammas <- params$gammas)
# Intrinsic growth rates
alphas <- data.frame(times = data[1:(nrow(data)-1),1], A1 = params$alphas[,1], A2 = params$alphas[,2], A3 = params$alphas[,3])

# Look at the intrinsic growth rates: 
fig <- plot_ly(alphas, x = ~times, y = ~A1, name = spNs[1], 
               type = 'scatter', mode = 'lines+markers', color='A')
fig <- fig %>% add_trace(x = ~times, y = ~A2, name = spNs[2], 
                         type = 'scatter', mode = 'lines+markers', color='B')
fig <- fig %>% add_trace(x = ~times, y = ~A3, name = spNs[3], 
                         type = 'scatter', mode = 'lines+markers', color='C')
fig <- fig %>% layout(
  title = "LV intrinsic growth rates",
  xaxis = list(title = "Period (months)"),
  yaxis = list (title = "Intrinsic growth rate"))
fig

#########################################################################
# Reconstructing time series
#########################################################################

reconstructPoints <- NULL
for(pi in 1:nrow(alphas)){
  theS2 <- getEvolLV3(as.matrix(alphas[pi,2:4]), gammas, data[pi,2:4], 1000, data[pi,1], data[pi+1,1])
  reconstructPoints <- rbind(reconstructPoints, theS2[1000,])
}
dataF <- data.frame(x = data[,1], SP1 = data[,2], SP2 = data[,3], SP3 = data[,4])
fig <- plot_ly(dataF, x = ~x, y = ~SP1, name = spNs[1], 
               type = 'scatter', mode = 'line', color = 'A')
fig <- fig %>% add_trace(y = ~SP2, name = spNs[2], color = 'B',
                         type = 'scatter', mode = 'line')
fig <- fig %>% add_trace(y = ~SP3, name = spNs[3], color = 'C', 
                         type = 'scatter', mode = 'line')
fig <- fig %>% add_trace(x = reconstructPoints[,1], y = reconstructPoints[,2], 
                         name = paste('LV',spNs[1]), color = 'A',
                         type = 'scatter', mode = 'markers')
fig <- fig %>% add_trace(x = reconstructPoints[,1], y = reconstructPoints[,3], 
                         name = paste('LV',spNs[2]), color='B',
                         type = 'scatter', mode = 'markers')
fig <- fig %>% add_trace(x = reconstructPoints[,1], y = reconstructPoints[,4], 
                         name = paste('LV',spNs[3]), color = 'C',
                         type = 'scatter', mode = 'markers')
fig <- fig %>% layout(
  title = "Real and reconstructed (LV) data",
  xaxis = list(title = "Period (months)"),
  yaxis = list (title = "Abundance"))
fig


#########################################################################
# Looking at a particular IS
#########################################################################

# Number of period
pi <- 140
# Build the IS
IS3 <- ISbuild(as.matrix(alphas[pi,2:4]), gammas)
# Path to the next period: 
path <- getEvolLV3(as.matrix(alphas[pi,2:4]), gammas, data[pi,2:4], 100, 1, 2)
# Path in 1000 more periods: 
path2 <- getEvolLV3(as.matrix(alphas[pi,2:4]), gammas, data[pi,2:4], 10000, 1, 1000)

# Figure
scene = list(xaxis = list(title = spNs[1], scaleratio = 1),
             yaxis = list(title = spNs[2], scaleratio = 1),
             zaxis = list(title = spNs[3], scaleratio = 1))
points <- data.frame(SP1 = IS3$points[,1], SP2 = IS3$points[,2], SP3 = IS3$points[,3])
fig <- plot_ly(points, x = ~SP1, y = ~SP2, z = ~SP3, mode = 'markers', 
               type = 'scatter3d', 
               marker = list(size = 3, color = 'black'),showlegend=FALSE)
fig <- fig %>% add_trace(x = path2[11:10000,2], y = path2[11:10000,3], z = path2[11:10000,4], mode = 'lines',
                         marker = list(size = 1, color = 'brown'))
fig <- fig %>% add_trace(x = path[,2], y = path[,3], z = path[,4],
                         marker = list(size = 1, color = 'blue'))
fig <- fig %>% add_trace(x = path[1,2], y = path[1,3], z = path[1,4], 
                         marker = list(size = 3, color = 'red'))
fig <- fig %>% layout(title=paste("IS at time",data[pi,1]),showlegend=FALSE, scene = scene)
fig

# Looking at the IS
gr3 <- ISgraph(IS3, 1:3)
# Draw the IS in different ways:
#ISgraphDrawLabels(IS3,gr3, ISgraphLayout(IS3, gr3, "tree"))
#ISgraphDrawLabels(IS3,gr3, ISgraphLayout(IS3, gr3, "3Dcube"))
#ISgraphDrawPie(IS3,gr3, ISgraphLayout(IS3, gr3, "tree"), c("#806000", "#002080", "#408000"))
ISgraphDrawPie(IS3,gr3, ISgraphLayout(IS3, gr3, "3Dcube"),  c("#66c2a5", "#fc8d61", "#8da0cb"))


## Sphere with max bio cone (green) and cones with 1 species (blue)
# Intrinsic growth rates are represented by a point

# We set n to divide each cone into n^2 triangles
# The higher value, the better visualization
n <- 120
# Raw mesh that will be instantiated in each cone
rawMesh <- BuildMeshPattern(n)
# Information of triangles:
it <- rawMesh$it-1

# Instantiate the points for several cones:
vb111 <- setPointsCone(-gammas,c(1,1,1),rawMesh$vb)
vb100 <- setPointsCone(-gammas,c(1,0,0),rawMesh$vb)
vb010 <- setPointsCone(-gammas,c(0,1,0),rawMesh$vb)
vb001 <- setPointsCone(-gammas,c(0,0,1),rawMesh$vb)

# Join all points and triangles.
vbAll <- cbind(vb111, vb100, vb010, vb001)
itAll <- rbind(it, dim(rawMesh$vb)[1]+it,2*dim(rawMesh$vb)[1]+it,3*dim(rawMesh$vb)[1]+it)
# Visualization
plot_ly(
  x = vbAll[1,], y = vbAll[2,], z = vbAll[3,],
  i = itAll[,1], j = itAll[,2], k = itAll[,3],
  facecolor = toRGB(rep(c("#0D9107", "#66c2a5", "#fc8d61", "#8da0cb"), each = n^2), alpha = 0.5),
  type = "mesh3d", 
  flatshading=TRUE) %>% 
add_markers(x = 1.03*normalize(alphas[pi,2:4])[1,1],
            y = 1.03*normalize(alphas[pi,2:4])[1,2], 
            z = 1.03*normalize(alphas[pi,2:4])[1,3],
            marker = list(size = 4, color = 'yellow')) %>% 
layout(showlegend=FALSE, scene = scene)


#### All alpha values

normAlphas <- as.data.frame(t(apply(alphas[, 2:4], 1, normalize)))
plot_ly(
  x = vbAll[1,], y = vbAll[2,], z = vbAll[3,],
  i = itAll[,1], j = itAll[,2], k = itAll[,3],
  facecolor = toRGB(rep(c("#0D9107", "#66c2a5", "#fc8d61", "#8da0cb"), each = n^2), alpha = 0.5),
  type = "mesh3d", 
  flatshading=TRUE) %>% 
add_trace(x = normAlphas[,1], y = normAlphas[,2], z = normAlphas[,3], 
     type = 'scatter3d', mode = 'markers',
     marker = list(size = 2, color = 'black'),
     opacity = 1, line = list(width = 1, color = 'blue', reverscale = FALSE)) %>% 
layout(showlegend=FALSE, scene = scene)






