source("R/ISbuildImproved.R")
source("R/ISgraph.R")
source("figures/circleCones.R")
source("figures/sphereCones.R")
source("figures/InformationField.R")
# library(RColorBrewer)
# 
# n <- 3
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# pie(rep(1,n), col=sample(col_vector, n))
# sample(col_vector, 2)

## Colors for species in the IS
color1 <- "#80B1D3"
color2 <- "#B2DF8A"
colors2d <- c(color1, color2)


## First 

# Alpha values
a <- as.data.frame(matrix(c(2,1.8),1,2))
# Gamma matrix
g <- matrix(c(-1,0.2,0.3,-1),2,2)
# Build the IS
IS <- ISbuildThird(a,g)
gr <- ISgraph(IS, 1:2)
# Draw the IS in two different ways:
pdf("figures/IS-1a.pdf")
ISgraphDrawLabels(IS,gr, ISgraphLayout(IS, gr, "tree"))
dev.off()
pdf("figures/IS-1b.pdf")
ISgraphDrawPie(IS, gr, ISgraphLayout(IS, gr, "tree"), colors2d)
dev.off()

# Plot the Lyapunov function
plotLyapFunc(IS, 120)

pdf("figures/Space-1.pdf")
# Create the canvas with appropriate limits
limX <- max(IS$points[,1])*1.3 # Limit for x
limY <- max(IS$points[,2])*1.3 # Limit for y
plot(0, 0, xlim = c(0, limX), ylim = c(0, limY),
     type = 'n',xlab = "u1", ylab = "u2",
     geom = "blank")
# Plot vectors showing the direction of each point
for(x in seq(0,limX,limX/15)){
  for(y in seq(0,limY,limY/15)){
    dir <- getDirection(a,g,c(x,y),.01)
    Arrows(
      x0 = x, y0 = y,
      x1 = dir[1], y1 = dir[2],
      arr.length = 0.2, arr.type = 'triangle'
    )
  }
}
pointsCols <- c("#F0027F", "#33A02C", "#A65628", "#8DA0CB", "#999999")
## Show the trayectories for some initial points
initPoints <- data.frame(y1 = c(0.01, 0.2,0.01,3, 0.3),  # u1 values
                         y2 = c(0.5, 0.01,2.3,0.25, 0.3),   # u2 values
                         cols = pointsCols) # Colors

# Solve gLV equations and draw trajectories
pars <- c(a1 = a[1,1], a2 = a[1,2], g12 = g[1,2], g21 = g[2,1])
times <- seq(0,1000,length.out = 2)
for(p in 1:nrow(initPoints)){
  y0 <- c(y1 = initPoints[p,1], y2 = initPoints[p,2])
  sol2 <- ode(y = y0, times = seq(0,20,length.out = 10000),
              parms = pars, func = LotkaVolterra, atol = 1e-12)
  drawLineArrow(sol2, initPoints[p,3], 0.5)
  points(sol2[1,'y1'], sol2[1,'y2'], pch = 19, col = initPoints[p,3], cex=1.5)
}
# Finally, plot the IS nodes
points(IS$points[,1], IS$points[,2], cex = 2, pch = 19, col = "black")
dev.off()



## Second

# Alpha values
a <- as.data.frame(matrix(c(2,.5),1,2))
# Gamma matrix
g <- matrix(c(-1,-0.4,0.2,-1),2,2)
# Build the IS
IS <- ISbuildThird(a,g)
gr <- ISgraph(IS, 1:2)
# Draw the IS in two different ways:
pdf("figures/IS-2a.pdf")
ISgraphDrawLabels(IS,gr, ISgraphLayout(IS, gr, "tree"))
dev.off()
pdf("figures/IS-2b.pdf")
ISgraphDrawPie(IS, gr, ISgraphLayout(IS, gr, "tree"), colors2d)
dev.off()

# Plot the Lyapunov function
plotLyapFunc(IS, 120)

pdf("figures/Space-2.pdf")
# Create the canvas with appropriate limits
limX <- max(IS$points[,1])*1.3 # Limit for x
limY <- max(IS$points[,2])*1.3 # Limit for y
plot(0, 0, xlim = c(0, limX), ylim = c(0, limY),
     type = 'n',xlab = "u1", ylab = "u2",
     geom = "blank")
# Plot vectors showing the direction of each point
for(x in seq(0,limX,limX/15)){
  for(y in seq(0,limY,limY/15)){
    dir <- getDirection(a,g,c(x,y),.01)
    Arrows(
      x0 = x, y0 = y,
      x1 = dir[1], y1 = dir[2],
      arr.length = 0.2, arr.type = 'triangle'
    )
  }
}
pointsCols <- c("#F0027F", "#33A02C", "#A65628", "#8DA0CB", "#999999")
## Show the trayectories for some initial points
initPoints <- data.frame(y1 = c(0.1, 0.01, 2.5, 0.01, 0.001),  # u1 values
                         y2 = c(0.001, 0.6, 0.6, 0.1, 0.3),   # u2 values
                         cols = pointsCols) # Colors

# Solve gLV equations and draw trajectories
pars <- c(a1 = a[1,1], a2 = a[1,2], g12 = g[1,2], g21 = g[2,1])
times <- seq(0,1000,length.out = 2)
for(p in 1:nrow(initPoints)){
  y0 <- c(y1 = initPoints[p,1], y2 = initPoints[p,2])
  sol2 <- ode(y = y0, times = seq(0,20,length.out = 10000),
              parms = pars, func = LotkaVolterra, atol = 1e-12)
  drawLineArrow(sol2, initPoints[p,3], 0.5)
  points(sol2[1,'y1'], sol2[1,'y2'], pch = 19, col = initPoints[p,3], cex=1.5)
}
# Finally, plot the IS nodes
points(IS$points[,1], IS$points[,2], cex = 2, pch = 19, col = "black")
dev.off()


