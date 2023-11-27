
color1 <- "#33A2FF"
color2 <- "#FFC433"
colors12 <- c(color1, color2)



# Draw the circle
drawCircleCones2 <- function(gamma,pointA){
  cone <- setMaxCone2D(gamma)
  # Colors can be changed here
  # The first list contains the colors of the cones
  # The second list, "light" versions of the first three colors
  colors1 <- c("#80FF00E0","#0080FFE0","#FF8000E0","#A0A0A0E0")
  colors2 <- c("#80FF0080","#0080FF80","#FF800080")
  plot(c(-1.5, 1.5), c(-1.5, 1.5), type = "n", axes = FALSE, ann = FALSE, asp = 1)
  coneLabelsSize <- matrix(0, ncol=2, nrow=6)
  # cone 00 (always 3rd quadrant)
  draw.sector(180, 270,clock.wise = FALSE, col = colors1[4], border = NA)
  coneLabelsSize[1,] <- c(180, 270)
  iss <- data.frame(180,270,"00",stringsAsFactors = FALSE)
  if(cone[1]>180){
    # Cone 10
    draw.sector(270, cone[1],clock.wise = FALSE, col = colors1[3], border = NA)
    coneLabelsSize[2,] <- c(270,cone[1])
    iss[nrow(iss)+1,] <- c(270,cone[1], "0010")
    draw.sector(cone[1], 360, clock.wise = FALSE, col = colors2[1], border = NA)
    coneLabelsSize[3,] <- c(cone[1],360)
    iss[nrow(iss)+1,] <- c(cone[1],360,"001011")
    startMainCone <- 0
  }
  if(cone[1]<=180){
    draw.sector(270, 360,clock.wise = FALSE, col = colors1[3], border = NA)
    coneLabelsSize[2,] <- c(270,360)
    iss[nrow(iss)+1,] <- c(270,360,"0010")
    draw.sector(360, cone[1], clock.wise = FALSE, col = colors2[3], border = NA)
    coneLabelsSize[3,] <- c(0,cone[1])
    iss[nrow(iss)+1,] <- c(0,cone[1], "000110")
    startMainCone <- cone[1]
  }
  if(cone[2]>90){
    draw.sector(cone[2], 180, clock.wise = FALSE, col = colors1[2], border = NA)
    coneLabelsSize[4,] <- c(cone[2], 180)
    iss[nrow(iss)+1,] <- c(cone[2], 180, "0001")
    draw.sector(90, cone[2], clock.wise = FALSE, col = colors2[1], border = NA)
    coneLabelsSize[5,] <- c(90, cone[2])
    iss[nrow(iss)+1,] <- c(90, cone[2], "000111")
    endMainCone <- 90
  }
  if(cone[2]<=90){
    draw.sector(90, 180, clock.wise = FALSE, col = colors1[2], border = NA)
    coneLabelsSize[4,] <- c(90, 180)
    iss[nrow(iss)+1,] <- c(90, 180, "0001")
    draw.sector(cone[2], 90, clock.wise = FALSE, col = colors2[2], border = NA)
    coneLabelsSize[5,] <- c(cone[2], 90)
    iss[nrow(iss)+1,] <- c(cone[2], 90, "001001")
    endMainCone <- cone[2]
  }
  draw.sector(startMainCone, endMainCone,clock.wise = FALSE, col = colors1[1], border = NA)
  coneLabelsSize[6,] <- c(startMainCone, endMainCone)
  iss[nrow(iss)+1,] <- c(startMainCone, endMainCone,"00011011")
  # Axes
  lines(c(-1.1,1.1), c(0,0), lwd=1)
  lines(c(0,0), c(-1.1,1.1), lwd=1)
  # Axes labels
  text(-1.3,0, expression('a'[1]),cex=1.3)
  text(0,-1.3, expression('a'[2]),cex=1.3)
  for(i in 1:6){
    writeConeSize(coneLabelsSize[i,1], coneLabelsSize[i,2])
    drawIScone(as.double(iss[i,1]), as.double(iss[i,2]), iss[i,3])
  }
  # Point
  p <- normalize(pointA)
  points(p[1], p[2], pch = 16, cex = 1, col = "red")
}

setColorVal <- function(n){
  if(n>=0){
    "blue"
  } else {
    "red"
  }
}

plot2Dsystem <- function(gmm, alp, point){
  # Crea el grafo
  g <- graph.empty() + vertices(letters[1:2])
  g <- add.edges(g, c('a','a', 'a','b', 'b','a', 'b','b'))
  
  # Define las posiciones de los nodos
  layout <- rbind(c(0,0), c(3,0))
  
  # Define los colores y valores de los nodos
  u1 <- point[1]  # Valor para el primer nodo
  u2 <- point[2] # Valor para el segundo nodo
  
  node_values <- c(u1, u2)
  colfunc <- colorRampPalette(brewer.pal(9, "YlGn"))
  node_colors <- colfunc(100)[cut(node_values, breaks = seq(0, 4, length.out = 100))]
  
  param_values <- c(alp[1], gmm[2,1], gmm[1,2], alp[2])
  
  colfunc <- colorRampPalette(brewer.pal(11, "RdBu"))
  for(i in c(1,4)){
    E(g)$color[i] <- colfunc(100)[cut(param_values[i], 
                                      breaks = seq(-3, 3, length.out = 100))]
  }
  for(i in c(2,3)){
    E(g)$color[i] <- colfunc(100)[cut(param_values[i], 
                                      breaks = seq(-0.5, 0.5, length.out = 100))]
  }
  
  V(g)$color <- node_colors
  V(g)$label <- node_values
  E(g)$width <- 4
  
  plot(g, layout = layout, 
       edge.arrow.size = 1, vertex.size = 40, edge.loop.angle=c(-pi/2,0, 0, -pi/2), 
       edge.curved = 0.4, edge.width=E(g)$width, margin=c(.1,0,-1.5,0))
  if(alp[1] != 0){
    text(-1,-0.35, alp[1], cex=1.3,col=setColorVal(alp[1]))
  }
  if(alp[2] != 0){
    text(1,-0.35, alp[2], cex=1.3,col=setColorVal(alp[2]))
  }
  if(gmm[2,1] != 0){
    text(0.3,-0.7, gmm[2,1], cex=1.3,col=setColorVal(gmm[2,1]))
  }
  if(gmm[1,2] != 0){
    text(-0.3,-1.3, gmm[1,2], cex=1.3,col=setColorVal(gmm[1,2]))
  }
  
}

evolGraphic <- function(g,a,initP,time,checkVecISp){
  
  IS <- ISbuild(a,g)
  u10 <- initP[1]
  u20 <- initP[2]
  limX <- max(IS$points[,1])*1.3 # Limit for x
  limY <- max(IS$points[,2])*1.3 # Limit for y
  
  pars <- c(a1 = a[1], a2 = a[2], g12 = g[1,2], g21 = g[2,1])
  times <- seq(0,time*3,length.out = 2)
  
  y0 <- c(y1 = u10, y2 = u20)
  sol2 <- ode(y = y0, times = seq(0,time,length.out = 10000), parms = pars, func = LotkaVolterra, atol = 1e-12)
  
  
  plot1 <- ggplot2::qplot(sol2[, "time"], sol2[, "y1"], geom = "line", xlab = "Time", ylab = "u1", colour = I(color1), size = I(1))
  plot2 <- ggplot2::qplot(sol2[, "time"], sol2[, "y2"], geom = "line", xlab = "Time", ylab = "u2", colour = I(color2), size = I(1))
  plot3 <- ggplot2::qplot(0, 0, xlim = c(0, limX), ylim = c(0, limY), geom = "blank", xlab = "u1", ylab = "u2")
  
  if(checkVecISp[1]){
    
    arrow_data <- data.frame(
      x = double(),
      y = double(),
      xend = double(),
      yend = double()
    )
    
    for(x in seq(0,limX,limX/20)){
      for(y in seq(0,limY,limY/20)){
        dir <- getDirection(a,g,c(x,y),.01)
        arrow_data <- rbind(arrow_data, data.frame(
          x = x, 
          y = y, 
          xend = dir[1], 
          yend = dir[2]
        ))
      }
    }
    
    plot3 <- plot3 + geom_segment(
      data = arrow_data, 
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.2, "cm"), angle=15, type = "closed")
    )
  }
  
  df_sol2 <- data.frame(y1 = sol2[,'y1'], y2 = sol2[,'y2'])
  plot3 <- plot3 + geom_path(data = df_sol2, 
                             aes(x = y1, y = y2), 
                             linetype = 1, color = "green", linewidth = 1)
  
  
  plot3 <- plot3 + geom_point(aes(x = sol2[1,'y1'], y = sol2[1,'y2']), 
                              color = "green", size = 3, shape = 19)
  
  if(checkVecISp[2]){
    df_IS_points <- data.frame(x = IS$points[,1], y = IS$points[,2])
    plot3 <- plot3 + geom_point(data = df_IS_points, aes(x = x, y = y), color = "black", size = 4, shape = 19)
  }
  grid.arrange(plot1, plot2, plot3, ncol = 2, nrow = 2, layout_matrix = rbind(c(1, 2), c(3, 3)), heights = c(1, 2))
}


## Display the IF

# Function to plot the Information Field of a given
# IS with 2 species. The argument 'np' sets
# the divisions of the grid. The higher value, the
# best resolution but taking more time.
plotLyapFunc2 <- function(IS, g, a, stateSys, time, np, checkIF){
  # Set x and y points:
  # Values are from 0 to 1.3 times the maximum in the IS points. The
  # coordinates of IS points are added for better visualization.
  x <- unique(sort(c(seq(0, max(c(1,IS$points[,1]))*1.3, length.out = np),
                     array(IS$points[,1]))))
  y <- unique(sort(c(seq(0, max(c(1,IS$points[,2]))*1.3, length.out = np),
                     array(IS$points[,2]))))
  
  # Value of the Lyapunov function for all (x,y):
  zetas <- matrix(0,length(x), length(y))
  for(i in 1:length(x)){
    for(j in 1:length(y)){
      zetas[i,j] <- lyapunovFunc(IS,c(x[i],y[j]))
    }
  }
  # Transpose for plotting:
  z <- t(zetas)
  
  # IS points:
  zetas <- rep(0, nrow(IS$points))
  for(r in 1:nrow(IS$points)){
    zetas[r] <- lyapunovFunc(IS,c(IS$points[r,1], IS$points[r,2]))
  }
  
  # Figure with the plot of the Lyapunov function:
  kd <- with(MASS::geyser, MASS::kde2d(duration, waiting, n = 50))
  fig <- plot_ly(x = x, y = y, z = z,
                 showscale = FALSE,
                 showticklabels = FALSE,
                 colorscale= 'Viridis',
                 reversescale = TRUE,
                 showgrid = TRUE
  ) %>% add_surface()
  # Add marks for the IS points:
  fig <- fig %>%
    add_trace(
      type = "scatter3d",
      mode = "markers",
      x = array(IS$points[,1]),
      y = array(IS$points[,2]),
      z = zetas,
      showlegend = FALSE,
      marker = list(color = 'black', size = 3)
    )
  fig <- fig %>%
    layout(
      scene = list(
        xaxis = list(title = "u1", showticklabels = FALSE),
        yaxis = list(title = "u2", showticklabels = FALSE),
        zaxis = list(title = "", showticklabels = FALSE)
      )
    )
  
  
  if(checkIF){
    ## Add the lines representing the evolution
    u10 <- stateSys[1]
    u20 <- stateSys[2]
    
    pars <- c(a1 = a[1], a2 = a[2], g12 = g[1,2], g21 = g[2,1])
    times <- seq(0,time*3,length.out = 2)
    
    y0 <- c(y1 = u10, y2 = u20)
    sol2 <- ode(y = y0, times = seq(0,time,length.out = 10000), parms = pars, func = LotkaVolterra, atol = 1e-12)
    
    zetas <- mapply(function(x, y) lyapunovFunc(IS, c(x, y)), sol2[, "y1"], sol2[, "y2"])
    
    fig <- fig %>%
      add_trace(
        type = "scatter3d",
        mode = "lines",
        x = sol2[, "y1"],
        y = sol2[, "y2"],
        z = rep(0, length(sol2[, "y1"])),
        line = list(color = "green", width = 4),
        showlegend = FALSE
      )
    
    fig <- fig %>%
      add_trace(
        type = "scatter3d",
        mode = "lines",
        x = sol2[, "y1"],
        y = sol2[, "y2"],
        z = zetas,
        line = list(color = "#FF4633", width = 6),
        showlegend = FALSE
      )
    
    for(i in seq(1, length(sol2[, "y1"]), 200)){
      fig <- fig %>% add_trace( 
        x = rep(sol2[i, "y1"], 2), y = rep(sol2[i, "y2"], 2),
        z = c(zetas[i],0), 
        type = "scatter3d",
        mode = "lines",
        color = I("black"), 
        line = list(color = 'black'), showlegend = F)
    }
  }
  fig
}
