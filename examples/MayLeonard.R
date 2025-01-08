#source("R/ISbuild.R")
#source("R/ISgraph.R")
library(plotly)
#library(geometry)
library(deSolve)

# r values
r <- c(1,1,1)
# alpha, beta
alpha <- 2
beta  <- 2
# Gamma matrix
(A <- matrix(c(-1, -beta, -alpha, -alpha, -1, -beta, -beta, -alpha, -1),3,3))

# Build the IS
#(IS <- ISbuild(r, A))
#gr <- ISgraph(IS, 1:3)
# Draw the IS in two different ways:
#ISgraphDrawLabels(IS,gr, ISgraphLayout(IS, gr, "3Dcube"))

points <- rbind(c(0,0,0),
                c(1,0,0),
                c(0,1,0),
                c(0,0,1),
                (1/(1-alpha*beta)) * c(1-alpha, 1-beta, 0),
                (1/(1-alpha*beta)) * c(0, 1-alpha, 1-beta),
                (1/(1-alpha*beta)) * c(1-beta, 0, 1-alpha),
                (1/(1+alpha+beta)) * c(1, 1, 1)
                )


### Iniciar gráfico en blanco
fig <- plot_ly() %>% layout(
  scene = list(
    xaxis = list(
      showgrid = FALSE,
      zeroline = TRUE,     # Línea del eje en el origen
      showline = FALSE,    # No mostrar línea adicional
      range = c(0, 1.2),     # Longitud del eje desde 0
      title = list(
        text = "sp1",
        font = list(size = 12),
        xanchor = "left"   # Posicionar etiqueta cerca del eje
      )
    ),
    yaxis = list(
      showgrid = FALSE,
      zeroline = TRUE,
      showline = FALSE,
      range = c(0, 1.2),
      title = list(
        text = "sp2",
        font = list(size = 12),
        yanchor = "bottom"
      )
    ),
    zaxis = list(
      showgrid = FALSE,
      zeroline = TRUE,
      showline = FALSE,
      range = c(0, 1.2),
      title = list(
        text = "sp3",
        font = list(size = 12),
        zanchor = "bottom"
      )
    ),
    aspectratio = list(x = 1, y = 1, z = 1)  # Proporciones iguales
  )
)

# Add points
fig <- fig %>%
  add_trace(x = points[,1], 
            y = points[,2], 
            z = points[,3],
            type = 'scatter3d', mode = 'markers',
            showlegend = FALSE,
            name = 'IS points',
            hovertemplate = paste('<b>sp1</b>: %{x:.2f}',
                                  '<br><b>sp2</b>: %{y:.2f}',
                                  '<br><b>sp3</b>: %{z:.2f}'),
            marker = list(size = 4, color = "black", opacity = 1, symbol = 'circle'))

# Add triangle

fig <- fig %>%
  add_trace(
    type = 'mesh3d',
    x = c(1, 0, 0),  
    y = c(0, 1, 0),  
    z = c(0, 0, 1),  
    i = c(0),  
    j = c(1),  
    k = c(2),  
    facecolor = rep("lightgrey", 1),  # Asigna un color fijo a la cara del triángulo
    #flatshading = TRUE,
    opacity = 0.5,
    lighting = list(),
    hoverinfo = 'none'
  )




initpoint <- c(1,1.1,1.2)
time <- 200

# Equations of the gLV system
gLV <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    list(state * (r + A %*% state))
  })
}
pars <- c(r = r, A = A)
#y0 <- c(y1 = initpoint[1], y2 = initpoint[2], y3 = initpoint[3])
sol2 <- ode(y = initpoint, times = seq(0, time, length.out = 10000),
            parms = pars, func = gLV, atol = 1e-12)


fig <- fig %>% 
  add_trace(
    type = "scatter3d",
    mode = "lines",
    x = sol2[,2],
    y = sol2[,3],
    z = sol2[,4],
    #z = zetas+shift[3],
    line = list(color = "orange", width = 2),
    showlegend = FALSE,
    name = paste('Init sp1:', round(initpoint[1],2),
                 '<br>Init sp2:', round(initpoint[2],2),
                 '<br>Init sp3:', round(initpoint[3],2)),
    hovertemplate = paste('<b>sp1</b>: %{x:.2f}',
                          '<br><b>sp2</b>: %{y:.2f}',
                          '<br><b>sp3</b>: %{z:.2f}')
  )

fig


########
# Diagram with alpha/beta positions
par(pty = "s")
plot(0, 0,xlim = c(0, 3), ylim = c(0, 3),
  type = 'n',
  xlab = expression(alpha),  # Etiqueta del eje x como letra griega α
  ylab = expression(beta),   # Etiqueta del eje y como letra griega β
  asp = 1,                   # Aspect ratio igual para x e y
  xaxs = "i",                # Evita que se extienda el eje x
  yaxs = "i" 
)
segments(1, 0, 1, 3, col = "blue", lwd = 2)
segments(0, 1, 3, 1, col = "blue", lwd = 2)
segments(2, 0, 0, 2, col = "blue", lwd = 2)
points(alpha, beta, pch = 19, col = "red", cex = 1.5)
