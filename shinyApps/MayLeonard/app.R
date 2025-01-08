#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(plotly)
library(deSolve)
library(mathjaxr)
library(shiny)

###############################################################
# general functions for the May-Leonard app
###############################################################

# Get equilibrium points given alpha and beta
getPoints <- function(alpha, beta){
  points <- rbind(c(0,0,0),
                  c(1,0,0),
                  c(0,1,0),
                  c(0,0,1),
                  (1/(1+alpha+beta)) * c(1, 1, 1)
  )
  if(alpha<1 & beta < 1){
    points <- rbind(points,
                    (1/(1-alpha*beta)) * c(1-alpha, 1-beta, 0),
                    (1/(1-alpha*beta)) * c(0, 1-alpha, 1-beta),
                    (1/(1-alpha*beta)) * c(1-beta, 0, 1-alpha))
  }
  points
}

# gLV model
gLV <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    list(state * (r + A %*% state))
  })
}

# Draw the figure with the relevant points
createAttFig <- function(points){
  # start blank plotly fig
  fig <- plot_ly() %>% layout(
    scene = list(
      xaxis = list(
        showgrid = FALSE,
        zeroline = TRUE,     
        showline = FALSE,    
        range = c(0, 1.2),   
        title = list(
          text = "sp1",
          font = list(size = 12),
          xanchor = "left"   
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
      aspectratio = list(x = 1, y = 1, z = 1) 
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
      facecolor = rep("lightgrey", 1), 
      #flatshading = TRUE,
      opacity = 0.5,
      lighting = list(),
      hoverinfo = 'none'
    )
  fig  
}

# add solution
addSolution <- function(fig, initpoint, time, alpha, beta){
  mA <- matrix(c(-1, -beta, -alpha, -alpha, -1, -beta, -beta, -alpha, -1),3,3)
  pars <- list(r = c(1,1,1), A = mA)
  solLV <- ode(y = initpoint, times = seq(0, time, length.out = 10000),
               parms = pars, func = gLV, atol = 1e-12)
  fig <- fig %>%
    add_trace(x = initpoint[1], 
              y = initpoint[2], 
              z = initpoint[3], 
              type = 'scatter3d', mode = 'markers',
              showlegend = FALSE,
              name = 'Init point',
              hovertemplate = paste('<b>sp1</b>: %{x:.2f}',
                                    '<br><b>sp2</b>: %{y:.2f}',
                                    '<br><b>sp3</b>: %{z:.2f}'),
              marker = list(size = 2, color = "orange", opacity = 1, symbol = 'circle'))
  fig <- fig %>% 
    add_trace(
      type = "scatter3d",
      mode = "lines",
      x = solLV[,2],
      y = solLV[,3],
      z = solLV[,4],
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
}




##################################################################
# shiny app
##################################################################


# Define UI for application that draws a histogram
ui <- fluidPage(
  
  tags$head(
    tags$style(
      HTML("
        .align-items-center {
          display: flex;
          align-items: center;
        }
      ")
    )
  ),
  withMathJax(),
  tags$div(HTML("
            <script type='text/x-mathjax-config' >
            MathJax.Hub.Config({
            tex2jax: {inlineMath: [['$','$']]}
            });
            </script >
            ")),
  
  # Application title
  titlePanel("The May-Leonard model"),
  
  # Sidebar 
  sidebarLayout(
    sidebarPanel(
      sliderInput("alpha", strong("Value of $\\alpha$"),
                  min = 0.01, max = 3, value = 1.3, step = 0.01
      ),
      sliderInput("beta", strong("Value of $\\beta$"),
                  min = 0.01, max = 3, value = 0.75, step = 0.01
      ),
      numericInput("sp1", strong("Initial value of 1st species"),
                   min = 0, max = 3, value = .7, step = 0.00001,
      ),
      numericInput("sp2", strong("Initial value of 2nd species"),
                   min = 0, max = 3, value = .6, step = 0.00001,
      ),
      numericInput("sp3", strong("Initial value of 3rd species"),
                   min = 0, max = 3, value = .9, step = 0.00001,
      ),
      sliderInput("time", strong("Time"),
                  min = 100, max = 10000, value = 300, step = 100,
      ),
      width = 3),
    
    mainPanel(
      
      fluidRow(
        column(7, plotlyOutput("plotMayLeonard")),
        column(5, plotOutput("plotParams"))
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$plotParams <- renderPlot({
    
    alpha <- input$alpha
    beta <- input$beta
    
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
    
    
  })
  
  output$plotMayLeonard <- renderPlotly({
    
    r <- c(1,1,1)
    alpha <- input$alpha
    beta <- input$beta
    initpoint <- c(input$sp1, input$sp2, input$sp3)
    time <- input$time
    
    points <- getPoints(alpha, beta)
    fig <- createAttFig(points)
    fig <- addSolution(fig, initpoint, time, alpha, beta)
    fig
    
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
