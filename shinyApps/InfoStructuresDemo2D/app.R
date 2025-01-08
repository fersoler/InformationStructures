library(shiny)
library(plotly)
library(gridExtra)
library(geometry)
library(readr)
library(mathjaxr)
source("R/ISbuild.R")
source("R/ISgraph.R")
source("R/figs/circleCones.R")
source("R/figs/InformationField.R")
source("shinyApps/InfoStructuresDemo2D/R/modFuncs.R")


# Define the app UI
ui <- fluidPage(
  title = "Information Structures and Information Fields",
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
  titlePanel("Information Structures and Information Fields"),
  helpText("Modelling the interaction between two species using the gLV model."),
  sidebarPanel(
    sliderInput("g12", strong("Influence of 2 over 1 ($\\gamma_{12}$)"),
                min = -0.45, max = 0.45, value = 0.2, step = 0.01
    ),
    sliderInput("g21", strong("Influence of 1 over 2 ($\\gamma_{21}$)"),
                min = -0.45, max = 0.45, value = 0.3, step = 0.01
    ),
    sliderInput("a1", strong("Intrinsic grow rate of 1 ($\\alpha_1$)"),
                min = -3.00, max = 3.00, value = 1, step = 0.1,
    ),
    sliderInput("a2", strong("Intrinsic grow rate of 2 ($\\alpha_2$)"),
                min = -3, max = 3, value = 1.4, step = 0.1
    ),
    sliderInput("u10", strong("Initial abundance of 1 ($u_1(0)$)"),
                min = 0, max = 4, value = 1, step = 0.001
    ),
    sliderInput("u20", strong("Initial abundance of 2 ($u_2(0)$)"),
                min = 0, max = 4, value = 1.6, step = 0.001
    ),
    width = 3
  ),
  mainPanel(
    tabsetPanel(
      id = "tabs",
      tabPanel(
        title = "System",
        icon = icon("arrows-to-circle"),
        helpText(HTML(paste("This graph represents a system of two components.", 
                            "Each component is given by a variable $u_i$", 
                            "whose value appears inside the corresponding node.", 
                            "There are parameters determining how a component", 
                            "influences itself ($\\alpha_i$) and the interactions", 
                            "between components ($\\gamma_{ij}$).", 
                            "The parameters can be modified in the left settings.", 
                            "The system behaves by means of the", 
                            "generalized Lotka-Volterra model (below)", 
                            "which determine the derivative of each", 
                            "component in time."))),
        column(width=12, 
               plotOutput("plotSystem"),uiOutput("dynamicText1"))
      ),
      tabPanel(
        title = "Evolution",
        icon = icon("chart-line"),
        helpText(HTML(paste("In this tab you can see the evolution of the system",
                       "over time with the selected parameters, starting",
                       "at the defined starting point. In the <em>Settings</em> tab",
                       "you can specify the value of <em>Time</em>. You can", 
                       "also display points corresponding to IS nodes and", 
                       "vectors indicating the directions at different",
                       "points of the phase space. "))),
        plotOutput("compPlotEvol")
      ),
      tabPanel(
        title = "Info Struc",
        icon = icon("diagram-project"),
        helpText(paste("These graphs represent the Information Structure", 
                       "containing all stationary communities. Edges of the",
                       "graph represent solutions going from one community",
                       "to another one.",
                       "Colors represent relative abundances.")),
        fluidRow(
          class = "align-items-center",
          column(5, plotOutput("plotIS")),
          column(5, plotOutput("plotIS2")),
          column(2, tableOutput('ISpoints'))
        ),
      ),
      tabPanel(
        title = "Info Field",
        icon = icon("mountain-sun"),
        helpText(paste("This is the informational field,",
                       "which represents the energy level of each point in", 
                       "the phase space. It is built by means of a Lyapunov function.",
                       "In the setting tab you can select the option to display",
                       "lines representing the evolution of the system.")),
        column(12, plotlyOutput("plotIF"))
      ),
      tabPanel(
        title = "Cones",
        icon = icon("chart-pie"),
        helpText(paste("The circle represents the biodiversity cones of",
                       "the system, which depend on the $\\gamma_{ij}$ values.",
                       "The red dot indicates where is the system,",
                       "depending on the $\\alpha_i$ values.",
                       "By modifying the $\\alpha_i$, the system can",
                       "change of cone and the IS is modified.")),
        column(6, plotOutput("plotCones")),
        column(6, plotOutput("plotISb"))
      ),
      tabPanel(
        title = "Settings",
        icon = icon("sliders-h"),
        sliderInput("varTime", strong("Time"),
                    min = 5, max = 100, value = 20, step = 5
        ),
        helpText("Set the Time to display the system evolution"),
        checkboxInput("checkVec", "Show the directions in the phase space", value = FALSE),
        checkboxInput("checkISp", "Mark the IS points in the phase space", value = FALSE),
        sliderInput("ifRes", strong("Resolution of the Info Field"),
                    min = 30, max = 100, value = 40, step = 10
        ),
        helpText("Set the resolution of the mesh representing the IF"),
        checkboxInput("checkIF", "Show the evolution in the IF", value = FALSE)
      )
    )
  )
)

# define the app logic:
server <- function(input, output) {
  
  output$compPlotEvol <- renderPlot({
    
    a1 <- input$a1
    a2 <- input$a2
    g12 <- input$g12
    g21 <- input$g21
    u1 <- input$u10
    u2 <- input$u20
    
    g <- matrix(c(-1,g21,g12,-1),2,2)
    a <-  c(a1,a2)
    
    return(evolGraphic(g,a,c(u1,u2), 
                       input$varTime, 
                       c(input$checkVec, input$checkISp)))
  })
  
  output$plotSystem <- renderPlot({
    
    a1 <- input$a1
    a2 <- input$a2
    g12 <- input$g12
    g21 <- input$g21
    u1 <- input$u10
    u2 <- input$u20
    
    gmm <- matrix(c(0,g21,g12,0),2,2)
    
    plot2Dsystem(gmm,c(a1,a2),c(u1,u2))
    
  })
  
  output$plotCones <- renderPlot({
    
    a1 <- input$a1
    a2 <- input$a2
    g12 <- input$g12
    g21 <- input$g21
    
    # Gamma matrix
    g <- matrix(c(-1,g21,g12,-1),2,2)
    
    drawCircleCones2(g,c(a1,a2))
    
  })
  
  output$plotIF <- renderPlotly({
    
    a1 <- input$a1
    a2 <- input$a2
    g12 <- input$g12
    g21 <- input$g21
    
    u1 <- input$u10
    u2 <- input$u20
    
    a <- c(a1,a2)
    # Gamma matrix
    g <- matrix(c(-1,g21,g12,-1),2,2)
    # Build the IS
    IS <- ISbuild(a,g)
    # Plot the Lyapunov function
    plotLyapFunc2(IS, g, a, c(u1,u2), input$varTime, input$ifRes, input$checkIF)
  
  })
  
  
  output$plotIS <- renderPlot({
    
    a1 <- input$a1
    a2 <- input$a2
    g12 <- input$g12
    g21 <- input$g21
    
    a <- c(a1,a2)
    # Gamma matrix
    g <- matrix(c(-1,g21,g12,-1),2,2)
    # Build the IS
    IS <- ISbuild(a,g)
    gr <- ISgraph(IS, 1:2)
    ISgraphDrawPie(IS, gr, ISgraphLayout(IS, gr, "tree"), colors12)
  })
  
  
  output$plotISb <- renderPlot({
    
    a1 <- input$a1
    a2 <- input$a2
    g12 <- input$g12
    g21 <- input$g21
    
    a <- c(a1,a2)
    # Gamma matrix
    g <- matrix(c(-1,g21,g12,-1),2,2)
    # Build the IS
    IS <- ISbuild(a,g)
    gr <- ISgraph(IS, 1:2)
    ISgraphDrawPie(IS, gr, ISgraphLayout(IS, gr, "tree"), colors12)
  })
  
  output$plotIS2 <- renderPlot({
    
    a1 <- input$a1
    a2 <- input$a2
    g12 <- input$g12
    g21 <- input$g21
    
    a <- c(a1,a2)
    # Gamma matrix
    g <- matrix(c(-1,g21,g12,-1),2,2)
    # Build the IS
    IS <- ISbuild(a,g)
    gr <- ISgraph(IS, 1:2)
    ISgraphDrawLabels(IS,gr, ISgraphLayout(IS, gr, "tree"))
  })
  
  output$ISpoints <- renderTable({
    
    a1 <- input$a1
    a2 <- input$a2
    g12 <- input$g12
    g21 <- input$g21
    
    a <- c(a1,a2)
    # Gamma matrix
    g <- matrix(c(-1,g21,g12,-1),2,2)
    # Build the IS
    IS <- ISbuild(a,g)
  
    pp <- IS$points
    
    pp <- rbind(pp[-IS$gassInd, ], pp[IS$gassInd, ])
    
    pp <- round(pp,3)
    colnames(pp) <- c("u1","u2")
    
    pp
  
  })
  
  
  output$dynamicText1 <- renderUI({
    
    a1 <- input$a1
    a2 <- input$a2
    g12 <- input$g12
    g21 <- input$g21
    u1 <- input$u10
    u2 <- input$u20
    
    withMathJax(paste(
      "$$",
      "u_1' = u_1 \\cdot (\\alpha_1 - u_1 + \\gamma_{12}\\cdot u_2) = ", 
      round(u1*(a1 - u1 + u2*g12),3), "$$",
      "$$",
      "u_2' = u_2 \\cdot (\\alpha_2 - u_2 + \\gamma_{21}\\cdot u_1) = ", 
      round(u2*(a2 - u2 + u1*g21),3), "$$"))
  })
  
  
}



# run the app:
shinyApp(ui = ui, server = server)

