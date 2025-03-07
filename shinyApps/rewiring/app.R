source("R/ISbuild.R")
source("R/ISgraph.R")
source("R/rewiring.R")
source("R/figs/circleCones.R")
source("R/figs/sphereCones.R")
source("R/figs/sphereSubCones.R")
source("R/figs/InformationField.R")


library(shiny)


# Function to draw 2d cones with a point representing r_i
drawCircleCones_with_R <- function(a, r){
  drawCircleCones(a)
  point <- r #r/sqrt(sum(r^2))
  points(point[1],point[2], pch = 21, bg = "yellow", cex = 1)
}


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
  titlePanel("Rewiring LV parameters"),

  # Sidebar 
  sidebarLayout(
    sidebarPanel(
      sliderInput("a11", strong("$a_{11}$"),
                  min = -3, max = -.1, value = -1, step = 0.01
      ),
      sliderInput("a22", strong("$a_{22}$"),
                  min = -3, max = -.1, value = -1, step = 0.01
      ),
      sliderInput("a12", strong("$a_{12}$"),
                  min = -.99, max = .99, value = .4, step = 0.01
      ),
      sliderInput("a21", strong("$a_{21}$"),
                  min = -.99, max = .99, value = .4, step = 0.01
      ),
      sliderInput("r1", strong("$r_{1}$"),
                  min = -3, max = 3, value = 1, step = 0.01
      ),
      sliderInput("r2", strong("$r_{2}$"),
                  min = -3, max = 3, value = 1, step = 0.01
      ),
      width = 2),
    
    mainPanel(
      fluidRow(#style = "margin: 0px; padding: 0px;",
        column(6, 
               #style = "margin: 0px; padding: 0px;",
               div("Initial system", style = "text-align: center; margin: 0;"),
               fluidRow(#style = "margin: 0px; padding: 0px;",
                        column(5,
                   plotOutput("plotInitial")),
                   column(4,plotOutput("isInitial")),
                   column(3,uiOutput("plotParamsInit"))),
               div("Increase inter-cooperation", style = "text-align: center; margin: 0;"),
               fluidRow(#style = "margin: 0px; padding: 0px;",
                        column(5,
                               plotOutput("plotInter")),
                        column(4,plotOutput("isInter")),
                        column(3,uiOutput("plotParamsInter")))), 
        column(6, 
               #style = "margin: 0px; padding: 0px;",
               div("Increase intra-competition", style = "text-align: center; margin: 0;"),
               fluidRow(#style = "margin: 0px; padding: 0px;",
                        column(5,
                               plotOutput("plotIntra")),
                        column(4,plotOutput("isIntra")), 
                        column(3,uiOutput("plotParamsIntra"))), 
               div("Both effects together", style = "text-align: center; margin: 0;"),
               fluidRow(#style = "margin: 0px; padding: 0px;",
                        column(5,
                               plotOutput("plotBoth")),
                        column(4,plotOutput("isBoth")),
                        column(3,uiOutput("plotParamsBoth")))),
        
      )
    )
))

# Define server logic required to draw a histogram
server <- function(input, output) {

  output$plotInitial <- renderPlot({
    
    par(mar = c(0, 0, 0, 0))
    a <- matrix(c(input$a11, input$a21, input$a12, input$a22),2,2)
    r <- c(input$r1, input$r2)
    
    drawCircleCones_with_R(a, r/2)
    
  })
  
  
  output$isInitial <- renderPlot({
    
    par(mar = c(0, 0, 0, 0))
    a <- matrix(c(input$a11, input$a21, input$a12, input$a22),2,2)
    r <- c(input$r1, input$r2)
    
    IS <- ISbuild(r,a)
    gr <- ISgraph(IS, 1:2)
    ISgraphDrawPie(IS, gr, ISgraphLayout(IS, gr, "tree"), c("#FF8000E0","#0080FFE0"))
  })
  
  
  output$plotIntra <- renderPlot({
    par(mar = c(0, 0, 0, 0))
    a <- matrix(c(input$a11, input$a21, input$a12, input$a22),2,2)
    r <- c(input$r1, input$r2)
    
    a2 <- rewiringIncreaseIntraComp(a,r)
    drawCircleCones_with_R(a2, r/2)

  })
  
  output$isIntra <- renderPlot({
    
    par(mar = c(0, 0, 0, 0))
    a <- matrix(c(input$a11, input$a21, input$a12, input$a22),2,2)
    r <- c(input$r1, input$r2)
    a2 <- rewiringIncreaseIntraComp(a,r)
    IS <- ISbuild(r,a2)
    gr <- ISgraph(IS, 1:2)
    ISgraphDrawPie(IS, gr, ISgraphLayout(IS, gr, "tree"), c("#FF8000E0","#0080FFE0"))
  })
  
  output$plotInter <- renderPlot({
    par(mar = c(0, 0, 0, 0))
    a <- matrix(c(input$a11, input$a21, input$a12, input$a22),2,2)
    r <- c(input$r1, input$r2)
    
    a2 <- keepDiagDomRows(a, rewiringIncreaseCoop(a,r))
    drawCircleCones_with_R(a2, r/2)

  })

  output$isInter <- renderPlot({
    
    par(mar = c(0, 0, 0, 0))
    a <- matrix(c(input$a11, input$a21, input$a12, input$a22),2,2)
    r <- c(input$r1, input$r2)
    a2 <- keepDiagDomRows(a, rewiringIncreaseCoop(a,r))
    IS <- ISbuild(r,a2)
    gr <- ISgraph(IS, 1:2)
    ISgraphDrawPie(IS, gr, ISgraphLayout(IS, gr, "tree"), c("#FF8000E0","#0080FFE0"))
  })
  
    
  output$plotBoth <- renderPlot({
    par(mar = c(0, 0, 0, 0))
    a <- matrix(c(input$a11, input$a21, input$a12, input$a22),2,2)
    r <- c(input$r1, input$r2)
    a2 <- rewiringIncreaseIntraComp(a,r)
    a3 <- keepDiagDomRows(a2, rewiringIncreaseCoop(a2,r))
    drawCircleCones_with_R(a3, r/2)

  })
  
  output$isBoth <- renderPlot({
    par(mar = c(0, 0, 0, 0))
    a <- matrix(c(input$a11, input$a21, input$a12, input$a22),2,2)
    r <- c(input$r1, input$r2)
    a2 <- rewiringIncreaseIntraComp(a,r)
    a3 <- keepDiagDomRows(a2, rewiringIncreaseCoop(a2,r))
    IS <- ISbuild(r,a3)
    gr <- ISgraph(IS, 1:2)
    ISgraphDrawPie(IS, gr, ISgraphLayout(IS, gr, "tree"), c("#FF8000E0","#0080FFE0"))
  })
  
  output$plotParamsInit <- renderUI({
    par(mar = c(0, 0, 0, 0))
    
    a <- matrix(c(input$a11, input$a21, input$a12, input$a22),2,2)
    a <- round(a,3)
    withMathJax(paste(
      "$$",
      "\\begin{bmatrix}",
      a[1,1], "&", a[1,2], "\\\\",
      a[2,1], "&", a[2,2],
      "\\end{bmatrix}", "$$"
    ))
  })
  
  
  output$plotParamsIntra <- renderUI({
    par(mar = c(0, 0, 0, 0))
    
    a <- matrix(c(input$a11, input$a21, input$a12, input$a22),2,2)
    r <- c(input$r1, input$r2)
    a2 <- rewiringIncreaseIntraComp(a,r)
    a2 <- round(a2,3)
    withMathJax(paste(
      "$$",
      "\\begin{bmatrix}",
      a2[1,1], "&", a2[1,2], "\\\\",
      a2[2,1], "&", a2[2,2],
      "\\end{bmatrix}", "$$"
    ))
  })
  
  output$plotParamsInter <- renderUI({
    par(mar = c(0, 0, 0, 0))
    
    a <- matrix(c(input$a11, input$a21, input$a12, input$a22),2,2)
    r <- c(input$r1, input$r2)
    a2 <- keepDiagDomRows(a, rewiringIncreaseCoop(a,r))
    a2 <- round(a2,3)
    withMathJax(paste(
      "$$",
      "\\begin{bmatrix}",
      a2[1,1], "&", a2[1,2], "\\\\",
      a2[2,1], "&", a2[2,2],
      "\\end{bmatrix}", "$$"
    ))
  })

  output$plotParamsBoth <- renderUI({
    par(mar = c(0, 0, 0, 0))
    a <- matrix(c(input$a11, input$a21, input$a12, input$a22),2,2)
    r <- c(input$r1, input$r2)
    a2 <- rewiringIncreaseIntraComp(a,r)
    a3 <- keepDiagDomRows(a2, rewiringIncreaseCoop(a2,r))
    a3 <- round(a3,3)
    withMathJax(paste(
      "$$",
      "\\begin{bmatrix}",
      a3[1,1], "&", a3[1,2], "\\\\",
      a3[2,1], "&", a3[2,2],
      "\\end{bmatrix}", "$$"
    ))
  })
  
    
}

# Run the application 
shinyApp(ui = ui, server = server)

