library(shiny)
library(plotly)
library(readr)
library(mathjaxr)
source("R/ISbuild.R")
source("R/ISgraph.R")
source("R/ISbuildImproved.R")


gMatrix <- t(matrix(c(-1,-0.1,-0.23,-0.3,-1,-0.12,-0.2,-0.17,-1),3,3))
  

ISgraphDrawPie2 <- function(IS, ISgr, ISlay, colors, maxAbund = 0){
  gr <- ISgr$graph
  isP <- t(IS$points)
  abund <- colSums(isP)
  maxAbundThis <- max(abund)
  for(i in 1:ncol(isP)){
    if(maxAbund == 0){
      V(gr)[i]$size = 20
    }
    else{
      if(maxAbund == -1){
        mm = maxAbundThis
      }else{
        mm = maxAbund
      }
      V(gr)[i]$size = 10+10*(abund[i]/mm)
    }
    if(i==1){
      V(gr)[i]$shape = "circle"
      V(gr)[i]$color = "white"
    }
    else {
      if(sum(as.integer(isP[,i]>0))==1){
        V(gr)[i]$shape = "circle"
        V(gr)[i]$color = colors[isP[,i]>0]
      } else {
        V(gr)[i]$shape = "pie"
        V(gr)[i]$pie = list(10*(isP[,i]))
        V(gr)[i]$pie.color = list(colors)
      }
    }
  }
  list(gr = gr, lay = ISlay, group = list(c(IS$gassInd)))
}


showIS <- function(aVals){
  IS3 <- ISbuildThird(as.data.frame(t(aVals)),gMatrix)
  gr3 <- ISgraph(IS3, 1:3)
  ISgraphDrawPie2(IS3,gr3, ISgraphLayout(IS3, gr3, "3Dcube"), c("#806000", "#002080", "#408000"))
}

# Genera los datos iniciales
load('shinyApps/CompetitiveSphere/datosComp.RData')


# Define la UI de la aplicación
ui <- fluidPage(
  
  titlePanel("The inner structure of biodiversity (competitive version)"),
  
  tags$p("Scroll over the sphere and click a point", style = "text-align:left;"),
  
  hr(),
  
  fluidRow(
    column(width = 6, plotlyOutput("plot1"),uiOutput("dynamicText1")),
    column(width = 6, plotOutput("plot2"),uiOutput("dynamicText2"))
  )
)

# Define la lógica de la aplicación
server <- function(input, output) {

  # Genera el primer plot
  output$plot1 <- renderPlotly({
    plot_ly(
      x = array(spherePoints[1,]),
      y = array(spherePoints[2,]),
      z = array(spherePoints[3,]),
      i = array(sphereTriang[,1]),
      j = array(sphereTriang[,2]),
      k = array(sphereTriang[,3]),
      facecolor = toRGB(sphereColors,alpha=1),
      type = "mesh3d",
      flatshading = TRUE,
      lighting = list(ambient = 0.7, # 0.5
                      diffuse = 0.8,
                      fresnel = 0.1,
                      specular = 0.5,
                      roughness = 0.5,
                      facenormalsepsilon = 0,
                      vertexnormalsepsilon = 0)
    )
  })
  
  # Genera el segundo plot
  output$plot2 <- renderPlot({
    event_data <- event_data("plotly_click")
    req(!is.null(event_data))
    
    x <- event_data$x
    y <- event_data$y
    z <- event_data$z
    
    datos <- showIS(c(x,y,z))

    x_lim <- c(-5, 15)
    y_lim <- c(-10, 10)
    plot(datos$gr,
         layout = (datos$lay/2.5),
         edge.width=1,
         mark.groups=datos$group,
         rescale=FALSE,
         mark.col="orange",
         mark.expand=4,
         mark.shape=-2,
         vertex.label=NA)
    graphics::plot.window(xlim=x_lim, ylim=y_lim)
  })
  
  output$dynamicText1 <- renderUI({
    
    withMathJax(paste(
      "$$",
      "\\Gamma = \\begin{bmatrix}",
      "-1 &", gMatrix[1,2], " & ", gMatrix[1,3], "\\\\",
      gMatrix[2,1], " & -1 & ", gMatrix[2,3], "\\\\",
      gMatrix[3,1], " &", gMatrix[3,2], " & -1",
      "\\end{bmatrix}", "$$"
    ))
  })
  
  output$dynamicText2 <- renderUI({
    
    d <- event_data("plotly_click")
    
    if(is.null(d)){
      x <- 0
      y <- 0
      z <- 0
    } else {
      x <- round(d$x,3)
      y <- round(d$y,3)
      z <- round(d$z,3)
    }
    
    withMathJax(paste(
      "\\begin{eqnarray}",
      "\\color{#806000}{\\alpha_1} & \\color{#806000}{=} & \\color{#806000}{", x,"} \\\\",
      "\\color{#002080}{\\alpha_2} & \\color{#002080}{=} & \\color{#002080}{",y,"} \\\\",
      "\\color{#408000}{\\alpha_3} & \\color{#408000}{=} & \\color{#408000}{",z,"} \\\\",
      "\\end{eqnarray}"
    ))
  })
  
}


# Ejecuta la aplicación
shinyApp(ui = ui, server = server)
