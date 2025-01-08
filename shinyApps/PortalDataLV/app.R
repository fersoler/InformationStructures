library(shiny)
library(tidyverse)
library(plotly)
source("R/ISgraph.R")
source("R/figs/sphereCones.R")
source("R/figs/sphereSubCones.R")
source("R/setParameters.R")

global_PI <- reactiveVal(400)

#### COMMON CODE

#### LOAD DATA
PortalR_abundances <- read.csv("data/PortalData_abundances.csv")
# Data table
data <- PortalR_abundances[,c(1,3:5)]
# Abbreviations of the three species
spNs <- colnames(data)[2:4]
# Avoid 0 values (because of the use of log)
data[data == 0] <- 0.1
# Divide by 50 to scale resulting parameters 
data[,2:4] <- (data[,2:4]/50)

#### SETTING PV PARAMS
# Obtain alpha and gamma parameters: 
winW <- 10 # Margins to set alpha parameters
params <- getLVparamsNonAuto(as.matrix(data[,2:4]),winW)
# Matrix of interespecific interactions: 
gammas <- params$gammas
# Intrinsic growt rates
alphas <- data.frame(times = data[1:448,1], A1 = params$alphas[,1], A2 = params$alphas[,2], A3 = params$alphas[,3])
normAlphas <- as.data.frame(t(apply(alphas[, 2:4], 1, normalize)))
#### RECONSTRUCT TIME SERIES
reconstructPoints <- NULL
for(pi2 in 1:nrow(alphas)){
  theS2 <- getEvolLV3(as.matrix(alphas[pi2,2:4]), gammas, data[pi2,2:4], 1000, data[pi2,1], data[pi2+1,1])
  reconstructPoints <- rbind(reconstructPoints, theS2[1000,])
}
dataF <- data.frame(x = data[,1], SP1 = data[,2], SP2 = data[,3], SP3 = data[,4])

# For 3d figures
scene = list(xaxis = list(title = spNs[1], scaleratio = 1),
             yaxis = list(title = spNs[2], scaleratio = 1),
             zaxis = list(title = spNs[3], scaleratio = 1))


#### RECONSTRUCTED FIGURE
figT1 <- plot_ly(dataF, x = ~x, y = ~SP1, name = spNs[1], 
               type = 'scatter', mode = 'line', color = 'A')
figT1 <- figT1 %>% add_trace(y = ~SP2, name = spNs[2], color = 'B',
                         type = 'scatter', mode = 'line')
figT1 <- figT1 %>% add_trace(y = ~SP3, name = spNs[3], color = 'C', 
                         type = 'scatter', mode = 'line')
figT1 <- figT1 %>% add_trace(x = reconstructPoints[,1], y = reconstructPoints[,2], 
                         name = paste('LV',spNs[1]), color = 'A',
                         type = 'scatter', mode = 'markers')
figT1 <- figT1 %>% add_trace(x = reconstructPoints[,1], y = reconstructPoints[,3], 
                         name = paste('LV',spNs[2]), color='B',
                         type = 'scatter', mode = 'markers')
figT1 <- figT1 %>% add_trace(x = reconstructPoints[,1], y = reconstructPoints[,4], 
                         name = paste('LV',spNs[3]), color = 'C',
                         type = 'scatter', mode = 'markers')
figT1 <- figT1 %>% layout(
  title = "Real and reconstructed (LV) data",
  xaxis = list(title = "Period (months)"),
  yaxis = list (title = "Abundance"))

#### FIGURE WITH ALPHAS
figT2 <- plot_ly(alphas, x = ~times, y = ~A1, name = spNs[1], 
               type = 'scatter', mode = 'lines+markers', color='A')
figT2 <- figT2 %>% add_trace(x = ~times, y = ~A2, name = spNs[2], 
                         type = 'scatter', mode = 'lines+markers', color='B')
figT2 <- figT2 %>% add_trace(x = ~times, y = ~A3, name = spNs[3], 
                         type = 'scatter', mode = 'lines+markers', color='C')
figT2 <- figT2 %>% layout(
  title = "LV intrinsic growth rates",
  xaxis = list(title = "Period (months)"),
  yaxis = list (title = "Intrinsic growth rate"))

#### FOR THE SPHERE
## Sphere with max bio cone (green) and cones with 1 species (blue)
# Intrinsic growth rates are represented by a point

# We set n to divide each cone into n^2 triangles
# The higher value, the better visualization
n <- 40
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


##### FOR THE IS


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
  IS3 <- ISbuild(aVals,gammas)
  gr3 <- ISgraph(IS3, 1:3)
  ISgraphDrawPie2(IS3,gr3, ISgraphLayout(IS3, gr3, "3Dcube"), c("#66c2a5", "#fc8d61", "#8da0cb"))
}

#############################################################################
# Shiny app
#############################################################################

# the app UI
ui <- fluidPage(
  
  titlePanel("Setting gLV parameters"),
  
  tags$p("Click a point on the top plots", style = "text-align:left;"),
  
  hr(),
  
  fluidRow(column(width=6, plotlyOutput("plotTop1")),
           column(width=6, plotlyOutput("plotTop2"))),
          
  fluidRow(
           column(width=4,plotlyOutput("plotSphere")), 
           column(width=4,plotOutput("plotIS")),
           column(width=4,plotlyOutput("plotEvo")))
 
)

# app logic
server <- function(input, output, session) {

  observeEvent(event_data("plotly_click"), {
    event_data <- event_data("plotly_click")
    req(!is.null(event_data))
    if(round(event_data$x) == event_data$x){
      pi <- which(data[,1] == event_data$x)
      global_PI(pi)
     }
  })
  
  # Top1 plot
  output$plotTop1 <- renderPlotly({figT1})
  
  # Top2 plot
  output$plotTop2 <- renderPlotly({figT2})

  # Sphere
  output$plotSphere <- renderPlotly({
    
    req(!is.null(global_PI))
    pi <- global_PI()
    
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
      add_markers(x = 1.03*normalize(alphas[pi,2:4])[1,1],
                  y = 1.03*normalize(alphas[pi,2:4])[1,2], 
                  z = 1.03*normalize(alphas[pi,2:4])[1,3],
                  marker = list(size = 4, color = 'yellow')) %>% 
      layout(showlegend=FALSE, scene = scene)
    
  })
  
  # IS
  output$plotIS <- renderPlot({
    
    req(!is.null(global_PI))
    pi <- global_PI()
    
    datos <- showIS(as.matrix(alphas[pi,2:4]))
    
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
  
  # Evolution
  output$plotEvo <- renderPlotly({
    
    req(!is.null(global_PI))
    pi <- global_PI()
    
    IS3 <- ISbuild(as.matrix(alphas[pi,2:4]), gammas)
    # Path to the next period: 
    path <- getEvolLV3(as.matrix(alphas[pi,2:4]), gammas, data[pi,2:4], 100, 1, 2)
    # Path in 1000 more periods: 
    path2 <- getEvolLV3(as.matrix(alphas[pi,2:4]), gammas, data[pi,2:4], 10000, 1, 1000)
    
    # Figure
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
  })
  
  
}

# Run the app:
shinyApp(ui = ui, server = server)

