library(plotly)
source("R/setParameters.R")

# Set matrix A
o_A <- 1
(A <- o_A * rbind(c(-1, 0.25),
                 c(0.5,  -1)))

# Set intrinsic growth rates
t_values <- 1:100   # Time
o <- 0.8
b1 <- o + 1 * cos(0.05*t_values + 1)
b2 <- o + 1.3 * sin(0.05*t_values)

# Dataframe with bs
bData <- data.frame(time = t_values, b1 = b1, b2 = b2)


# Look at the intrinsic growth rates: 
fig <- plot_ly(bData, x = ~time, y = ~b1, name = 'b1', 
               type = 'scatter', mode = 'line')
fig <- fig %>% add_trace(x = ~time, y = ~b2, name = 'b2', 
                         type = 'scatter', mode = 'line')
fig <- fig %>% layout(
  title = "LV intrinsic growth rates",
  xaxis = list(title = "Time"),
  yaxis = list (title = "Intrinsic growth rate (b1, b2)"))
fig


# Data (simulated time series with given A and b)
data <- data.frame(time=c(1), 
                   u1 = c(1.5),  # Starting value for sp1 
                   u2 = c(2))  # Starting value for sp2
for(t in t_values){
  #l <- 1000
  sol <- getEvolLV2(as.matrix(bData[t,2:3]), A, data[t,2:3], 2, t, t+1)
  data <- rbind(data,sol[2,])
}

# Look at the time series
fig <- plot_ly(data, x = ~time, y = ~u1, name = 'u1', 
               type = 'scatter', mode = 'line')
fig <- fig %>% add_trace(x = ~time, y = ~u2, name = 'u2', 
                         type = 'scatter', mode = 'line')
fig <- fig %>% layout(
  title = "Time series",
  xaxis = list(title = "Time"),
  yaxis = list (title = "Value"))
fig

########################################################################
# Estimate parameters
########################################################################


# Get the estimated matrix (Rafael's idea):
timeSelect <- 1:101 # 1:8  # 1:30  # 31:1000
windW <- 10
params <- getLVparamsNonAuto_weighted_A(as.matrix(data[timeSelect,2:3]), 10)
(A2 <- params$gammas)
# The real A:  
A
# # LV_map
# timeSelect <- 1:101  # 1:30  # 31:1000
# theta <- 0.2
# #LV_map <- LV_map_A_const(as.matrix(data[timeSelect,2:3]),theta)
# LV_map2 <- LV_map(as.matrix(data[timeSelect,2:3]),theta)
# 
# (A2 <- LV_map2$A[[10]])
# # The real A:  
# A
# 
# #bCalculated <- t(LV_map$B)
# bCalculated <- LV_map2$B
bCalculated <- params$alphas

# Get b values
#windW <- 10  # Window size to left and right
#bCalculated <- get_b_given_A(as.matrix(data[,2:3]), A2, windW)
#bCalculated

# Look at the obtained intrinsic growth rates: 
fig <- plot_ly(bData, x = ~time, y = ~b1, name = 'b1', 
               type = 'scatter', mode = 'line', color='A')
fig <- fig %>% add_trace(x = ~time, y = ~b2, name = 'b2', 
                         type = 'scatter', mode = 'line', color='B')
fig <- fig %>% add_trace(x = ~time, y = bCalculated[,1], name = 'b1 calculated', 
                         type = 'scatter', mode = 'markers', color='A', marker = list(size = 4))
fig <- fig %>% add_trace(x = ~time, y = bCalculated[,2], name = 'b2 calculated', 
                         type = 'scatter', mode = 'markers',color='B',  marker = list(size = 4))
fig <- fig %>% layout(
  title = "LV intrinsic growth rates",
  xaxis = list(title = "Time"),
  yaxis = list (title = "Intrinsic growth rates (b1, b2)"))
fig

# Reconstruct the timeseries
data2 <- data.frame(time=c(1), 
                    u1 = c(4),  # Starting value for sp1 
                    u2 = c(4))  # Starting value for sp2
for(t in t_values){
  sol <- getEvolLV2(as.matrix(bCalculated[t,]), A, data[t,2:3], 2, t, t+1)
  data2 <- rbind(data2,sol[2,])
}

fig <- plot_ly(data, x = ~time, y = ~u1, name = 'u1', 
               type = 'scatter', mode = 'line', color='A')
fig <- fig %>% add_trace(x = ~time, y = ~u2, name = 'u2', 
                         type = 'scatter', mode = 'line',color='B')
fig <- fig %>% add_trace(x = ~time, y = data2[,2], name = 'u1 calculated', 
                         type = 'scatter', mode = 'markers', color='A', marker = list(size = 4))
fig <- fig %>% add_trace(x = ~time, y = data2[,3], name = 'u2 calculated', 
                         type = 'scatter', mode = 'markers',color='B',  marker = list(size = 4))
fig <- fig %>% layout(
  title = "Time series",
  xaxis = list(title = "Time"),
  yaxis = list (title = "Value"))
fig

