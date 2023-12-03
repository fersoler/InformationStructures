library(circlize)

#################################################################
# Functions to draw the 2D cones
# Author: Fernando Soler-Toscano - fsoler@us.es
#################################################################

# Normalize a vector
normVec <- function(vec){
  vec/sqrt(sum(vec^2))
}

# Given gammas, get the limits of the max. bio. cone
setMaxCone2D <- function(gammas){
  m <- -gammas
  v1 <- normVec(m[,1])
  v2 <- normVec(m[,2])
  c(getAnglePoint(v1), getAnglePoint(v2))
}

# Get the dimension of the cone 11
dimMaxCone2D <- function(gammas){
  cone <- setMaxCone2D(gammas)
  ((cone[2]-cone[1])%%360)/360
}

# Get the angle of a given points (degrees)
getAnglePoint <- function(p){
  if(p[2]>0){  # between 0 and pi
    angle <- acos(p[1])
  } else {
    angle <- 2*pi - acos(p[1])
  }
  360*angle/(2*pi)
}

# Write the laben indicating the cone size
writeConeSize <- function(init,end){
  if(init%%360 != end%%360){
    ang <- 2*pi*(init+end)/720
    # The color and size of the label can be changed here:
    text(1.2*cos(ang),1.2*sin(ang), round((end-init)/360, digits=2),cex=1,col="brown")
  }
}

# Write a label of the IS point:
writeLabel <- function(angle,sep,label){
  ang <- 2*pi*angle/360
  # The color and size of the label can be changed here:
  text(sep*cos(ang),sep*sin(ang), label, cex=.7,col="black")
}

# Draw the circle
# gamma: Matrix of interactions
# allCones: indicates if all cones are drawn or just the one with 2 species
# deawLabels: indicates if IS and labels are drawn.
drawCircleCones <- function(gamma, allCones = TRUE, drawLabels = TRUE){
  cone <- setMaxCone2D(gamma)
  if(allCones){
    labelsDraw <- 1:6
  } else {
    labelsDraw <- c(6)
  }
  # Colors can be changed here
  # The first list contains the colors of the cones
  # The second list, "light" versions of the first three colors
  colors1 <- c("#80FF00E0","#0080FFE0","#FF8000E0","#A0A0A0E0")
  colors2 <- c("#80FF0080","#0080FF80","#FF800080")
  plot(c(-1.5, 1.5), c(-1.5, 1.5), type = "n", axes = FALSE, ann = FALSE, asp = 1)
  coneLabelsSize <- matrix(0, ncol=2, nrow=6)
  # cone 00 (always 3rd quadrant)
  if(allCones) draw.sector(180, 270,clock.wise = FALSE, col = colors1[4], border = NA)
  coneLabelsSize[1,] <- c(180, 270)
  iss <- data.frame(180,270,"00",stringsAsFactors = FALSE)
  if(cone[1]>180){
    # Cone 10
    if(allCones) draw.sector(270, cone[1],clock.wise = FALSE, col = colors1[3], border = NA)
    coneLabelsSize[2,] <- c(270,cone[1])
    iss[nrow(iss)+1,] <- c(270,cone[1], "0010")
    draw.sector(cone[1], 360, clock.wise = FALSE, col = colors2[1], border = NA)
    labelsDraw <- c(labelsDraw, 3)
    coneLabelsSize[3,] <- c(cone[1],360)
    iss[nrow(iss)+1,] <- c(cone[1],360,"001011")
    startMainCone <- 0
  }
  if(cone[1]<=180){
    if(allCones) draw.sector(270, 360,clock.wise = FALSE, col = colors1[3], border = NA)
    coneLabelsSize[2,] <- c(270,360)
    iss[nrow(iss)+1,] <- c(270,360,"0010")
    if(allCones) draw.sector(360, cone[1], clock.wise = FALSE, col = colors2[3], border = NA)
    coneLabelsSize[3,] <- c(0,cone[1])
    iss[nrow(iss)+1,] <- c(0,cone[1], "000110")
    startMainCone <- cone[1]
  }
  if(cone[2]>90){
    if(allCones) draw.sector(cone[2], 180, clock.wise = FALSE, col = colors1[2], border = NA)
    coneLabelsSize[4,] <- c(cone[2], 180)
    iss[nrow(iss)+1,] <- c(cone[2], 180, "0001")
    draw.sector(90, cone[2], clock.wise = FALSE, col = colors2[1], border = NA)
    labelsDraw <- c(labelsDraw, 5)
    coneLabelsSize[5,] <- c(90, cone[2])
    iss[nrow(iss)+1,] <- c(90, cone[2], "000111")
    endMainCone <- 90
  }
  if(cone[2]<=90){
    if(allCones) draw.sector(90, 180, clock.wise = FALSE, col = colors1[2], border = NA)
    coneLabelsSize[4,] <- c(90, 180)
    iss[nrow(iss)+1,] <- c(90, 180, "0001")
    if(allCones) draw.sector(cone[2], 90, clock.wise = FALSE, col = colors2[2], border = NA)
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
  if(drawLabels){
    # Axes labels
    text(-1.3,0, expression('b'[1]),cex=1.3)
    text(0,-1.3, expression('b'[2]),cex=1.3)
    for(i in labelsDraw){
      writeConeSize(coneLabelsSize[i,1], coneLabelsSize[i,2])
      drawIScone(as.double(iss[i,1]), as.double(iss[i,2]), iss[i,3])
    }
  }
}

# Draw the IS graphs:
drawIScone <- function(begin, end, type){
  mid <- (end+begin)/2
  if(type == "00"){
    writeLabel(mid, .3, "00")
  }
  if(type == "0010"){
    writeLabel(mid, .3, "00")
    writeLabel(mid, .9, "10")
    drawArrow(mid,.3, mid, .9, "long")
  }
  if(type == "0001"){
    writeLabel(mid, .3, "00")
    writeLabel(mid, .9, "01")
    drawArrow(mid,.3, mid, .9, "long")
  }
  if(type == "000110" && begin != end){
    writeLabel(mid, .3, "00")
    writeLabel(mid, .6, "01")
    writeLabel(mid, .9, "10")
    drawArrow(mid,.3, mid, .6)
    drawArrow(mid,.6, mid, .9)
  }
  if(type == "001001" && begin != end){
    writeLabel(mid, .3, "00")
    writeLabel(mid, .6, "10")
    writeLabel(mid, .9, "01")
    drawArrow(mid,.3, mid, .6)
    drawArrow(mid,.6, mid, .9)
  }
  if(type == "000111" && begin != end){
    writeLabel(mid, .3, "00")
    writeLabel(mid, .6, "01")
    writeLabel(mid, .9, "11")
    drawArrow(mid,.3, mid, .6)
    drawArrow(mid,.6, mid, .9)
  }
  if(type == "001011" && begin != end){
    writeLabel(mid, .3, "00")
    writeLabel(mid, .6, "10")
    writeLabel(mid, .9, "11")
    drawArrow(mid,.3, mid, .6)
    drawArrow(mid,.6, mid, .9)
  }
  if(type == "00011011" && begin != end){
    writeLabel(mid, .3, "00")
    writeLabel((mid+3*end)/4, .6, "01")
    writeLabel((mid+3*begin)/4, .6, "10")
    writeLabel(mid, .9, "11")
    drawArrow(mid, .3,(mid+3*end)/4, .6)
    drawArrow(mid, .3,(mid+3*begin)/4, .6)
    drawArrow((mid+3*end)/4, .6, mid, .9)
    drawArrow((mid+3*begin)/4, .6, mid, .9)
  }
}

# Draw an arrow:
drawArrow <- function(angleSt, sepSt, angleEnd, sepEnd,type="short"){
  sep = 0.3
  if(type=="long"){
    sep = 0.2
  }
  aS <- 2*pi*angleSt/360
  aE <- 2*pi*angleEnd/360
  xS <- cos(aS)*sepSt
  yS <- sin(aS)*sepSt
  xE <- cos(aE)*sepEnd
  yE <- sin(aE)*sepEnd
  x0 <- xS+sep*(xE-xS)
  y0 <- yS+sep*(yE-yS)
  x1 <- xE-sep*(xE-xS)
  y1 <- yE-sep*(yE-yS)
  arrows(x0,y0,x1,y1, angle = 10, length=0.05,lwd = 1.5, code = 2)
}
