###############################################################################
# Auxiliary functions of the file subCones3d for converting subcones into     #
# stl files for 3d printing.                                                  #
###############################################################################


# Given a table of points 'initPoints' which are approx. along a geodesic, this 
# function finds the extreme points and the normal vector of the geodesic 
# approximating the points. 
adjustPoints2 <- function(all_points){
  # Take any points
  firstPoint <- all_points[1,]
  # Distances (angles) to all other points
  distances1 <- apply(all_points, 1, function(x) angle_between_vectors(x, firstPoint))
  # Greatest distante
  maxd <- max(abs(distances1))
  pos1 <- which(abs(distances1) == maxd[1])
  # Set first extreme point
  point1 <- all_points[pos1[1],]
  # All distances to point1
  distances2 <- apply(all_points, 1, function(x) angle_between_vectors(x, point1))
  maxd <- max(abs(distances2))
  pos2 <- which(abs(distances2) == maxd[1])
  point2 <- all_points[pos2[1],]
  
  # Normal vector (cross product)
  normal1 <-  normalize(pracma::cross(point1, point2))
  # Normal vector (using PCA)
  normal2 <- get_normal_vector(all_points)
  
  if(max(normal1 - normal2) > 0.000001){
    # Calculate angle between vectors
    angle <- acos(sum(normal1 * normal2))
    # Determine which is used
    if (abs(angle) < 0.01) {
      normal1 <- normal2
    }
  }
  list(points = rbind(point1, point2), normal = normal1)
}


# Calculate normal vector using PCA
get_normal_vector <- function(points_matrix) {
  # Center points
  centered_points <- scale(points_matrix, center = TRUE, scale = FALSE)
  
  # Apply PCA to get values and vectors
  pca <- svd(centered_points)
  # Normal vector
  normal_vector <- pca$v[, 3]
  
  # Ensure normal vector has correct orientation
  if (normal_vector[3] < 0) {
    normal_vector <- -normal_vector
  }
  normalize(normal_vector)
}

# Normalize a vector
normalize <- function(v) {
  v / sqrt(sum(v^2))
}

# Function to sort the edges of a subcone
sort_edges <- function(tableEdges) {
  
  # Create a copy of the table to select the edges
  rest_table <- tableEdges
  
  # Start with the first row as the initial edge
  current_edge <- rest_table[1, , drop = FALSE]
  sorted_edges <- list(current_edge)
  rest_table <- rest_table[-1, , drop = FALSE]
  
  # Set current final point (p2)
  p2_current <- normalize(current_edge[6:8])
  
  while (length(rest_table) > 0) {
    
    if (nrow(rest_table) == 0) break
    
    # Calculate distances of p2_current to limits p1 and p2 of each other edge
    all_distances <- apply(rest_table, 1, function(row) {
      p1 <- as.numeric(row[3:5])
      p2 <- as.numeric(row[6:8])
      min(distance(p2_current, p1), distance(p2_current, p2))
    })
    
    # Find the index of the edge with the closer endpoint
    idx <- which.min(all_distances)
    next_edge <- rest_table[idx, , drop = FALSE]
    
    # Reorder limits if necessary
    p1_next <- normalize(next_edge[3:5])
    p2_next <- normalize(next_edge[6:8])
    
    d1 <- distance(p2_current, p2_next)
    d2 <- distance(p2_current, p1_next)
    
    if (d1 > d2) {
      p2_current <- p2_next
    } else {
      # Change p1 and p2 if p1 is closer to current extreme
      next_edge[3:5] <- p2_next
      next_edge[6:8] <- p1_next
      p2_current <- p1_next
    }
    
    # Add the edge to the sorted list and remove from the table
    sorted_edges <- append(sorted_edges, list(next_edge))
    rest_table <- rest_table[-idx, , drop = FALSE]
  }
  
  # Combine sorted edges in a table
  do.call(rbind, sorted_edges)
}

# Euclidean distance
distance <- function(p1, p2) {
  sqrt(sum((p1 - p2)^2))
}

# This function adjusts the approximated limits of each subcone to the instersection
# of the geodesic defining its limits
adjustIntersections <- function(tableLimits0){
  tableLimits <- tableLimits0
  for(i in 1:nrow(tableLimits0)){
    row1 <- i
    row2 <- (i%%nrow(tableLimits0))+1
    # Normal of the geodesic corresponding the one side
    normal1 <- normalize(tableLimits0[row1,9:11])
    # Normal of the geodesic of the following side
    normal2 <- normalize(tableLimits0[row2,9:11])
    # Intersection of the geodesics
    crossP <- normalize(pracma::cross(normal1, normal2))
    # Change to the opposite point if closer
    crossAlt <- -crossP
    if(distance(normalize(tableLimits0[row1,6:8]), crossAlt) <
       distance(normalize(tableLimits0[row1,6:8]), crossP)){
      crossP <- crossAlt
    }
    # In case it doesn't work, take the mean of the previous limits
    if(distance(normalize(tableLimits0[row1,6:8]), crossP) > 0.05){
      crossP <- (tableLimits0[row1,6:8]+tableLimits0[row2,3:5])/2
    }
    # Set the extremes to the intersection
    tableLimits[row1,6:8] <- crossP
    tableLimits[row2,3:5] <- crossP
  }
  tableLimits  
}

# This function removes non-significant zones with less than minT triangles
# It receives the list of triangles, their colors and the minimum number of 
# triangles with any color to be kept. It removes from the list the
# triangles of a color with less than 'minT' triangles. 
cleanNonSignificantZones <- function(triangles, colorCodes, minT){
  
  colorCodesNew <- colorCodes
  tableColors <- table(colorCodes)
  
  for(code in rownames(tableColors)){
    nT <- tableColors[code]
    if(nT<minT){
      colorCodesNew[which(colorCodes == code)] <- "0"
    }
  }
  trianglesNew <- triangles[which(colorCodesNew != "0"),]
  colorCodesNew <- colorCodesNew[which(colorCodesNew != "0")]
  
  list(triangles = trianglesNew, colorCodes = colorCodesNew)
  
}

# This function adjusts the limits of the subcones given relevant points which
# are known to be intersections of subcones. Points closed than 'thr' to a 
# known limit are set to thar limit. Also, points closed then 'thr' among them
# are changed to the mean of all them. 
adjustWithRelevantPoints <- function(limitsTable, relevantPoints, thr){
  
  # Adjust to relevant points
  limitsNew <- limitsTable
  for(nl in 1:nrow(limitsNew)){
    
    distancesInit <- apply(relevantPoints, 1, 
                           function(x) distance(x,limitsNew[nl,3:5]))
    minDInit <- min(distancesInit)
    if(minDInit < thr){
      limitsNew[nl, 3:5] <- relevantPoints[which(distancesInit == minDInit),]
    }
    
    
    distancesEnd <- apply(relevantPoints, 1, 
                          function(x) distance(x,limitsNew[nl,6:8]))
    minDEnd <- min(distancesEnd)
    if(minDEnd < thr){
      limitsNew[nl, 6:8] <- relevantPoints[which(distancesEnd == minDEnd),]
    }
  }
  
  # Equate close points
  for(nl in 1:nrow(limitsNew)){
    
    point1 <- limitsNew[nl, 3:5]
    point2 <- limitsNew[nl, 6:8]
    close1Left <- which(apply(limitsNew[,3:5], 1,
                              function(x) distance(x, point1)) < thr)
    close1Right <- which(apply(limitsNew[,6:8], 1,
                               function(x) distance(x, point1)) < thr)
    close2Left <- which(apply(limitsNew[,3:5], 1,
                              function(x) distance(x, point2)) < thr)
    close2Right <- which(apply(limitsNew[,6:8], 1,
                               function(x) distance(x, point2)) < thr)
    
    pointsClose1 <- rbind(limitsNew[close1Left, 3:5],
                          limitsNew[close1Right, 6:8])
    newp1 <- colSums(pointsClose1)/nrow(pointsClose1)
    
    pointsClose2 <- rbind(limitsNew[close2Left, 3:5],
                          limitsNew[close2Right, 6:8])
    newp2 <- colSums(pointsClose2)/nrow(pointsClose2)
    
    limitsNew[close1Left, 3:5] <- matrix(newp1, nrow = length(close1Left), ncol = 3, byrow = TRUE)
    limitsNew[close1Right, 6:8] <- matrix(newp1, nrow = length(close1Right), ncol = 3, byrow = TRUE)
    limitsNew[close2Left, 3:5] <- matrix(newp2, nrow = length(close2Left), ncol = 3, byrow = TRUE)
    limitsNew[close2Right, 6:8] <- matrix(newp2, nrow = length(close2Right), ncol = 3, byrow = TRUE)
    
  }
  
  limitsNew
  
}

# Function to calculate the angle between two vectors
angle_between_vectors <- function(u, v) {
  # dot product
  dot_product <- sum(u * v)
  # Magnitudes
  norm_u <- sqrt(sum(u^2))
  norm_v <- sqrt(sum(v^2))
  # Cosine of the angle
  cos_theta <- dot_product / (norm_u * norm_v)
  # ensure it's in the range [-1,1]
  cos_theta <- max(min(cos_theta, 1), -1)
  # Angle
  acos(cos_theta)
}

# Create a 3d mesh of a subcone starting at its limit points. 
# - pointsCone. Points representing the (sorted) limits of the subcone. They are
# on the sphere of radius 1. 
# - TopRad. Radius of the sphere with the top surface of the cone. 
# - BotRad. Radius of the sphere with the bottom surface of the cone. 
# - resol. Division of the cone edges. 
# - top (bool). Indicates if the top surface is included in the mest. 
meshSubConeFromPoints <- function(pointsCone, TopRad, BotRad, resol, top = TRUE){
  
  rawMesh <- BuildMeshPattern(resol)
  # Scheme of the points
  pointSch <- rawMesh$vb
  # Number of points of each triangle
  oneTriP <- nrow(pointSch)
  
  # Sides of the scheme
  side1 <- which(pointSch[,1] == 0)
  side2 <- which(pointSch[,2] == 0)
  side3 <- which(pointSch[,3] == 0)
  
  # Scheme of side triangles
  nPSide <- length(side1)
  triangScheme <- NULL
  for(i in 1:(nPSide-1)){
    triangScheme <- rbind(triangScheme, 
                          c(i, i+1, i+nPSide))
    triangScheme <- rbind(triangScheme, 
                          c(i+nPSide, i+1+nPSide, i+1))
  }
  
  nPoints <- 0
  points <- NULL
  triangles <- NULL
  
  for(v in 2:(nrow(pointsCone)-1)){
    v1 <- pointsCone[1,]
    v2 <- pointsCone[v,]
    v3 <- pointsCone[v+1,]
    
    gammas <- cbind(v1,v2,v3)
    
    newPoints <- setPointsCone(gammas,c(1,1,1),rawMesh$vb)
    
    triangles  <- rbind(triangles, 
                        rawMesh$it+(nPoints+oneTriP))
    
    if(top){
      triangles <- rbind(triangles, rawMesh$it+nPoints)
    }
    
    # In the first case, add side 3
    if(v == 2){
      sidePoints <- append(side3, side3+oneTriP)
      triangles <- rbind(triangles,
                         nPoints+matrix(sidePoints[triangScheme], ncol = 3))
    }
    
    # Add schemes to side 1
    sidePoints <- append(side1, side1+oneTriP)
    triangles <- rbind(triangles, 
                       nPoints+matrix(sidePoints[triangScheme], ncol = 3))
    
    # In the last case, add side 2
    if(v == (nrow(pointsCone)-1)){
      sidePoints <- append(side2, side2+oneTriP)
      triangles <- rbind(triangles,
                         nPoints+matrix(sidePoints[triangScheme], ncol = 3))
    }
    
    points <- rbind(points, t(newPoints*TopRad), 
                    t(newPoints*BotRad))
    nPoints <- nrow(points)
  }
  list(points = points, triangles = triangles)
}
