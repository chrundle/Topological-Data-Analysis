# ------ Install and reference required packages ------ #
# Install TDA package for Topological Data Analysis (comment out after installed)
#install.packages('TDA')

# Install package for Delaunay and Voronoi package (comment out after installed)
#install.packages('deldir')

# Reference installed packages
library(TDA)
library(deldir)

# ------ Functions used in lab problems ------ #
# Function that samples points from annulus
sample.annulus <- function(num.pts, inner.radius, outer.radius){
  # Sample point from uniform distribution to generate angle
  theta <- runif(num.pts) * 2 * pi
  # Compute radius of point on annulus
  radius <- sqrt(runif(num.pts, inner.radius^2, outer.radius^2))
  # Create x and y coordinate of point on annulus
  x <- radius * cos(theta)
  y <- radius * sin(theta)
  # Combine x and y coordinates into single point
  cbind(x,y)
}

# ------ Generate and plot sample points for lab ------ #
# Call sample function to generate points
X <- sample.annulus(num.pts = 100, inner.radius = 1, outer.radius = 2)
# Generate plot of sampled points
plot(X, pch=20, col='blue', asp=1)

# ------ Delaunay complex (also called alpha complex) ------ #
DelVor <- deldir(X[,1], X[,2], suppressMsge = TRUE)
# Generate Voronoi cells for sampled points
plot(DelVor, pch=20, col=c('black', 'red', 'blue'), wlines=('tess'))
# Delaunay complex
plot(DelVor, pch=20, col=c('black', 'red', 'blue'), wlines=('triang'))
# Voronoi cells and their dual (the Delaunay complex)
plot(DelVor, pch=20, col=c('black', 'red', 'blue'), wlines=('both'))
# Compute Persistent Homology (PH)
PH.output <- alphaComplexDiag(X)
# Get Persistence Diagram (PD) from output of persistent homology computation
PD <- PH.output[["diagram"]]
# For some reason, birth and death values have been squared; take square root
PD[,2] <- sqrt(PD[,2])
PD[,3] <- sqrt(PD[,3])
# Plot the persistence diagram
plot(PD, asp=1, diagLim=c(0,1.3))
legend(0.6, 0.4, c('Homology in degree 0', 'Homology in degree 1'),
       col = c(1,2), pch=c(19,2), cex=0.8, pt.lwd=2)
# Plot the bar code
plot(PD, diagLim=c(0,1.3), barcode=TRUE)

# ------ Vietoris-Rips complex ------ #
# Define euclidean distance function
eucl.dist <- function(u,v) sqrt(sum((u-v)^2))
# Set max filtration
max.filtration <- 2.2
# Generate plot
plot(X, pch=20, col='blue', asp=1)
# Set number of points
num.pts <- dim(X)[1]
# Loop through all points
for(i in 1:num.pts)
  for(j in 1:num.pts)
    # Check if euclidean distance is less than max filtration
    if (eucl.dist(X[i,], X[j,]) < max.filtration)
      # Plot connected line segment of two points that are within max filtration distance
      lines(rbind(X[i,], X[j,]))
# Compute persistent homology
PH.output <- ripsDiag(X, maxdimension=1, maxscale=max.filtration)
# Get persistence diagram from output of persistent homology computation
PD <- PH.output[["diagram"]]
# Plot the persistence diagram
plot(PD, asp=1, diagLim=c(0, max.filtration))
legend(1, 0.6, c('Homology in degree 0', 'Homology in degree 1'),
       col=c(1,2), pch=c(19,2), cex=0.8, pt.lwd=2)
# Plot the bar code
plot(PD, diagLim=c(0, max.filtration), barcode=TRUE)

# ------ Plot representative cycles for Delaunay comples
PH.output <- alphaComplexDiag(X, maxdimension=1, library=c("GUDHI", "Dionysus"),
                              location=TRUE)
PD <- PH.output[["diagram"]]
ones <- which(PD[,1] == 1)
persistence <- PD[ones,3] - PD[ones,2]
cycles <- PH.output[["cycleLocation"]][ones[order(persistence)]]
for(i in 1:length(cycles)){
  plot(X, pch=20, col='blue', asp=1)
  for(j in 1:dim(cycles[[i]])[1])
    lines(cycles[[i]][j,,])
}

# ------ Plot representative cycles for Vietoris-Rips comples
PH.output <- ripsDiag(X, maxdimension=1, maxscale=max.filtration,
                      library=c("GUDHI", "Dionysus"), location=TRUE)
PD <- PH.output[["diagram"]]
ones <- which(PD[, 1] == 1)
persistence <- PD[ones,3] - PD[ones,2]
cycles <- PH.output[["cycleLocation"]][ones[order(persistence)]]
for(i in 1:length(cycles)){
  plot(X, pch=20, col='blue', asp=1)
  for(j in 1:dim(cycles[[i]])[1])
    lines(cycles[[i]][j,,])
}
