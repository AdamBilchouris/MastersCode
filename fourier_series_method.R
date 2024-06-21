library(RandomFields)
library(plot3D)
library(doParallel)
library(foreach)
library(glue)
library(dplyr)

fourierSeries_direct <- function(x, y) {
  N <- length(x)
  ans <- c()
  bns <- c()
  for(i in 0:N) {
    ans <- c(ans, sum(y * cos(i*x)))
    bns <- c(bns, sum(y * sin(i*x)))
  }
  
  ans <- ans / N
  bns <- bns / N
  
  fSeries <- c()
  for(i in 1:length(x)) {
    temp <- c()
    for(j in 2:length(ans)) {
      cosS <- ans[j]*cos(x[i]*(j-1))
      sinS <- bns[j]*sin(x[i]*(j-1))
      temp <- c(temp, cosS + sinS)
    }
    fSeries <- c(fSeries, (ans[1]/2) + sum(temp))
  }
  
  return(list('x'=x, 'y'=y, 'fittedY'=fSeries, 'a0'=ans[1], 'ans'=ans[-1], 'bns'=bns[-1]))
}

midpointCircle <- function(xC, yC, r) {
  theta <- seq(0, 2*pi, length.out=6*r)
  circlePointsX <- r*cos(theta) + xC
  circlePointsY <- r*sin(theta) + yC
  
  points <- data.frame(x=circlePointsX, y=circlePointsY)
  return(points)
}


getCircleSample <- function(xC, yC, r, z, xRes) {
  points <- midpointCircle(xC, yC, r)
  # Restrict the points to within the size of matrix.
  points <- points %>% dplyr::filter(x >= 1, y >= 1, x <= nrow(z), y <= ncol(z))
  pointsSample <- c()
  for(i in 1:nrow(points)) {
    pointsSample <- c(pointsSample, z[points[i, 'x'], points[i, 'y']])
  }

  # upper half circle
  isOdd <- T
  if(length(pointsSample) %% 2 == 0) {
    lps <- length(pointsSample) / 2
    isOdd <- F
  } else {
    lps <- (1+length(pointsSample))/2
  }
  x_0_pi <- seq(0, pi, length.out=lps) 
  fs <- fourierSeries_direct(x_0_pi, pointsSample[1:lps])
  fs$x <- 2*(r*xRes)*sin(x_0_pi / 2)
  
  # lower half circle
  x_pi_2pi <- seq(pi, 2*pi, length.out=lps)

  # check if odd
  if(isOdd) {
    fs2 <- fourierSeries_direct(x_pi_2pi, pointsSample[lps:length(pointsSample)])
  } else {
    fs2 <- fourierSeries_direct(x_pi_2pi, pointsSample[(lps + 1):length(pointsSample)])
  }
  fs2$x <- 2*(r*xRes)*sin(x_pi_2pi / 2)
  
  return(list('upper'=fs, 'lower'=fs2))
}

registerDoParallel(detectCores())

# Generate a random field, or read one from memory
# x <- seq(-10, 10, 0.1)
# xLen <- length(x)
# y <- x

# xRes <- abs(x[1] - x[2])
# yRes <- abs(y[1] - y[2])
# cbind(xRes, yRes)

# model <- RMgauss(var=1, scale=1)
# sim <- RFsimulate(model, x=x, y=y)
# zMatrix <- RFspDataFrame2conventional(sim)$data
# image2D(x=x, y=y, z=zMatrix)

# Radius
# R <- 10
# If we want a radius of 1, we need
# R <- 1/xRes
R <- 30/xRes
#image2D(x=x, y=y, z=zMatrix)
maxCircleX <- round( length(x) / (2*R) )
maxCircleY <- round( length(y) / (2*R) )

# [-10, 10, 0.1]x[-10, 10, 0.1] grid
# Start at the top left corner (1, 201), end in bottom right corner (201, 1).
initialCentreX <- 1+R
initialCentreY <- 1+R
lastCentreX <- length(x) - R
lastCentreY <- length(x) - R

k <- 10
stepX <- floor((lastCentreX - initialCentreX)/k)
stepY <- floor((lastCentreY - initialCentreY)/k)

#image2D(x=1:length(x), y=1:length(y), z=zMatrix)
#lines(midpointCircle(initialCentreX + 2*1*R, initialCentreY - 2*1*R, R), col='red', lwd=2)

# image2D(x=1:length(x), y=1:length(y), z=zMatrix)
fsList <- list()
sampleNo <- 1
# Iterate over rows first (y is 0 initally)
startT <- Sys.time()
for(j in 0:k) {
  for(i in 0:k){
    # lines(midpointCircle(initialCentreX +i*stepX, initialCentreY+j*stepY, R), col='red', lwd=2)
    fs <- getCircleSample(initialCentreX +i*stepX, initialCentreY+j*stepY, R, zMatrix, xRes)
    fsList[[sampleNo]] <- fs[['upper']]
    sampleNo <- sampleNo + 1
    fsList[[sampleNo]] <- fs[['lower']]
    sampleNo <- sampleNo + 1
  }
}

# Compute C0 and Cl
a0s <- c()
for(i in 1:length(fsList)) {
  a0s <- c(a0s, fsList[[i]]$a0)
}
c0 <- mean(a0s^2) / 4
c0

ls <- length(fsList[[1]]$ans)
cls <- list()
cls2 <- list()
for(l in 1:ls) {
  ans <- c()
  bns <- c()
  for(i in 1:length(fsList)) {
    ans <- c(ans, fsList[[i]]$ans[l])
    bns <- c(bns, fsList[[i]]$bns[l])
    
    if(NA %in% ans) {
      print(c(i, l))
    }
  }
  cls[[l]] <- (1/4) * (mean(ans^2) + mean(bns^2))

}

thetas <- seq(0, pi, length.out = length(fsList[[1]]$x))

rho <- c()
cls2 <- cls[1:(length(cls))]
for(t in thetas) { 
  temp <- c()
  for(l in 1:length(cls2)) {
    temp <- c(temp, 2*cls2[[l]] * cos(l*t))
  }
  rho <- c(rho, c0 + sum(temp))
}

endT <- Sys.time()
print(endT - startT)

seqRF <- fsList[[1]]$x
plot(seqRF, rho/max(rho))
