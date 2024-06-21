library(RandomFields)
library(plot3D)
library(doParallel)
library(glue)
library(dplyr)
library(sf) # to construct a polygon and check if it is inside the polygon.
library(pracma) # orthogonal project
library(akima) # Library for interploation
library(spdep)

registerDoParallel(detectCores())

rot2D <- function(xy_vec, angle) {
  rot <- Rotation(xy_vec, angle)
  
  return(rot)
}

getPointsInTube <- function(start_vec, end_vec, angle, epsilon, x, y, e.type='y', parallel.line=F, y.shift=NA, just.line=F) {
  if(!e.type %in% c('y', 'box')) {
    errStr <- paste0("e.type is not one of 'y', 'box' | given: ", e.type)
    print(errStr)
    return(F)
  }
  
  rotLine <- rot2D(matrix(c(start_vec, end_vec), ncol=2, byrow=T), angle)
  m <- (rotLine[2, 2] - rotLine[1, 2]) / (rotLine[2, 1] - rotLine[1, 1])
  
  bVertical <- F
  bHorizontal <- F
  bNormal <- F
  
  # vertical line
  if(abs(m) == Inf | abs(m) > 1/.Machine$double.eps) {
    rotLineSeq <- data.frame(x=rotLine[1, 1], y=y)
    bVertical <- T
  }
  
  # Horizontal line
  else if(abs(m) < .Machine$double.eps | m == 0) {
    rotLineSeq <- data.frame(x=x, y=rotLine[1, 2])
    bHorizontal <- T
  }
  
  else {
    yLine <- m*x
    rotLineSeq <- data.frame(x=x, y=yLine)
    bNormal <- T
  }
  
  if(parallel.line) {
    rotLineSeq$y <- rotLineSeq$y + y.shift
  }
  rotLineSeq <- rotLineSeq %>% filter(x >= min(x), x <= max(x), y >= min(y), y <= max(y))
  rotLine_p1 <- rotLineSeq
  rotLine_m1 <- rotLineSeq
  
  if(just.line == T) {
    return(list('baseLine'=rotLineSeq, 'line_m1'=rotLine_m1, 'line_p1'=rotLine_p1))
  }
  
  if(e.type == 'box') {
    # Construct box around every point?
    boxPoints <- rotLineSeq
    # TL: top left, TR: top right, BL: bottom left, BR: bottom right
    # TL: (-epsilon, epsilon), TR: (epsilon, epsilon), BL: (-epsilon, -epsilon), BR: (epsilon, -epsilon)
    # xme: X - epsilon, xpe: X + epsilon, yme, ype are the same but for the y coordinates
    boxPoints$xme <- plyr::round_any(boxPoints$x - epsilon, abs(x[1] - x[2]))
    boxPoints$xpe <- plyr::round_any(boxPoints$x + epsilon, abs(x[1] - x[2]))
    boxPoints$yme <- plyr::round_any(boxPoints$y - epsilon, abs(y[1] - y[2]))
    boxPoints$ype <- plyr::round_any(boxPoints$y + epsilon, abs(y[1] - y[2]))
    
    withinEpsilon <- data.frame(x=numeric(0), y=numeric(0))
    for(i in 1:nrow(boxPoints)) {
      possibleXY <- expand.grid(seq(boxPoints[i, 'xme'], boxPoints[i, 'xpe'], by=abs(x[1] - x[2])),
                                seq(boxPoints[i, 'yme'], boxPoints[i, 'ype'], by=abs(y[1] - y[2])))
      colnames(possibleXY) <- c('x', 'y')
      for(j in 1:nrow(possibleXY)) {
        withinEpsilon[nrow(withinEpsilon) + 1, ] <- possibleXY[j, ]
      }
    }
    
    return(list('withinTube'=withinEpsilon, 'baseLine'=rotLineSeq))
  }
  
  if(e.type == 'y') {
    if(bVertical) {
      rotLine_p1$x <- rotLine_p1$x + epsilon
      rotLine_m1$x <- rotLine_m1$x - epsilon
    }
    
    else {
      rotLine_p1$y <- rotLine_p1$y + epsilon
      rotLine_m1$y <- rotLine_m1$y - epsilon
    }
    
    minXx <- min(x)
    maxXx <- max(x)
    minYy <- min(y)
    maxYy <- max(y)
    
    rotLine_p1 <- rotLine_p1 %>% filter(x > minXx, x < maxXx, y > minYy, y < maxYy)
    rotLine_m1 <- rotLine_m1 %>% filter(x > minXx, x < maxXx, y > minYy, y < maxYy)
    
    minY <- ifelse(min(rotLine_m1$y) < min(rotLine_p1$y), min(rotLine_m1$y), min(rotLine_p1$y))
    maxY <- ifelse(max(rotLine_m1$y) > max(rotLine_p1$y), max(rotLine_m1$y), max(rotLine_p1$y))
    
    candidateY <- c()
    for(yy in y) {
      if(yy >= minY && yy <= maxY) {
        candidateY <- c(candidateY, yy)
      }
    }
    
    minX <- ifelse(min(rotLine_m1$x) < min(rotLine_p1$x), min(rotLine_m1$x), min(rotLine_p1$x))
    maxX <- ifelse(max(rotLine_m1$x) > max(rotLine_p1$x), max(rotLine_m1$x), max(rotLine_p1$x))
    candidateX <- c()
    for(xx in x) {
      if(xx >= minX && xx <= maxX) {
        candidateX <- c(candidateX, xx)
      }
    }
    
    if(length(candidateY) == 0) {
      candidateY <- rep(plyr::round_any(minY, abs(y[1] - y[2])), length(candidateX))
    }
    
    if(length(candidateX) == 0) {
      candidateX <- rep(plyr::round_any(minX, abs(x[1] - x[2])), length(candidateY))
    }
    
    possiblePoints <- expand.grid(candidateX, candidateY)
    colnames(possiblePoints) <- c('x', 'y')
    
    A <- c(rotLine_m1$x[1], rotLine_m1$y[1])
    B <- c(rotLine_m1$x[nrow(rotLine_m1)], rotLine_m1$y[nrow(rotLine_m1)])
    C <- c(rotLine_p1$x[1], rotLine_p1$y[1])
    D <- c(rotLine_p1$x[nrow(rotLine_p1)], rotLine_p1$y[nrow(rotLine_p1)])
    
    pol.x <- c(A[1], C[1], D[1], B[1], A[1])
    pol.y <- c(A[2], C[2], D[2], B[2], A[2])
    
    withinEpsilon <- data.frame(x=numeric(0), y=numeric(0))
    for(i in 1:nrow(possiblePoints)) {
      isIn <-  point.in.polygon(possiblePoints[i, 'x'], possiblePoints[i, 'y'], pol.x, pol.y)
      if(isIn) {
        withinEpsilon[nrow(withinEpsilon) + 1, ] <- c(possiblePoints[i, 'x'], possiblePoints[i, 'y'])
      }
    } 
    return(list('withinTube'=withinEpsilon, 'baseLine'=rotLineSeq, 'line_m1'=rotLine_m1, 'line_p1'=rotLine_p1, 'corners'=list('A'=A, 'B'=B, 'C'=C, 'D'=D)))
  }
}

convertToIndex <- function(points, x, y) {
  outDf <- data.frame(x=numeric(0), y=numeric(0))
  for(i in 1:nrow(points)) {
    # https://stackoverflow.com/a/39251375
    possibleX <- which(abs(x - points[i, 1]) - min(abs(x - points[i, 1])) < sqrt(.Machine$double.eps))
    possibleY <- which(abs(y - points[i, 2]) - min(abs(y - points[i, 2])) < sqrt(.Machine$double.eps))
    if(length(possibleX) == 0 | length(possibleY) == 0) {
      next
    }
    outDf[nrow(outDf) + 1, ] <- c(possibleX, possibleY)
  }
  return(outDf)
}

euclideanDist <- function(x, y) {
  return( sqrt(sum((x - y)^2)) )
}

weightedCardinal <- function(vec, x, y, z, xRes, yRes) {
  up <- c(plyr::round_any(vec[1], xRes), plyr::round_any(vec[2] + yRes, yRes))
  down <- c(plyr::round_any(vec[1], xRes), plyr::round_any(vec[2] - yRes, yRes))
  left <- c(plyr::round_any(vec[1] - xRes, xRes), plyr::round_any(vec[2], yRes))
  right <- c(plyr::round_any(vec[1] + xRes, xRes), plyr::round_any(vec[2], yRes))
  
  d_up <- euclideanDist(vec, up)
  d_down <- euclideanDist(vec, down)
  d_left <- euclideanDist(vec, left)
  d_right <- euclideanDist(vec, right)
  
  weights <- 1/c(d_up, d_down, d_left, d_right) / sum(d_up, d_down, d_left, d_right)
  weights <- weights / sum(weights)
  idxArr <- convertToIndex(matrix(c(up, down, left, right), ncol=2, byrow=T), x, y)
  
  vals <- c()
  for(i in 1:nrow(idxArr)) {
    vals <- c(vals, z[idxArr[i, 1], idxArr[i, 2]])
  }
  
  return(sum(weights * vals))
}

interpolate <- function(x_coord, y_coord, x, y, z, method='bilinear') {
  if(method == 'bilinear') {
    return(bilinear(x, y, zMatrix, x_coord, y_coord))
  }
  else if(method == 'bicubic') {
    return(bicubic(x, y, zMatrix, x_coord, y_coord))
  }
  else {
    return("Method is not  'bilinear' or 'bicubic'")
  }
}

computeTube <- function(x, y, z, xRes, yRes, epsilon, start_vec, end_vec, angle, method='bicubic', type='y', parallel.line=F, y.shift=NA) {
  # Compute the points within the epsilon tube
  pointsInTube <- getPointsInTube(start_vec, end_vec, angle, epsilon, x, y, type, parallel.line, y.shift)
  if(typeof(pointsInTube) == 'logical') {
    if(!pointsInTube) {
      return('ERROR: type unknown!')
    }
  }
  
  # Convert the points in the tube to indices that can be used with the z matrix
  indexDf <- convertToIndex(pointsInTube$withinTube, x, y)
  
  # Compute the projections
  # Project all points onto main line
  # Assume orthogonal projection since the line goes through the origin ?
  base_vec <- as.numeric(pointsInTube$baseLine[nrow(pointsInTube$baseLine), ])
  projected <- data.frame(t(linearproj(matrix(base_vec, 2, 1), t(as.matrix(pointsInTube$withinTube)))$Q))
  colnames(projected) <- c('x', 'y')
  
  # Get 'estimated' Z values
  # This method uses a weighted sum (closer => more weight) of four points in the cardinal (compass) directions.
  zVals <- c()
  if(method == 'cardinal') {
    for(i in 1:nrow(projected)) {
      zVals <- c(zVals, weightedCardinal(as.numeric(projected[i, ]), x, y, z, xRes, yRes))
    }
  }
  else if(method %in% c('bicubic', 'bilinear')) {
    zVals <- interpolate(projected[, 1], projected[, 2], x, y, z, method)$z
  }
  else {
    return(message("Error: Method is not one of 'cardinal', 'bilinear', 'bicubic'."))
  }
  
  # Do MDS
  dMat <- dist(projected)
  fit1D <- cmdscale(dMat, k=1)

  # Get min and max distance between two neighbouring times
  maxDist <- -1000
  minDist <- 1000
  for(i in 2:length(fit1D)) {
    d <- abs(fit1D[i-1] - fit1D[i])
    if(d > maxDist) {
      maxDist <- d
    }
    if(d < minDist) {
      minDist <- d
    }
  }
  
  binSize <- 0.1
  # Add a small number to ensure the min/max are in the bins
  bins <- seq(min(fit1D) - 0.1, max(fit1D) + 0.1, by=binSize)
  binGap <- abs(bins[1] - bins[2])
  
  # Construct new DF based off MDS and estimated Z vals
  mdsEst <- data.frame(x=fit1D, z=zVals)
  mdsEst <- mdsEst %>% mutate(bin=cut(x, bins))
  mdsEst <- mdsEst %>% mutate(bin_no=ntile(x, n=length(bins)))
  
  mdsEstGroup <- mdsEst[, !colnames(mdsEst) %in% c("x")]
  mdsEstGroup <- mdsEst %>% group_by(bin_no) %>% summarise_at("z", mean)
  mdsEstGroup$bin_no <- mdsEstGroup$bin_no*binGap

  
  estCorr <- acf(mdsEstGroup$z, type='covariance', lag.max=(nrow(mdsEstGroup)), plot=F)
  return(list('tubeDf'=pointsInTube, 'indexDf'=indexDf, 'baseVec'=base_vec, 'projected'=projected,
              'zVals'=zVals, 'dMat'=dMat, 'mds1D'=fit1D, 'minMaxDist'=c(minDist, maxDist), 'binSize'=binSize,
              'bins'=bins, 'binGap'=binGap, 'mdsEst'=mdsEst, 'mdsEstGroup'=mdsEstGroup,
              'estCorr'=estCorr, 'lags'=as.vector(estCorr$lag))) #'lags'=mdsEstGroup$bin_no))
}

# Generate a random field, or read one from memory
# x <- seq(-100, 100, 0.1)
# xLen <- length(x)
# y <- x
# 
# xRes <- abs(x[1] - x[2])
# yRes <- abs(y[1] - y[2])
# cbind(xRes, yRes)
# 
# model <- RMgauss(var=1, scale=1)
# model <- RMbessel(var=1, nu=1, scale=1)
# sim <- RFsimulate(model, x=x, y=y)
# zMatrix <- RFspDataFrame2conventional(sim)$data
# image2D(x=x, y=y, z=zMatrix)

# Loop over many angles.
thetas <- seq(0, pi, length.out=23)

startT <- Sys.time()
aa_par <- foreach(i=1:length(thetas), .combine='list', .multicombine=T, .packages=c("dplyr", "magrittr", "sf", "sp", "pracma", "akima", "spdep")) %dopar% {
  aa <- computeTube(x, y, zMatrix, xRes, yRes, 0.3, c(-10, 0), c(10, 0), thetas[i], 'cardinal', 'y')
  aa
}


# find the total set of lags in all lists
set_of_lags <- c()
for(i in 1:length(aa_par)) {
  set_of_lags <- c(set_of_lags, aa_par[[i]]$lags)
}

maxLag <- max(set_of_lags)

set_of_lags <- unique(set_of_lags)
set_of_lags

lagSep <- abs(set_of_lags[1] - set_of_lags[2])

corrDf <- data.frame(matrix(nrow=length(set_of_lags), ncol=2))
colnames(corrDf) <- c('lags', 'corr')
corrDf$lags <- set_of_lags
corrDf$corr <- NA

avg_corr_est <- list()


# Construct a df with columns as lags
acfDf <- data.frame(matrix(NA, nrow=length(aa_par), ncol=length(set_of_lags)))
colnames(acfDf) <- 1:length(set_of_lags)
# loop over each list in the list of lists
for(i in 1:length(aa_par)) {
  # get the acf
  acfDf[i, 1:length(aa_par[[i]]$lags)] <- aa_par[[i]]$estCorr$acf
}

avgAcfDf <- colMeans(acfDf, na.rm=T)
endT <- Sys.time()
print(endT - startT)

plot(set_of_lags, avgAcfDf / max(avgAcfDf), type='l')