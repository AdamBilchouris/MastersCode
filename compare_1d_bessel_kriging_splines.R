load("1d_bessel_n10_50_3001points_ests_NEWER_BESSEL_10.RData")

areaBetween <- function(xVals, trueY, estY) {
  if(length(xVals) != length(trueY)) {
    return("length of x != length of true y")
  }
  
  if(length(xVals) != length(estY)) {
    return("length of x != length of estimated y")
  }
  
  yDiff <- abs(trueY - estY)
  # (b-a) / n 
  weight <- (xVals[length(xVals)] - xVals[1]) / length(xVals)
  # \sigma_{i=1}^{n} f(x_{i}) \delta x
  return(sum(yDiff * weight))
}

maxDist <- function(trueY, estY) {
  if(length(trueY) != length(estY)) {
    return("length of true y != length of estimated y")
  }
  maxDist <- -1000
  for(i in 1:length(trueY)) {
    d <- abs(trueY[i] - estY[i])
    if(d > maxDist) {
      maxDist <- d
    }
  }
  return(maxDist)
}

create_cyclic_matrix <- function(vec) {
  n <- length(vec)
  mat <- matrix(NA, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      mat[i, j] <- (vec[1+abs(i-j)])
    }}
  return(mat)
}

detDist <- function(trueY, estY) {
  if(length(trueY) != length(estY)) {
    return("length of true y != length of estimated y")
  }
  det1 <- norm(create_cyclic_matrix((trueY-estY)),type = "2")
  detDist <- det1
  return(detDist)
}

compute_dists <- function(x, upperX, randomPoint, spatDf, spatGrid) {
  dists <- spDists(coordinates(spatDf), coordinates(spatGrid[randomPoint, ]))
  dists <- data.frame(dist=dists, idx=1:nrow(spatDf))
  rownames(dists) <- dists$idx
  dists_within <- dists[dists$dist <= x[upperX] / 2, ]
  
  dists_within_order <- dists_within[order(dists_within$dist), ]
  
  if(nrow(dists_within_order) > 512) {
    dists_within_order <- dists_within_order[1:512, ]
  }
  
  # Compute the distances for each pair of the filtered distances.
  candidatePoints <- spatDf[dists_within_order$idx, ]
  candidateDists <- spDists(coordinates(candidatePoints))
  colnames(candidateDists) <- dists_within_order$idx
  rownames(candidateDists) <- dists_within_order$idx
  
  candidateDists[candidateDists > x[upperX]] <- 0
  
  return(list('dists_within_order'=dists_within_order, 'candidateDists'=candidateDists))
}

generate_cov_matrix_vector <- function(x, uppperX, randomPoint, spatDf, spatGrid, estY, dists_within_order, candidateDists, zeroIdx) {
  # Find the covariance values for the above matrix.
  idxDist <- as.matrix(dist(dists_within_order$idx)) + 1
  covMat2 <- matrix(estY[idxDist], ncol=ncol(idxDist))
  covMat2[is.na(covMat2)] <- 0
  covMat <- covMat2
  
  v1 <- rep(1, nrow(covMat))
  v2 <- (c(v1, 0))
  covMat <- rbind(covMat, v1)
  covMat <- cbind(covMat, v2)
  
  seqIdx <- abs((dists_within_order$idx - (randomPoint - zeroIdx + 1))) + 1
  estY_val <- estY[seqIdx]
  dists_within_order_idx <- data.frame(dists_within_order, 'xIdx'=seqIdx, 'estY_val'=estY_val)
  
  return(list('covMat'=covMat, 'covVec'=dists_within_order_idx))
}

# x = seqRF
bessel_kriging <- function(estY, x, upperX, randomPoints, dists_withinList, spatDf, spatGrid, z, zOrig, zeroIdx) {
  predVec <- rep(NA, length(randomPoints))
  
  for(i in 1:length(randomPoints)) {
    covList <- generate_cov_matrix_vector(x, upperX, randomPoints[i], spatDf, spatGrid, estY, 
                                          dists_withinList[[i]][['dists_within_order']], dists_withinList[[i]][['candidateDists']], zeroIdx)
    covMat <- covList[['covMat']]
    dists_within_order_idx <- covList[['covVec']]
    
    result_matrix_inv <- ginv(covMat)
    
    vec_1 <- dists_within_order_idx$estY_val
    vec1 <- c(vec_1, 1)
    lamvec <- result_matrix_inv %*% vec1
    estZ <- lamvec[1:length(vec_1)] %*% z[dists_within_order_idx$idx]
    (zOrig[randomPoints[i]] - estZ)^2
    predVec[i] <- (zOrig[randomPoints[i]] - estZ)^2
  }
  
  return( mean(predVec) )
}

bessel_kriging_gstat <- function(estY, lags, x, upperX, randomPoints, spatDf, spatGrid, zOrig) {
  nlsDf <- data.frame(x=lags, vals=estY)
  nlsControl <- nls.control(maxiter=400)
  fitModel <- nls(vals ~ ifelse(x==0, 1,  (2^(a) * gamma(a+1) * besselJ(x, nu=a)) / (x^(a))), data=nlsDf, start=list(a=0), weights=rep(1/nrow(nlsDf), nrow(nlsDf)), control=nlsControl)
  # Ensure it is a valid covariance function!
  if(coef(fitModel)[[1]] <= -0.5) {
    fitModel$m$setPars(-0.499)
  }
  m <- vgm(1, 'Be2', 1, kappa=coef(fitModel)[[1]] , nugget=0.0000001)
  zz <- c()
  for(r in randomPoints) {
    zz <- c(zz, krige(z~1, spatDf, spatGrid[r, ], model=m, maxdist=x[upperX] / 2, nmax=512, debug.level=0)$var1.pred)
  }
  return( mean( abs(zOrig[randomPoints] - zz)^2 ) )
}

seqRF <- xNew[1:256]
trueY <- c(1, ((2^(model@par.model$nu) * gamma(model@par.model$nu+1) * besselJ(seqRF, nu=model@par.model$nu)) / (seqRF^(model@par.model$nu)))[-1])

y <- rep(0, length(zero_forty_x))
spatDf <- data.frame(x=zero_forty_x, y=y, z=z)
rownames(spatDf) <- 1:nrow(spatDf)
coordinates(spatDf) <- ~x+y
proj4string(spatDf) <- CRS("NA")

spatGrid <- data.frame(x=x, y=rep(0, length(x)))
rownames(spatGrid) <- 1:nrow(spatGrid)
coordinates(spatGrid) <- ~x+y
proj4string(spatGrid) <- CRS("NA")

seqRandom <- outsideXIdx
randomPoints <- sample(seqRandom, size=MSPE_count, replace=F, prob=rep(1/length(seqRandom), length(seqRandom)))

dists_withinList <- list()
for(i in 1:length(randomPoints)) {
  dists_withinList[[i]] <- compute_dists(seqRF, upperX, randomPoints[i], spatDf, spatGrid)
}

estY <- vals[1:length(seqRF)] / vals[1]
yaglom1_area <- areaBetween(seqRF, trueY, estY)
yaglom1_area
yaglom1_dist <- maxDist(trueY, estY)
yaglom1_dist
yaglom1_kriging <- bessel_kriging(vals, seqRF, upperX, randomPoints, dists_withinList, spatDf, spatGrid, z, zOrig, zeroIdx)
yaglom1_kriging
yaglom1_kriging_gstat <- bessel_kriging_gstat(vals, lags=xNew[1:256], xNew, 256, randomPoints, spatDf, spatGrid, zOrig)
yaglom1_kriging_gstat
yaglom1_det <- detDist(trueY, estY)
yaglom1_det

estY <- vals2[1:length(seqRF)] / vals2[1]
yaglom2_area <- areaBetween(seqRF, trueY, estY)
yaglom2_area
yaglom2_dist <- maxDist(trueY, estY)
yaglom2_dist
yaglom2_kriging <- bessel_kriging(vals2, seqRF, upperX, randomPoints, dists_withinList, spatDf, spatGrid, z, zOrig, zeroIdx)
yaglom2_kriging
yaglom2_kriging_gstat <- bessel_kriging_gstat(vals2, lags=xNew[1:256], xNew, 256, randomPoints, spatDf, spatGrid, zOrig)
yaglom2_kriging_gstat
yaglom2_det <- detDist(trueY, estY)
yaglom2_det


estY <- vals3[1:length(seqRF)] / vals3[1]
yaglom3_area <- areaBetween(seqRF, trueY, estY)
yaglom3_area
yaglom3_dist <- maxDist(trueY, estY)
yaglom3_dist
yaglom3_kriging <- bessel_kriging(vals3, seqRF, upperX, randomPoints, dists_withinList, spatDf, spatGrid, z, zOrig, zeroIdx)
yaglom3_kriging
yaglom3_kriging_gstat <- bessel_kriging_gstat(vals3, lags=xNew[1:256], xNew, 256, randomPoints, spatDf, spatGrid, zOrig)
yaglom3_kriging_gstat
yaglom3_det <- detDist(trueY, estY)
yaglom3_det


estY <- vals_hall_2_idct[1:length(seqRF)] / vals_hall_2_idct[1]
hall_area <- areaBetween(seqRF, trueY, estY)
hall_area
hall_dist <- maxDist(trueY, estY)
hall_dist
hall_kriging <- bessel_kriging(vals_hall_2_idct, seqRF, upperX, randomPoints, dists_withinList, spatDf, spatGrid, z, zOrig, zeroIdx)
hall_kriging
hall_kriging_gstat <- bessel_kriging_gstat(vals_hall_2_idct, lags=xNew[1:256], xNew, 256, randomPoints, spatDf, spatGrid, zOrig)
hall_kriging_gstat
hall_det <- detDist(trueY, estY)
hall_det


estY <- rob_est_cor_vals[1:length(seqRF)]
rob_area <- areaBetween(seqRF, trueY, estY)
rob_area
rob_dist <- maxDist(trueY, estY)
rob_dist
rob_kriging <- bessel_kriging(rob_est_vals, seqRF, upperX, randomPoints, dists_withinList, spatDf, spatGrid, z, zOrig, zeroIdx)
rob_kriging
rob_kriging_gstat <- bessel_kriging_gstat(rob_est_vals, lags=xNew[1:256], xNew, 256, randomPoints, spatDf, spatGrid, zOrig)
rob_kriging_gstat
rob_det <- detDist(trueY, estY)
rob_det


estY <- tapered[1:length(seqRF)] / tapered[1]
tapered_area <- areaBetween(seqRF, trueY, estY)
tapered_area
tapered_dist <- maxDist(trueY, estY)
tapered_dist
tapered_kriging <- bessel_kriging(tapered, seqRF, upperX, randomPoints, dists_withinList, spatDf, spatGrid, z, zOrig, zeroIdx)
tapered_kriging
tapered_kriging_gstat <- bessel_kriging_gstat(tapered, lags=xNew[1:256], xNew, 256, randomPoints, spatDf, spatGrid, zOrig)
tapered_kriging_gstat
tapered_det <- detDist(trueY, estY)
tapered_det


estY <- splineVals[1:length(seqRF)] / splineVals[1]
splines_area <- areaBetween(seqRF, trueY, estY)
splines_area
splines_dist <- maxDist(trueY, estY)
splines_dist
splines_kriging <- bessel_kriging(splineVals, seqRF, upperX, randomPoints, dists_withinList, spatDf, spatGrid, z, zOrig, zeroIdx)
splines_kriging
splines_kriging_gstat <- bessel_kriging_gstat(splineVals, lags=xNew[1:256], xNew, 256, randomPoints, spatDf, spatGrid, zOrig)
splines_kriging_gstat
splines_det <- detDist(trueY, estY)
splines_det


ests <- data.frame('area'=c(yaglom1_area, yaglom2_area, yaglom3_area, hall_area, rob_area, tapered_area, splines_area),
                   'max_dist'=c(yaglom1_dist, yaglom2_dist, yaglom3_dist, hall_dist, rob_dist, tapered_dist, splines_dist),
                   'det'=c(yaglom1_det, yaglom2_det, yaglom3_det, hall_det, rob_det, tapered_det, splines_det),
                   'mspe'=c(yaglom1_kriging, yaglom2_kriging, yaglom3_kriging, hall_kriging, rob_kriging, tapered_kriging, splines_kriging),
                   'mspeGstat'=c(yaglom1_kriging_gstat, yaglom2_kriging_gstat, yaglom3_kriging_gstat, hall_kriging_gstat, rob_kriging_gstat, tapered_kriging_gstat, splines_kriging_gstat))
rownames(ests) <- c('C^{*}', 'C^{**}', 'C^{(a)}', '\\widetidle{C}', '\\widehat{C}_{Q}', '\\widehat{C}_{N}^{a}', '\\widehat{C}^{B}')
ests

kableExtra::kable(ests, format='latex', align='c', booktabs=F, escape=F)
