library(gsignal)
library(RandomFields)
library(plot3D)
library(evmix)
library(tsqn)
library(pracma)
library(glue)
library(doParallel)
library(gstat)
registerDoParallel(detectCores() - 1)

model <- RMcauchy(gamma=0.2, var=1, scale=1)
x <- seq(-10, 50, length.out=3001)
xLen <- length(x)
y <- x
xRes <- abs(x[1] - x[2])
yRes <- abs(y[1] - y[2])
cbind(xRes, yRes)
plot(model, xlim=c(0, 5.2))

# Construct new x sets variables to remove the [-10, 0) and (40, 50] values.
zeroIdx <- which(x == 0)
fortyIdx <- which(x == 40)
belowX <- x[1:(zeroIdx-1)]
aboveX <- x[(fortyIdx + 1):length(x)]
outsideX <- c(belowX, aboveX)
xNew <- seq(0, 40, length.out=2001)
zero_forty_x <- x[zeroIdx:fortyIdx]
y <- rep(0, length(zero_forty_x))
# 750 for [0, 25], 256 for [0, 5.1]
# upperX <- length(xNew) - 750
upperX <- 256

# Estimator functions
# Yaglom ones
{
  kern_gaussian_yaglom <- function(tau, theta) {
    if(theta <= 0) {
      return("theta must be > 0")
    }
    return( exp(-(tau / theta)^2) )
  }
  
  kern_wave_yaglom <- function(tau, theta) {
    if(theta == 0) { return(1) }
    if(tau == 0) { return(1) }
    return( (theta/tau) * sin(tau/theta)  )
  }
  
  kern_rational_yaglom <- function(tau, theta) {
    return( 1 - (tau^2 / (tau^2 + theta)) )
  }
  
  # Estimator (1.68) (page 68) B_{T}^{*}(\tau) = \frac{1}{T - \tau} \sum_{t=1}^{T - \tau} x(t + \tau})x(t)
  yaglom_1.68 <- function(X, limT, tau, meanX) {
    vals <- c()
    if((limT - tau) <= 0) {
      return(0)
    }
    
    vals <- sapply(seq(1, limT - tau, by=1), function(t) (X[t+tau] - meanX)*(X[t] - meanX))
    return(sum(vals) / (limT - tau))
  }
  
  # Estimator (1.75) (page 74) B^{**}_{T}(\tau) = \frac{1}{T} \sum_{t=1}^{T - \tau} x(t + \tau)x(t)
  yaglom_1.75 <- function(X, limT, tau, meanX) {
    vals <- c()
    if((limT - tau) <= 0) {
      return(0)
    }
    vals <- sapply(seq(1, limT - tau, by=1), function(t) (X[t+tau] - meanX)*(X[t] - meanX))
    return(sum(vals) / limT)
  }
  
  
  # Estimator (1.76) (page 76) B_{T}^{(a)} = a_{T}(\tau) B^{**}_{T}(\tau)
  yaglom_1.76 <- function(X, limT, tau, trueTau, k_t, meanX) {
    vals <- c()
    if((limT - tau) <= 0) {
      return(0)
    }
    if(length(trueTau) == 0) { trueTau=0 }
    vals <- sapply(seq(1, limT - tau, by=1), function(t) (X[t+tau] - meanX)*(X[t] - meanX))
    return((sum(vals) / limT) * kern_rational_yaglom(tau, k_t))
  }
}

# Hall
{
  kern_gaussian_hall <- function(tau, theta) {
    if(theta <= 0) {
      return("theta must be > 0")
    }
    # Make its total area equal to 1 (i.e. a probably density)
    return( exp(-(tau^2) / theta) / (sqrt(theta) * sqrt(pi)) )
  }
  
  kern_rational_hall <- function(tau, theta) {
    return( (1 - (tau^2 / (tau^2 + theta))) / (pi*sqrt(theta)) )
  }
  
  acf_at_t2 <- function(X, meanX, t, x, X_ij_mat) {
    numerators <- c()
    denominators <- c()
    
    for(i in 1:length(x)) {
      tij <- x[i] - x
      t_tij <- t - tij
      X_ij <- X_ij_mat[i, ]
      K_ij <- kern_gaussian_hall(t_tij, 0.2)
      numerator <- K_ij * X_ij
      denominator <- K_ij
      
      numerators <- c(numerators, sum(numerator))
      denominators <- c(denominators, sum(denominator))
    }
    
    return( sum(numerators) / sum(denominators) )
  }
  
  
  acf_at_t2_trunc <- function(X, meanX, t, x, T1, T2, h, X_ij_mat) {
    numerators <- c()
    denominators <- c()
    
    # Case 1: 0 <= t <= T1
    # \hat{\rho}(t)
    if(t >= 0 && t <= T1) {
      for(i in 1:length(x)) {
        tij <- x[i] - x
        t_tij <- t - tij
        X_ij <- X_ij_mat[i, ]
        K_ij <- kern_rational_hall(t_tij, h)
        numerator <- K_ij * X_ij
        denominator <- K_ij
        
        numerators <- c(numerators, sum(numerator))
        denominators <- c(denominators, sum(denominator))
      }
      return( sum(numerators) / sum(denominators) )
    }
    
    # Case 2: T1 < t <= T2
    # \hat{\rho}(T1) (T2 - t) / (T2 - T1)
    else if(T1 < t && t <= T2) {
      for(i in 1:length(x)) {
        tij <- x[i] - x
        t_tij <- T1 - tij
        X_ij <- X_ij_mat[i, ]
        K_ij <- kern_rational_hall(t_tij, h)
        numerator <- K_ij * X_ij
        denominator <- K_ij
        
        numerators <- c(numerators, sum(numerator))
        denominators <- c(denominators, sum(denominator))
      }
      
      rho_T1 <- ( sum(numerators) / sum(denominators) )
      linear_part <- (T2 - t) * (T2 - T1)^(-1)
      return( rho_T1 * linear_part )
    }
    
    else {
      return(0)
    }
  }
  
}

# Tapered
{
  window <- function(x) {
    return( (1/2) * (1 - cos(pi*x)) )
  }
  
  computeTaper_single <- function(x, rho) {
    if(x >= 0 & x < rho/2) {
      return(window(2*x / rho))
    }
    else if(x >= rho/2 & x <= 1/2) {
      return(1)
    }
    # Doing this to prevent recusive calls.
    else if(x > 1/2 & x <= 1) { 
      newX <- 1 - x
      if(newX >= 0 && newX < rho/2) {
        return(window(2*newX / rho))
      }
      else if(newX >= rho/2 && newX <= 1/2) {
        return(1)
      }
    }
    else {
      return(NaN)
    }
  }
  
  taper <- function(x, rho) {
    retTaper <- c()
    for(xx in x) {
      retTaper <- c(retTaper, computeTaper_single(xx, rho))
    }
    return(retTaper)
  }
  
  Hfunction <- function(j, n, a, rho) {
    sSeq <- 1:n
    hSeq <- taper( ((sSeq - 1/2) / n), rho)^j
    return(sum(hSeq))
    
  }
  
  compute_taper <- function(X, k, n, rho, meanX) {
    acfVals <- sapply(seq(1, n - k, by=1), function(t) (X[t+k] - meanX)*(X[t] - meanX) * (taper((t - 1/2)/n, rho) * taper((t + k - 1/2)/n, rho)))
    sumAcf <- sum(acfVals)
    return(sumAcf / Hfunction(2, n, 0, rho))
  }
}

# Splines
{
  generate_knots <- function(m) {
    knotVec <- c(0)
    for(i in 1:m) {
      knotVec <- c(knotVec, i / (m + 1))
    }
    knotVec <- c(knotVec, 1)
    return(knotVec)
  }
  
  get_tau <- function(i, p, m, kVec) {
    if((i %in% -p:-1) | (i %in% (m + 2):(m + p + 1))) {
      return(i / (m + 1))
    }
    else if(i %in% 0:(m + 1)) {
      return(kVec[i + 1])
    }
    else {
      return(NA)
    }
  }
  
  get_all_tau <- function(p, m) { 
    kVec <- generate_knots(m)
    tauVec <- list()
    for(i in -p:(m + p + 1)) {
      tauVec[[glue("{i}")]] <- get_tau(i, p, m, kVec)
    }
    return(tauVec)
  }
  
  # Cox-de Boor recursion formula as per Choi2013.pdf
  f_j_l <- function(x, j, l, m, p, tau) {
    # base case
    # l == 0
    if(l == 0) {
      tau1 <- tau[[glue("{j - p}")]]
      tau2 <- tau[[glue("{j - p + 1}")]]
      if((tau1 >= 0 & tau1 <= 1) & (tau2 >= 0 & tau2 <= 1) ) {
        constant <- (m + 1) / (x + 1)
        inner <- tau2^(x + 1) - tau1^(x + 1)
        return(constant * inner)
      }
      return(0)
    }
    else {
      constant <- (m + 1) / l
      firstInner <- f_j_l(x + 1, j, l - 1, m, p, tau)
      tau_jp <- tau[[glue("{j - p}")]]
      secondInner <- tau_jp * f_j_l(x, j, l - 1, m, p, tau)
      
      tau_jpl1 <- tau[[glue("{j - p + l + 1}")]]
      
      thirdInner <- tau_jpl1 * f_j_l(x, j + 1, l - 1, m, p, tau)
      fourthInner <- f_j_l(x + 1, j + 1, l - 1, m, p, tau)
      
      retVal <- constant * (firstInner - secondInner + thirdInner - fourthInner)
      return(retVal)
    }
  }
  
  solveSpline <- function(par, wlsDf, weights) {
    if(any(par < 0)) {
      return(10^6)
    }
    return( sum(weights * (wlsDf$estCov - ( (par[1] * wlsDf$j1) + (par[2] * wlsDf$j2) + (par[3] * wlsDf$j3) + (par[4] * wlsDf$j4) + (par[5] * wlsDf$j5) ))^2) )
  }
  
  
  m <- 2
  p <- 3
  generate_knots(m)
  taus <- get_all_tau(p, m)
  
  valsDf <- data.frame(x=xNew[1:upperX])
  for(j in 1:(m+p)) {
    jStr <- glue('j{j}') 
    valsDf[, jStr] <- rep(NA, nrow(valsDf))
  }
  
  l <- p-1
  for(i in 1:nrow(valsDf)) {
    for(j in 1:(m+p)) {
      valsDf[i, j + 1] <- f_j_l(valsDf[i, 1]^2, j, l, m, p, taus)
    }
  }
}

# Error metrics
{
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
  cauchy_kriging <- function(estY, x, upperX, randomPoints, dists_withinList, spatDf, spatGrid, z, zOrig, zeroIdx) {
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
  
  cauchy_kriging_gstat <- function(estY, lags, x, upperX, randomPoints, spatDf, spatGrid, zOrig) {
    nlsDf <- data.frame(x=lags, vals=estY)
    nlsControl <- nls.control(maxiter=2000)
    fitModel <- nls(vals ~ (1 + x^2)^(-a), data=nlsDf, start=list(a=0.3), weights=rep(1/nrow(nlsDf), nrow(nlsDf)), control=nlsControl)
    m <- vgm(1, 'Cau', 1, kappa=coef(fitModel)[[1]], nugget=0.0000001)
    zz <- c()
    for(r in randomPoints) {
      zz <- c(zz, krige(z~1, spatDf, spatGrid[r, ], model=m, maxdist=x[upperX] / 2, nmax=512, debug.level=0)$var1.pred)
    }
    return( mean( abs(zOrig[randomPoints] - zz)^2 ) )
  }
}

seqRF <- xNew[1:upperX]
trueY <- (1 + seqRF^2)^(-model@par.model$gamma)

lowerIdx <- zeroIdx - upperX+130
upperIdx <- fortyIdx + upperX-130
if(lowerIdx < 1) {
  lowerIdx <- 1
}
if(upperIdx > length(x)) {
  upperIdx <- length(x)
}

outsideXIdx <- c(lowerIdx:(zeroIdx-1), (fortyIdx+1):upperIdx)

MSPE_count <- 50
iters <- 20
resList <- list("Cs"=list("errDf"=data.frame('area'=rep(NA, iters), 'dist'=rep(NA, iters), 'det'=rep(NA, iters), 'mspe'=rep(NA, iters), 'mspe_gstat'=rep(NA, iters)), "estCorr"=list(), "estCov"=list()),
                "Css"=list("errDf"=data.frame('area'=rep(NA, iters), 'dist'=rep(NA, iters), 'det'=rep(NA, iters), 'mspe'=rep(NA, iters), 'mspe_gstat'=rep(NA, iters)), "estCorr"=list(), "estCov"=list()),
                "C(a)"=list("errDf"=data.frame('area'=rep(NA, iters), 'dist'=rep(NA, iters), 'det'=rep(NA, iters), 'mspe'=rep(NA, iters), 'mspe_gstat'=rep(NA, iters)), "estCorr"=list(), "estCov"=list()),
                "hall"=list("errDf"=data.frame('area'=rep(NA, iters), 'dist'=rep(NA, iters), 'det'=rep(NA, iters), 'mspe'=rep(NA, iters), 'mspe_gstat'=rep(NA, iters)), "estCorr"=list(), "estCov"=list()),
                "CQ"=list("errDf"=data.frame('area'=rep(NA, iters), 'dist'=rep(NA, iters), 'det'=rep(NA, iters), 'mspe'=rep(NA, iters), 'mspe_gstat'=rep(NA, iters)), "estCorr"=list(), "estCov"=list()),
                "CNa"=list("errDf"=data.frame('area'=rep(NA, iters), 'dist'=rep(NA, iters), 'det'=rep(NA, iters), 'mspe'=rep(NA, iters), 'mspe_gstat'=rep(NA, iters)), "estCorr"=list(), "estCov"=list()),
                "CB"=list("errDf"=data.frame('area'=rep(NA, iters), 'dist'=rep(NA, iters), 'det'=rep(NA, iters), 'mspe'=rep(NA, iters), 'mspe_gstat'=rep(NA, iters)), "estCorr"=list(), "estCov"=list()),
                "areaRanks"=list(), "distRanks"=list(), "detRanks"=list(), "mspeRanks"=list(), 'mspeGstatRanks'=list()
)
# Total time
startT <- Sys.time()
fe <- foreach(it=1:iters, .packages=c("RandomFields", "gstat", "gsignal", "tsqn", "pracma", "minpack.lm", "glue", "dfoptim"), export=ls(globalenv())) %dopar% {
  resListInner <- list("Cs"=list("errDf"=data.frame('area'=NA, 'dist'=NA, 'det'=NA, 'mspe'=NA, 'mspe_gstat'=NA), "estCorr"=c(), "estCov"=c()),
                       "Css"=list("errDf"=data.frame('area'=NA, 'dist'=NA, 'det'=NA,'mspe'=NA, 'mspe_gstat'=NA), "estCorr"=c(), "estCov"=c()),
                       "C(a)"=list("errDf"=data.frame('area'=NA, 'dist'=NA, 'det'=NA, 'mspe'=NA, 'mspe_gstat'=NA), "estCorr"=c(), "estCov"=c()),
                       "hall"=list("errDf"=data.frame('area'=NA, 'dist'=NA, 'det'=NA, 'mspe'=NA, 'mspe_gstat'=NA), "estCorr"=c(), "estCov"=c()),
                       "CQ"=list("errDf"=data.frame('area'=NA, 'dist'=NA, 'det'=NA, 'mspe'=NA, 'mspe_gstat'=NA), "estCorr"=c(), "estCov"=c()),
                       "CNa"=list("errDf"=data.frame('area'=NA, 'dist'=NA, 'det'=NA, 'mspe'=NA, 'mspe_gstat'=NA), "estCorr"=c(), "estCov"=c()),
                       "CB"=list("errDf"=data.frame('area'=NA, 'dist'=NA, 'det'=NA, 'mspe'=NA, 'mspe_gstat'=NA), "estCorr"=c(), "estCov"=c()),
                       'areaRank'=c(), 'distRank'=c(), 'detRank'=c(), 'mspeRank'=c(), 'mspeGstatRank'=c()
  )
  ssT <- Sys.time()
  sim <- RFsimulate(model, x=x)
  z <- RFspDataFrame2conventional(sim)$data
  # Usable area should be [0, 40], not [-10, 50]
  z <- z / sd(z)
  zOrig <- z
  z <- z[zeroIdx:fortyIdx]
  resListInner[["z"]] <- z
  resListInner[["zOrig"]] <- zOrig
  
  zMean <- mean(z)
  zMean
  
  # Random points
  # Points outside the domain that are know to use, but not to the estimation process.
  randomPoints <- sample(outsideXIdx, size=MSPE_count, replace=F, prob=rep(1/length(outsideXIdx), length(outsideXIdx)))
  resListInner[["randomPoints_outsideX"]] <- randomPoints
  
  # Construct spatial dataframe for Kriging
  spatDf <- data.frame(x=zero_forty_x, y=y, z=z)
  rownames(spatDf) <- 1:nrow(spatDf)
  coordinates(spatDf) <- ~x+y
  proj4string(spatDf) <- CRS("NA")
  
  spatGrid <- data.frame(x=x, y=rep(0, length(x)))
  rownames(spatGrid) <- 1:nrow(spatGrid)
  coordinates(spatGrid) <- ~x+y
  proj4string(spatGrid) <- CRS("NA")
  
  dists_withinList <- list()
  for(i in 1:length(randomPoints)) {
    dists_withinList[[i]] <- compute_dists(seqRF, upperX, randomPoints[i], spatDf, spatGrid)
  }
  
  # C*: (1.68)
  vals <- c()
  for(i in 0:(upperX - 1)) {
    vals <-c(vals, yaglom_1.68(z, length(z), i, zMean))
  }
  
  estY <- vals[1:length(seqRF)] / vals[1]
  yaglom1_area <- areaBetween(seqRF, trueY, estY)
  resListInner[["Cs"]][["errDf"]][1, 'area'] <- yaglom1_area
  yaglom1_dist <- maxDist(trueY, estY)
  yaglom1_det <- detDist(trueY, estY)
  resListInner[["Cs"]][["errDf"]][1, 'dist'] <- yaglom1_dist
  resListInner[["Cs"]][["errDf"]][1, 'det'] <- yaglom1_det
  resListInner[["Cs"]][["errDf"]][1, 'mspe'] <- cauchy_kriging(vals, seqRF, upperX, randomPoints, dists_withinList, spatDf, spatGrid, z, zOrig, zeroIdx)
  resListInner[["Cs"]][["errDf"]][1, 'mspe_gstat'] <- cauchy_kriging_gstat(vals, lags=xNew[1:upperX], xNew, upperX, randomPoints, spatDf, spatGrid, zOrig)
  
  resListInner[["Cs"]][["estCov"]] <- vals
  resListInner[["Cs"]][["estCorr"]] <- vals / vals[1]
  print("DONE C*")
  
  # C**: 1.75
  vals2 <- c()
  zMean <- mean(z)
  for(i in 0:(upperX - 1)) {
    vals2 <- c(vals2, yaglom_1.75(z, length(z), i, zMean))
  }
  
  estY <- vals2[1:length(seqRF)] / vals2[1]
  yaglom2_area <- areaBetween(seqRF, trueY, estY)
  resListInner[["Css"]][["errDf"]][1, 'area'] <- yaglom2_area
  yaglom2_dist <- maxDist(trueY, estY)
  yaglom2_det <- detDist(trueY, estY)
  resListInner[["Css"]][["errDf"]][1, 'dist'] <- yaglom2_dist
  resListInner[["Css"]][["errDf"]][1, 'det'] <- yaglom2_det
  resListInner[["Css"]][["errDf"]][1, 'mspe'] <- cauchy_kriging(vals2, seqRF, upperX, randomPoints, dists_withinList, spatDf, spatGrid, z, zOrig, zeroIdx)
  resListInner[["Css"]][["errDf"]][1, 'mspe_gstat'] <- cauchy_kriging_gstat(vals2, lags=xNew[1:upperX], xNew, upperX, randomPoints, spatDf, spatGrid, zOrig)
  
  resListInner[["Css"]][["estCov"]] <- vals2
  resListInner[["Css"]][["estCorr"]] <- vals2 / vals2[1]
  print("DONE C**")
  
  # C(a): 1.76
  vals3 <- c()
  zMean <- mean(z)
  for(i in 0:(upperX - 1)) {
    vals3 <- c(vals3, yaglom_1.76(z, length(z), i, xNew[i+1],  1000*length(z), zMean))
  }
  
  estY <- vals3[1:length(seqRF)] / vals3[1]
  yaglom3_area <- areaBetween(seqRF, trueY, estY)
  resListInner[["C(a)"]][["errDf"]][1, 'area'] <- yaglom3_area
  yaglom3_dist <- maxDist(trueY, estY)
  yaglom3_det <- detDist(trueY, estY)
  resListInner[["C(a)"]][["errDf"]][1, 'dist'] <- yaglom3_dist
  resListInner[["C(a)"]][["errDf"]][1, 'det'] <- yaglom3_det
  resListInner[["C(a)"]][["errDf"]][1, 'mspe'] <- cauchy_kriging(vals3, seqRF, upperX, randomPoints, dists_withinList, spatDf, spatGrid, z, zOrig, zeroIdx)
  resListInner[["C(a)"]][["errDf"]][1, 'mspe_gstat'] <- cauchy_kriging_gstat(vals3, lags=xNew[1:upperX], xNew, upperX, randomPoints, spatDf, spatGrid, zOrig)
  
  resListInner[["C(a)"]][["estCov"]] <- vals3
  resListInner[["C(a)"]][["estCorr"]] <- vals3 / vals3[1]
  
  print("DONE C(a)")
  
  #===== Hall
  # Compute X_ijs
  Xij <- matrix(nrow=length(xNew), ncol=length(xNew))
  for(i in 1:length(xNew)) {
    Xij[i, ] <- sapply(1:length(xNew), function(Xj) ( (z[i] - zMean) * (z[Xj] - zMean)))
  }
  
  T1 <- 41
  T2 <- 42
  h <- 0.01
  tVals <- xNew[1:upperX]
  vals_hall <- c()
  sT <- Sys.time()
  for(i in 1:length(tVals)) {
    vals_hall <- c(vals_hall, acf_at_t2_trunc(z, zMean, tVals[i], xNew, T1, T2, h, Xij))
  }
  print(Sys.time() - sT)
  
  vals_hall_2 <- vals_hall
  # Perform a Fourier transform.
  vals_hall_2_dct <- dct(vals_hall_2)
  
  # Find the first negative value
  firstMin <- which(vals_hall_2_dct[-1] < 0)[1] + 1
  if(is.na(firstMin)) {
    print("All values are > 0!")
  }
  vals_hall_2_dct[firstMin:length(vals_hall_2_dct)] <- 0
  
  # Invert FT
  vals_hall_2_idct <- idct(vals_hall_2_dct)
  
  estY <- vals_hall_2_idct / vals_hall_2_idct[1]
  hall_area <- areaBetween(seqRF, trueY, estY)
  resListInner[["hall"]][["errDf"]][1, 'area'] <- hall_area
  hall_dist <- maxDist(trueY, estY)
  hall_det <- detDist(trueY, estY)
  resListInner[["hall"]][["errDf"]][1, 'dist'] <- hall_dist
  resListInner[["hall"]][["errDf"]][1, 'det'] <- hall_det
  resListInner[["hall"]][["errDf"]][1, 'mspe'] <- cauchy_kriging(vals_hall_2_idct, seqRF, upperX, randomPoints, dists_withinList, spatDf, spatGrid, z, zOrig, zeroIdx)
  resListInner[["hall"]][["errDf"]][1, 'mspe_gstat'] <- cauchy_kriging_gstat(vals_hall_2_idct, lags=xNew[1:upperX], xNew, upperX, randomPoints, spatDf, spatGrid, zOrig)
  
  resListInner[["hall"]][["estCov"]] <- vals_hall_2_idct
  resListInner[["hall"]][["estCorr"]] <- vals_hall_2_idct / vals_hall_2_idct[1]
  
  print("DONE hall")
  
  # Genton Qn covariance
  rob_est <- robacf(z, type='covariance', lag.max=upperX, plot=F)
  rob_est_vals <- rob_est$acf[, , 1]
  rob_est_cor_vals <- robacf(z, type='correlation', lag.max=upperX, plot=F)$acf[ , , 1]
  
  estY <- rob_est_cor_vals
  cq_area <- areaBetween(seqRF, trueY, estY)
  resListInner[["CQ"]][["errDf"]][1, 'area'] <- cq_area
  cq_dist <- maxDist(trueY, estY)
  cq_det <- detDist(trueY, estY)
  resListInner[["CQ"]][["errDf"]][1, 'dist'] <- cq_dist
  resListInner[["CQ"]][["errDf"]][1, 'det'] <- cq_det
  resListInner[["CQ"]][["errDf"]][1, 'mspe'] <- cauchy_kriging(rob_est_vals, seqRF, upperX, randomPoints, dists_withinList, spatDf, spatGrid, z, zOrig, zeroIdx)
  resListInner[["CQ"]][["errDf"]][1, 'mspe_gstat'] <- cauchy_kriging_gstat(rob_est_vals, lags=xNew[1:upperX], xNew, upperX, randomPoints, spatDf, spatGrid, zOrig)
  
  resListInner[["CQ"]][["estCov"]] <- rob_est_vals
  resListInner[["CQ"]][["estCorr"]] <- rob_est_cor_vals
  print("DONE CQ")
  
  # Tapered
  sT <- Sys.time()
  tapered <- c()
  zMean <- mean(z)
  for(i in 1:upperX) {
    tapered <- c(tapered, compute_taper(z, i-1, length(z), 1, zMean))
  }
  print(Sys.time() - sT)
  
  estY <- tapered / tapered[1]
  tapered_area <- areaBetween(seqRF, trueY, estY)
  resListInner[["CNa"]][["errDf"]][1, 'area'] <- tapered_area
  tapered_dist <- maxDist(trueY, estY)
  tapered_det <- detDist(trueY, estY)
  resListInner[["CNa"]][["errDf"]][1, 'dist'] <- tapered_dist
  resListInner[["CNa"]][["errDf"]][1, 'det'] <- tapered_det
  resListInner[["CNa"]][["errDf"]][1, 'mspe'] <- cauchy_kriging(tapered, seqRF, upperX, randomPoints, dists_withinList, spatDf, spatGrid, z, zOrig, zeroIdx)
  resListInner[["CNa"]][["errDf"]][1, 'mspe_gstat'] <- cauchy_kriging_gstat(tapered, lags=xNew[1:upperX], xNew, upperX, randomPoints, spatDf, spatGrid, zOrig)
  
  resListInner[["CNa"]][["estCov"]] <- tapered
  resListInner[["CNa"]][["estCorr"]] <- tapered / tapered[1]
  print("DONE CNa")
  
  # Splines
  weights <- c()
  for(i in 0:(upperX - 1)) {
    weights <- c(weights, (length(z) - i) / ( (1 - vals[i + 1])^2 ))
  }
  
  wlsDf <- data.frame(lags=xNew[1:upperX], 'j1'=valsDf$j1, 'j2'=valsDf$j2, 'j3'=valsDf$j3, 'j4'=valsDf$j4, 'j5'=valsDf$j5, estCov=vals)
  testLm <- hjkb(par=c(0.5, 0.5, 0.5, 0.5, 0.5), fn=solveSpline, wlsDf=wlsDf, weights=weights, lower=c(0, 0, 0, 0, 0))
  optimalBetaLM <- testLm$par
  optimalValsLM <- c()
  for(i in 1:upperX) {
    tempSum <- 0
    for(j in 2:ncol(valsDf)) {
      tempSum <- tempSum + (optimalBetaLM[j - 1] * valsDf[i, j])
    }
    optimalValsLM <- c(optimalValsLM, tempSum)
  }
  
  splineVals <- optimalValsLM
  
  estY <- splineVals[1:length(seqRF)] / splineVals[1]
  splines_area <- areaBetween(seqRF, trueY, estY)
  resListInner[["CB"]][["errDf"]][1, 'area'] <- splines_area
  splines_dist <- maxDist(trueY, estY)
  splines_det <- detDist(trueY, estY)
  resListInner[["CB"]][["errDf"]][1, 'dist'] <- splines_dist
  resListInner[["CB"]][["errDf"]][1, 'det'] <- splines_det
  resListInner[["CB"]][["errDf"]][1, 'mspe'] <- cauchy_kriging(splineVals, seqRF, upperX, randomPoints, dists_withinList, spatDf, spatGrid, z, zOrig, zeroIdx)
  resListInner[["CB"]][["errDf"]][1, 'mspe_gstat'] <- cauchy_kriging_gstat(splineVals, lags=xNew[1:upperX], xNew, upperX, randomPoints, spatDf, spatGrid, zOrig)
  
  resListInner[["CB"]][["estCov"]] <- splineVals
  resListInner[["CB"]][["estCorr"]] <- splineVals / splineVals[1]
  
  resListInner[["areaRank"]] <- rank(c(yaglom1_area, yaglom2_area, yaglom3_area, hall_area, cq_area, tapered_area, splines_area))
  resListInner[["distRank"]] <- rank(c(yaglom1_dist, yaglom2_dist, yaglom3_dist, hall_dist, cq_dist, tapered_dist, splines_dist))
  resListInner[["detRank"]] <- rank(c(yaglom1_det, yaglom2_det, yaglom3_det, hall_det, cq_det, tapered_det, splines_det))
  resListInner[["mspeRank"]] <- rank(c(resListInner[["Cs"]][["errDf"]][1, 'mspe'], resListInner[["Css"]][["errDf"]][1, 'mspe'],
                                       resListInner[["C(a)"]][["errDf"]][1, 'mspe'], resListInner[["hall"]][["errDf"]][1, 'mspe'],
                                       resListInner[["CQ"]][["errDf"]][1, 'mspe'], resListInner[["CNa"]][["errDf"]][1, 'mspe'], resListInner[["CB"]][["errDf"]][1, 'mspe']))
  resListInner[["mspeGstatRank"]] <- rank(c(resListInner[["Cs"]][["errDf"]][1, 'mspe_gstat'], resListInner[["Css"]][["errDf"]][1, 'mspe_gstat'],
                                            resListInner[["C(a)"]][["errDf"]][1, 'mspe_gstat'], resListInner[["hall"]][["errDf"]][1, 'mspe_gstat'],
                                            resListInner[["CQ"]][["errDf"]][1, 'mspe_gstat'], resListInner[["CNa"]][["errDf"]][1, 'mspe_gstat'], resListInner[["CB"]][["errDf"]][1, 'mspe_gstat']))
  
  
  print("Time taken for loop:")
  print(Sys.time() - ssT)
  return( resListInner )
}

print(glue("Total time for ", iters, " iterations: ", (Sys.time() - startT)))

# Convert to resList format.
for(i in 1:length(fe)) {
  for(n in names(resList)) {
    resList[[n]][["errDf"]][i, ] <- fe[[i]][[n]][["errDf"]]
    resList[[n]][["estCorr"]][[i]] <- fe[[i]][[n]][["estCorr"]]
    resList[[n]][["estCov"]][[i]] <- fe[[i]][[n]][["estCov"]]
  }
  resList[["areaRanks"]][[i]] <- fe[[i]][["areaRank"]]
  resList[["distRanks"]][[i]] <- fe[[i]][["distRank"]]
  resList[["detRanks"]][[i]] <- fe[[i]][["detRank"]]
  resList[["mspeRanks"]][[i]] <- fe[[i]][["mspeRank"]]
  resList[["mspeGstatRanks"]][[i]] <- fe[[i]][["mspeGstatRank"]]
}

# Find the minimum/maximum of each ACF to act as bands, plot them, and save the plots.
count <- 1
colours <- c(2, 3, 4, 6, 7, 8, 13)
adjustedColours <- sapply(colours, function(c) adjustcolor(c, alpha=0.2))
nEsts <- 7
estNames <- names(resList)[-((nEsts + 1):length(resList))]
for(n in estNames) {
  minACF_C_star <- c()
  maxACF_C_star <- c()
  corrList <- resList[[n]][["estCorr"]]
  for(i in 1:upperX) {
    minVal <- Inf
    maxVal <- -Inf
    for(j in 1:length(corrList)) {
      if(corrList[[j]][i] < minVal) {
        minVal <- corrList[[j]][i]
      }
      if(corrList[[j]][i] > maxVal) {
        maxVal <- corrList[[j]][i]
      }
    }
    minACF_C_star <- c(minACF_C_star, minVal)
    maxACF_C_star <- c(maxACF_C_star, maxVal)
  } 
  
  # PLOT
  plot(xNew[1:upperX], (1 + xNew[1:upperX]^2)^(-model@par.model$gamma), col=1, lty=1, lwd=2, type='l', xlab=expression(h), ylab=expression(rho(h)), ylim=c(-1, 1))
  lines(xNew[1:upperX], minACF_C_star, lwd=2, lty=1, col=colours[count])
  lines(xNew[1:upperX], maxACF_C_star, lwd=2, lty=1, col=colours[count])
  
  polygon(c(xNew[1:upperX], rev(xNew[1:upperX])), c(minACF_C_star, rev(maxACF_C_star)), col=adjustedColours[count])
  count <- count + 1
}

# All bands.
# Basically the above code, but only calls plot once.
plotName <- glue("cauchy_corr/bands.pdf")
pdf(plotName, width=10, height=9.6)
par(mar=c(4,4.5,0.25,0.25)+.1)
plot(xNew[1:upperX], (1 + xNew[1:upperX]^2)^(-model@par.model$gamma), col=1, lty=1, lwd=2, type='l', xlab=expression(h), ylab=expression(hat(rho)(h)*' bands'), ylim=c(-0.5, 1))

count <- 1
colours <- c(2, 3, 4, 6, 7, 8, 13)
for(n in estNames) {
  minACF_C_star <- c()
  maxACF_C_star <- c()
  corrList <- resList[[n]][["estCorr"]]
  for(i in 1:upperX) {
    minVal <- Inf
    maxVal <- -Inf
    for(j in 1:length(corrList)) {
      if(corrList[[j]][i] < minVal) {
        minVal <- corrList[[j]][i]
      }
      if(corrList[[j]][i] > maxVal) {
        maxVal <- corrList[[j]][i]
      }
    }
    minACF_C_star <- c(minACF_C_star, minVal)
    maxACF_C_star <- c(maxACF_C_star, maxVal)
  } 
  
  # PLOT
  lines(xNew[1:upperX], minACF_C_star, lwd=2, lty=1, col=colours[count])
  lines(xNew[1:upperX], maxACF_C_star, lwd=2, lty=1, col=colours[count])
  
  polygon(c(xNew[1:upperX], rev(xNew[1:upperX])), c(minACF_C_star, rev(maxACF_C_star)), col=adjustcolor(colours[count], alpha=0.2), lty=0)
  count <- count + 1
}
adjustedColours <- sapply(colours, function(c) adjustcolor(c, alpha=0.2))
legend('bottomleft', c('True', expression('C'^'*'*'(h)'), expression('C'^'**'*'(h)'), expression('C'^'(a)'*'(h)'),
                       expression(tilde('C')*'(h)'), expression(hat('C')[Q]*'(h)'), expression(hat('C')[N]^'a'*'(h)'), expression(hat('C')^'B'*('h'))),
       col=c(1, colours), lty=c(1, rep(NA, 7)), lwd=c(2, rep(NA, 7)), density=c(0, rep(NA, 7)), fill=c(NA, adjustedColours), border=c(NA, colours), y.intersp=1.2, cex=1.2)

dev.off()

# Compute average errors
nEsts <- 7
estNames <- names(resList)[-((nEsts + 1):length(resList))]
resDf <- data.frame('area'=rep(NA, nEsts), 'dist'=rep(NA, nEsts), 'det'=rep(NA, nEsts), 'mspe'=rep(NA, nEsts), 'mspe_gstat'=rep(NA, nEsts))
rownames(resDf) <- estNames

for(n in estNames) {
  resDf[n, ] <- colMeans(resList[[n]][["errDf"]])
}

rownames(resDf) <- c('C^{*}', 'C^{**}', 'C^{(a)}', '\\widehat{\\rho}', '\\widehat{C}_{Q}', 'C_{N}^{a}', 'C^{B}')
resDf

# Ranks
rankNames <- names(resList)[(nEsts + 1):length(resList)]
rankDf <- data.frame('areaRanks'=rep(NA, nEsts), 'distRanks'=rep(NA, nEsts), 'detRanks'=rep(NA, nEsts), 'mspeRanks'=rep(NA, nEsts), 'mspeGstatRanks'=rep(NA, nEsts))
for(r in rankNames) {
  rankDf[, r] <- rowMeans(matrix(unlist(resList[[r]]), nrow=nEsts, ncol=iters))
}

rownames(rankDf) <- c('C^{*}', 'C^{**}', 'C^{(a)}', '\\widehat{\\rho}', '\\widehat{C}_{Q}', 'C_{N}^{a}', 'C^{B}')
rankDf

kableExtra::kable(resDf, format='latex', align='c', booktabs=F, escape=F)

resDf2 <- unlist(resDf)
resDf2[which(resDf2 < sqrt(.Machine$double.eps))] <- 0
resDf2 <- matrix(resDf2, nrow=7, ncol=5)
rownames(resDf2) <- c('C^{*}', 'C^{**}', 'C^{(a)}', '\\widehat{\\rho}', '\\widehat{C}_{Q}', 'C_{N}^{a}', 'C^{B}')
colnames(resDf2) <- colnames(resDf)
resDf2 <- as.data.frame(resDf2)
knitr::kable(resDf2, format='latex', align='c', booktabs=F, escape=F, digits=32)

doParallel::stopImplicitCluster()
