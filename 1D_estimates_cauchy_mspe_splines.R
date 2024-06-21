library(gsignal)
library(RandomFields)
library(plot3D)
library(evmix)
library(tsqn)
library(dfoptim)

load("1d_cauchy_gamma_0.2_n10_50_3001points_origVar_1.017_newVar_0.984.RData")

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

outsideXIdx <- c(lowerIdx:(zeroIdx-1), (fortyIdx+1):upperIdx)

# Run if you aren't loading from an Rdata file.
while((origVar < 0.98 | origVar > 1.02) | (newVar < 0.98 | newVar > 1.02)) {
  sim <- RFsimulate(model, x=x)
  
  z <- RFspDataFrame2conventional(sim)$data
  origVar <- var(z)
  # Usable area should be [0, 40], not [-10, 50]
  zOrig <- z
  z <- z[zeroIdx:fortyIdx]
  newVar <- var(z)
}
zMean <- mean(z)
zMean

# For [0, 25]:
# upperX <- length(x) - 750
# For [0, 5.1]:
upperX <- 256
lowerIdx <- zeroIdx - upperX+130
upperIdx <- fortyIdx + upperX-130
if(lowerIdx < 1) {
  lowerIdx <- 1
}
if(upperIdx > length(x)) {
  upperIdx <- length(x)
}

kern_gaussian_yaglom <- function(tau, theta) {
  if(theta <= 0) {
    return("theta must be > 0")
  }
  return( exp(-(tau / theta)^2) )
}

kern_rational_yaglom <- function(tau, theta) {
  return( 1 - (tau^2 / (tau^2 + theta)) )
}

kern_wave_yaglom <- function(tau, theta) {
  if(theta == 0) { return(1) }
  return( (theta/tau) * sin(tau/theta)  )
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


vals <- c()
zMean <- mean(z)
for(i in 0:255) {
  vals <-c(vals, yaglom_1.68(z, length(z), i, zMean))
}


xNew[256]
plot(xNew[1:256], vals, type='o', ylim=c(-0.1, 1))
lines(xNew[1:256], (1 + xNew[1:256]^2)^(-model@par.model$gamma), col=2, lty=2, lwd=2, type='l', xlab=expression(tau), ylab=expression(rho(tau)), ylim=c(-0.1, 1))


# Estimator (1.75) (page 74) B^{**}_{T}(\tau) = \frac{1}{T} \sum_{t=1}^{T - \tau} x(t + \tau)x(t)
yaglom_1.75 <- function(X, limT, tau, meanX) {
  vals <- c()
  if((limT - tau) <= 0) {
    return(0)
  }
  vals <- sapply(seq(1, limT - tau, by=1), function(t) (X[t+tau] - meanX)*(X[t] - meanX))
  return(sum(vals) / limT)
}

vals2 <- c()
for(i in 0:255) {
  vals2 <-c(vals2, yaglom_1.75(z, length(z), i, zMean))
}

xNew[256]
plot(xNew[1:256], vals2, type='o', ylim=c(-0.1, 1))
lines(xNew[1:256], (1 + xNew[1:256]^2)^(-model@par.model$gamma), col=2, lty=2, lwd=2, type='l', xlab=expression(tau), ylab=expression(rho(tau)), ylim=c(-0.1, 1))

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

vals3 <- c()
zMean <- mean(z)
for(i in 0:255) {
  vals3 <- c(vals3, yaglom_1.76(z, length(z), i, xNew[i+1], 1000*length(z), zMean))
}

xNew[256]
plot(xNew[1:256], vals3, type='o', ylim=c(-0.1, 1))
lines(xNew[1:256], (1 + xNew[1:256]^2)^(-model@par.model$gamma), col=2, lty=2, lwd=2, type='l', xlab=expression(tau), ylab=expression(rho(tau)), ylim=c(-0.1, 1))


#===== Hall
kern_gaussian_hall <- function(tau, theta) {
  if(theta <= 0) {
    return("theta must be > 0")
  }
  # Make its total area equal to 1 (i.e. a probably density)
  return( exp(-(tau^2) / theta) / (sqrt(theta) * sqrt(pi)) )
}

kern_wave_hall <- function(tau, theta) {
  if(tau == 0) { return( 1 / ( (pi*theta)/2 ) ) }
  return( ((theta/tau) * sin(tau/theta)) / (pi*theta / 2) )
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

# Compute X_ijs
Xij <- matrix(nrow=length(xNew), ncol=length(xNew))
for(i in 1:length(xNew)) {
  Xij[i, ] <- sapply(1:length(xNew), function(Xj) ( (z[i] - zMean) * (z[Xj] - zMean)))
}

T1 <- 6
T2 <- 7
h <- 0.01
tVals <- xNew[1:256]
vals_hall <- c()
sT <- Sys.time()
for(i in 1:length(tVals)) {
  vals_hall <- c(vals_hall, acf_at_t2_trunc(z, zMean, tVals[i], xNew, T1, T2, h, Xij))
}
print(Sys.time() - sT)

plot(tVals, vals_hall, type='l', ylim=c(-0.1, 1.1))
lines(xNew[1:256], (1 + xNew[1:256]^2)^(-model@par.model$gamma), col=2, lty=2, lwd=2, type='l', xlab=expression(tau), ylab=expression(rho(tau)), ylim=c(-0.1, 1))

vals_hall_2 <- vals_hall
# Perform a Fourier transform.
vals_hall_2_dct <- dct(vals_hall_2)
plot(vals_hall_2_dct, type='l')

# Find the first negative value
firstMin <- which(vals_hall_2_dct[-1] < 0)[1] + 1
if(is.na(firstMin)) {
  print("All values are > 0!")
}
vals_hall_2_dct[firstMin:length(vals_hall_2_dct)] <- 0
plot(vals_hall_2_dct, type='l')

# Invert FT
vals_hall_2_idct <- idct(vals_hall_2_dct)
plot(tVals, vals_hall_2_idct / vals_hall_2_idct[1], type='l')
lines(xNew[1:256], (1 + xNew[1:256]^2)^(-model@par.model$gamma), col=2, lty=2, lwd=2, type='l', xlab=expression(tau), ylab=expression(rho(tau)), ylim=c(-0.1, 1))

# Genton Qn covariance
rob_est <- robacf(z, type='covariance', lag.max=256, plot=F)
rob_est_vals <- rob_est$acf[, , 1]
rob_est_cor_vals <- robacf(z, type='correlation', lag.max=256, plot=F)$acf[,,1]

# Tapered
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

curve(taper(x, rho=1), from=0, to=1)

Hfunction <- function(j, n, a, rho) {
  sSeq <- 1:n
  hSeq <- taper( ((sSeq - 1/2) / n), rho)^j
  return(sum(hSeq))
}

Hfunction(2, 100, 0, 1)

compute_taper <- function(X, k, n, rho, meanX) {
  acfVals <- sapply(seq(1, n - k, by=1), function(t) (X[t+k] - meanX)*(X[t] - meanX) * (taper((t - 1/2)/n, rho) * taper((t + k - 1/2)/n, rho)))
  sumAcf <- sum(acfVals)
  return(sumAcf / Hfunction(2, n, 0, rho))
}


sT <- Sys.time()
tapered <- c()
zMean <- mean(z)
for(i in 1:256) {
  tapered <- c(tapered, compute_taper(z, i-1, length(z), 1, zMean))
}
print(Sys.time() - sT)

# Splines
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

solveSpline <- function(par, wlsDf, weights) {
  if(any(par < 0)) {
    return(10^6)
  }
  return( sum(weights * (wlsDf$estCov - ( (par[1] * wlsDf$j1) + (par[2] * wlsDf$j2) + (par[3] * wlsDf$j3) + (par[4] * wlsDf$j4) + (par[5] * wlsDf$j5) ))^2) )
}

weights <- c()
for(i in 0:255) {
  weights <- c(weights, (length(z) - i) / ( (1 - vals[i + 1])^2 ))
}

wlsDf <- data.frame(lags=xNew[1:upperX], 'j1'=valsDf$j1, 'j2'=valsDf$j2, 'j3'=valsDf$j3, 'j4'=valsDf$j4, 'j5'=valsDf$j5, estCov=vals)
testLm <- hjkb(par=c(0.5, 0.5, 0.5, 0.5, 0.5), fn=solveSpline, wlsDf=wlsDf, weights=weights, lower=c(0, 0, 0, 0, 0))
optimalBetaLM <- testLm$par
optimalValsLM <- c()
for(i in 1:length(xNew[1:upperX])) {
  tempSum <- 0
  for(j in 2:ncol(valsDf)) {
    tempSum <- tempSum + (optimalBetaLM[j - 1] * valsDf[i, j])
  }
  optimalValsLM <- c(optimalValsLM, tempSum)
}
splineVals <- optimalValsLM

# Plot!
# Autocorrelation
plotName <- glue("cauchy_autocorr.pdf")
pdf(plotName, width=10, height=9.6)
par(mar=c(4,4.5,0.25,0.25)+.1)
plot(xNew[1:256], (1 + xNew[1:256]^2)^(-model@par.model$gamma), col=1, lty=1, lwd=2, type='l', xlab=expression(h), ylab=expression(hat(rho)(h)), ylim=c(-0.35, 1))
lines(xNew[1:256], vals / vals[1], col=2, lty=2, lwd=2)
lines(xNew[1:256], vals2 / vals2[1], col=3, lty=3, lwd=2)
lines(xNew[1:256], vals3 / vals3[1], col=4, lty=4, lwd=2)
lines(xNew[1:256], vals_hall_2_idct / vals_hall_2_idct[1], col=6, lty=5, lwd=2)
lines(xNew[1:256], rob_est_cor_vals, col=7, lty=6, lwd=2)
lines(xNew[1:256], tapered / tapered[1], col=8, lty=8, lwd=2)
lines(xNew[1:256], splineVals / splineVals[1], col=13, lty=9, lwd=2)
legend('bottomleft', c('True', expression('C'^'*'*'(h)'), expression('C'^'**'*'(h)'), expression('C'^'(a)'*'(h)'),
                     expression(tilde('C')*'(h)'), expression(hat('C')[Q]*'(h)'), expression(hat('C')[N]^'a'*'(h)'), expression(hat('C')^'B'*('h'))),
       lty=c(1, 2, 3, 4, 5, 6, 8, 9), col=c(1, 2, 3, 4, 6, 7, 8, 13), lwd=c(2, 2, 2, 2, 2, 2, 2, 2), y.intersp=1.2, cex=1.2)
dev.off()


# Autocovariance
plotName <- glue("cauchy_autocov.pdf")
pdf(plotName, width=10, height=9.6)
par(mar=c(4,4.5,0.25,0.25)+.1)
plot(xNew[1:256], (1 + xNew[1:256]^2)^(-model@par.model$gamma), col=1, lty=1, lwd=2, type='l', xlab=expression(h), ylab=expression(hat(C)(h)), ylim=c(-0.4, 1.1))
ticks <- c(-0.4, -0.25, 0.25, 0.75)
axis(side=2, at=ticks)
lines(xNew[1:256], vals, col=2, lty=2, lwd=2)
lines(xNew[1:256], vals2, col=3, lty=3, lwd=2)
lines(xNew[1:256], vals3, col=4, lty=4, lwd=2)
lines(xNew[1:256], vals_hall_2_idct, col=6, lty=5, lwd=2)
lines(xNew[1:256], rob_est_vals, col=7, lty=6, lwd=2)
lines(xNew[1:256], tapered, col=8, lty=8, lwd=2)
lines(xNew[1:256], splineVals, col=13, lty=9, lwd=2)
legend('bottomleft', c('True', expression('C'^'*'*'(h)'), expression('C'^'**'*'(h)'), expression('C'^'(a)'*'(h)'),
                     expression(tilde('C')*'(h)'), expression(hat('C')[Q]*'(h)'), expression(hat('C')[N]^'a'*'(h)'), expression(hat('C')^'B'*('h'))),
       lty=c(1, 2, 3, 4, 5, 6, 8, 9), col=c(1, 2, 3, 4, 6, 7, 8, 13), lwd=c(2, 2, 2, 2, 2, 2, 2, 2), y.intersp=1.2, cex=1.2)
dev.off()

