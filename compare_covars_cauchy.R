rm(list=ls())
gc()
areaBetween <- function(xVals, trueY, estY) {
  if(length(xVals) != length(trueY)) {
    return("length of x != length of true y")
  }
  
  if(length(xVals) != length(estY)) {
    return("length of x != length of estimated y")
  }
  
  # Somehow?
  # This case shouldn't be reachable
  # if(length(trueY) != length(estY)) {
  #   return("length of true y != length of estimated y")
  # }
  
  # diffFun <- approxfun(xVals, trueY - estY)
  # absDiffFun <- function(x) { abs(diffFun(x)) }
  # int <- integrate(absDiffFun, min(xVals), max(xVals), rel.tol=.Machine$double.eps^.05)
  # 
  # if(int$message == 'OK') {
  #   return(list('error'=int$value, 'abs.error'=int$abs.error))
  # }
  # print("Something went wrong...")
  # print("Check the result.")
  # return(int)
  
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

load("from_my_pc/Cauchy_gamma_0.1_n10_p10_res_0.1_horizontal_vertical.RData")

seqRF <- seq(0, 6, by=0.1)

trueY <- (1 + seqRF^2)^(-0.1) # c(1, ((2^(model@par.model$nu) * gamma(model@par.model$nu+1) * besselJ(seqRF, nu=model@par.model$nu)) / (seqRF^(model@par.model$nu)))[-1])
estY <- estimatedAcf[1:length(seqRF)]

plot(seqRF, trueY, ylim=c(-0.1, 1.1))
lines(seqRF, estY, col='red', lty=2, lwd=2)

# Compute area between curves
estArea_hor_ver <- areaBetween(seqRF, trueY, estY)
estArea_hor_ver
maxDist_hor_ver <- maxDist(trueY, estY)
maxDist_hor_ver

load("from_my_pc/Cauchy_gamma_0.1_n10_p10_res_0.1_eps_0.1_cardinal_y.RData")
seqRF <- seq(0, 6, by=0.1)

trueY <- (1 + seqRF^2)^(-0.1) # c(1, ((2^(model@par.model$nu) * gamma(model@par.model$nu+1) * besselJ(seqRF, nu=model@par.model$nu)) / (seqRF^(model@par.model$nu)))[-1])
estY <- avgAcfDf / max(avgAcfDf)
estY <- estY[1:length(seqRF)]

plot(seqRF, trueY, ylim=c(-0.1, 1.1))
lines(seqRF, estY, col='red', lty=2, lwd=2)
estArea_0.1 <- areaBetween(seqRF, trueY, estY)
estArea_0.1
maxDist_0.1 <- maxDist(trueY, estY)
maxDist_0.1


load("from_my_pc/Cauchy_gamma_0.1_n10_p10_res_0.1_eps_0.2_cardinal_y.RData")
seqRF <- seq(0, 6, by=0.1)

trueY <- (1 + seqRF^2)^(-0.1) # c(1, ((2^(model@par.model$nu) * gamma(model@par.model$nu+1) * besselJ(seqRF, nu=model@par.model$nu)) / (seqRF^(model@par.model$nu)))[-1])
estY <- avgAcfDf / max(avgAcfDf)
estY <- estY[1:length(seqRF)]

plot(seqRF, trueY, ylim=c(-0.1, 1.1))
lines(seqRF, estY, col='red', lty=2, lwd=2)
estArea_0.2 <- areaBetween(seqRF, trueY, estY)
estArea_0.2
maxDist_0.2 <- maxDist(trueY, estY)
maxDist_0.2

load("from_my_pc/Cauchy_gamma_0.1_n10_p10_res_0.1_eps_0.25_cardinal_y.RData")
seqRF <- seq(0, 6, by=0.1)

trueY <- (1 + seqRF^2)^(-0.1) # c(1, ((2^(model@par.model$nu) * gamma(model@par.model$nu+1) * besselJ(seqRF, nu=model@par.model$nu)) / (seqRF^(model@par.model$nu)))[-1])
estY <- avgAcfDf / max(avgAcfDf)
estY <- estY[1:length(seqRF)]

plot(seqRF, trueY, ylim=c(-0.1, 1.1))
lines(seqRF, estY, col='red', lty=2, lwd=2)
estArea_0.25 <- areaBetween(seqRF, trueY, estY)
estArea_0.25
maxDist_0.25 <- maxDist(trueY, estY)
maxDist_0.25


load("from_my_pc/Cauchy_gamma_0.1_n10_p10_res_0.1_eps_0.3_cardinal_y.RData")
seqRF <- seq(0, 6, by=0.1)

trueY <- (1 + seqRF^2)^(-0.1) # c(1, ((2^(model@par.model$nu) * gamma(model@par.model$nu+1) * besselJ(seqRF, nu=model@par.model$nu)) / (seqRF^(model@par.model$nu)))[-1])
estY <- avgAcfDf / max(avgAcfDf)
estY <- estY[1:length(seqRF)]

plot(seqRF, trueY, ylim=c(-0.1, 1.1))
lines(seqRF, estY, col='red', lty=2, lwd=2)
estArea_0.3 <- areaBetween(seqRF, trueY, estY)
estArea_0.3
maxDist_0.3 <- maxDist(trueY, estY)
maxDist_0.3

load("from_my_pc/Cauchy_n10_p10_res_0.1_andriy_epsilon_2_kep_20.RData")
seqRF <- seq(0.1, 6, by=0.1)

trueY <- (1 + seqRF^2)^(-0.1) # c(1, ((2^(model@par.model$nu) * gamma(model@par.model$nu+1) * besselJ(seqRF, nu=model@par.model$nu)) / (seqRF^(model@par.model$nu)))[-1])
estY <- acf_fin
estY <- estY[1:length(seqRF)]

plot(seqRF, trueY, ylim=c(-0.1, 1.1))
lines(seqRF, estY, col='red', lty=2, lwd=2)
estArea_kep <- areaBetween(seqRF, trueY, estY)
estArea_kep
maxDist_kep <- maxDist(trueY, estY)
maxDist_kep

load("from_my_pc/Cauchy_gamma_0.1_n10_p10_res_0.1_Fourier.RData")
seqRF <- fsList[[1]]$x

# find closest to 6
fourierIndex <- which(round(seqRF, 1) == 6)
seqRF <- seqRF[1:fourierIndex]

trueY <- (1 + seqRF^2)^(-0.1) # c(1, ((2^(model@par.model$nu) * gamma(model@par.model$nu+1) * besselJ(seqRF, nu=model@par.model$nu)) / (seqRF^(model@par.model$nu)))[-1])
estY <- rho / max(rho)
estY <- estY[1:length(seqRF)]

plot(seqRF, trueY, ylim=c(-0.1, 1.1))
lines(seqRF, estY, col='red', lty=2, lwd=2)
estArea_fourier <- areaBetween(seqRF, trueY, estY)
estArea_fourier
maxDist_fourier <- maxDist(trueY, estY)
maxDist_fourier

load("from_my_pc/Cauchy_gamma_0.1_n10_p10_res_0.1_RFcov(bin=seq(0, 28, by=0.1)).RData")
seqRF <- seq(0, 6, by=0.1)
seqRF <- emp.cov@centers
empCovClosest <- which(plyr::round_any(seqRF, 0.1) == 6)[1]
seqRF <- seqRF[1:empCovClosest]

trueY <- (1 + seqRF^2)^(-0.1) # c(1, ((2^(model@par.model$nu) * gamma(model@par.model$nu+1) * besselJ(seqRF, nu=model@par.model$nu)) / (seqRF^(model@par.model$nu)))[-1])
# trueY <- c(1, ((2^(model@par.model$nu) * gamma(model@par.model$nu+1) * besselJ(emp.cov@centers, nu=model@par.model$nu)) / (emp.cov@centers^(model@par.model$nu)))[-1])
estY <- emp.cov@empirical
estY <- estY[1:length(seqRF)]

# Estimate mean and variance
mu.hat <- mean(zMatrix)
sample.var <- 1/(length(zMatrix) - 1) * sum((zMatrix - mu.hat)^2)

estY <- c(sample.var, estY[-1])

plot(seqRF, trueY, ylim=c(-0.1, 1.1))
lines(seqRF, estY, col='red', lty=2, lwd=2)
estArea_std <- areaBetween(seqRF, trueY, estY)
estArea_std
maxDist_std <- maxDist(trueY, estY)
maxDist_std

