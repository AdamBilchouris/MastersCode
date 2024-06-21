recursiveGegen <- function(n, d, u) {
  # Index 1 is actually C0
  n <- n+1
  if(abs(u) > 1) {
    print('|u| must be <= 1')
    return(-1)
  }
  
  if(n < 3) {
    print('n must be >= 3')
    return(-1)
  }
  # c_0, c_1, (c_2, ..., c_(n-1))
  coefs <- c(1, 2*d*u, rep(0, n-2))
  
  # 3='2'
  for(i in 3:n) {
    coefs[i] <- (1/(i-1)) * ( (2*u*((i-1) + d - 1)*coefs[i-1]) - ( ((i-1) + 2*d - 2)*coefs[i-2] ) )
  }
  
  return(coefs)
}

# k=1 only
gegenRealisation <- function(size, n, d, u, mu=0, sigma2=1) {
  # Formula (7) and (8) of Forecasting with k-factor Gegenbauer Processes: Theory and Applications, L. FERRARA AND D. GUEGAN (2001) pages 4-5.
  if(length(d) != 1) {
    print('d must be a number')
    return(-1)
  } 
  if(length(u) != 1) {
    print('u must be a number')
    return(-1)
  }
  if(abs(u) > 1) {
    print('|u| must be <= 1')
    return(-1)
  }
  psi <- recursiveGegen(n, d, u)
  
  epsilon <- rnorm(size + n, mu, sigma2)
  
  realisation <- c()
  
  # For the realisation t, sum over all j>=0, so all psi, multiplied by B^{j} epsilon t
  # So, we have the epsilons: epsilon t, t-1, t-2, ..., t-j.
  #
  # Always going to be a difference of n (the number of Gegenbauer coefficients).
  # In Formula (7), this is 'j'. j \in \{ 0, 1, ..., n \}
  for(j in 1:size) {
    realisation[j] <- sum(psi * epsilon[(n + j):j])
  }
  return(realisation)
}

# Stationarity
# |u| M<
u <- 0.3
d <- 0.2
n <- 200
nSize <- 10000

library(glue)

gegenReal <- gegenRealisation(nSize, n, d, u)

plot(gegenReal, type='l')
acf(gegenReal, type='correlation', lag.max=100)
specGegen <- spectrum(gegenReal, log='no', method='pgram')
specGegenAR <- spectrum(gegenReal, log='no', method='ar')

# True Gegenbauer frequency
trueGF <- acos(u) / (2*pi)
print(paste('True Gegenbauer Frequency: ', trueGF))
# Natural estimator of GF
idx_max <- which(specGegen$spec == max(specGegen$spec))
estGF <- specGegen$freq[idx_max]
print(paste('Esimated Gegenbauer Frequency: ', estGF))
diffGF <- trueGF - estGF
print(paste('Gegenbauer frequency difference: ', diffGF))
maxSpec <- specGegen$spec[idx_max]

# GF = arccos(u)
# Fourier transform is scaled by 2pi, so undo scaling.
estU <- cos(2*pi*estGF)
glue('Estimated u: {estU}')
diffU <- u - estU
glue('Difference u: {diffU}')


# b / |x-w0|^(-2*a)
f_1 <- function(param, gegenFreq, gegenSpec, GF) {
  denom <- abs(gegenFreq - GF)
  for(i in 1:length(denom)) {
    if(is.na(denom[i])) {
      denom[i] <- 1e10
    }
    if(is.infinite(denom[i])) {
      denom[i] <- 1e10
    }
    if(denom[i] == 0) {
      denom[i] <- 1e10
    }
  }
  f.x <- param[2] * (denom^(-2*param[1]))
  
  sum(abs(gegenSpec - f.x))
}

# f(w0)*b / |x-w0|^(a)
f_2 <- function(param, gegenFreq, gegenSpec, GF, maxSpec) {
  denom <- abs(gegenFreq - GF)
  for(i in 1:length(denom)) {
    if(is.na(denom[i])) {
      denom[i] <- 1e10
    }
    if(is.infinite(denom[i])) {
      denom[i] <- 1e10
    }
    if(denom[i] == 0) {
      denom[i] <- 1e10
    }
  }
  f.x <- (maxSpec*param[2]) / (denom^(2 * param[1]))
  
  sum(abs(gegenSpec - f.x))
}

# ( (f(w0)*b / |x-w0|) + c )^(2*a) <- not a very good function.
# [ f(w0)*b / |x-w0|^(2*a) ] + c <- Seems better for certain initial values of c.
f_3 <- function(param, gegenFreq, gegenSpec, GF, maxSpec) {
  denom <- abs(gegenFreq - GF)
  for(i in 1:length(denom)) {
    if(is.na(denom[i])) {
      denom[i] <- 1e10
    }
    if(is.infinite(denom[i])) {
      denom[i] <- 1e10
    }
    if(denom[i] == 0) {
      denom[i] <- 1e10
    }
  }
  
  f.x <- ( (maxSpec*param[2]) * (denom)^(-2*param[1]) ) + param[3]
  sum(abs(gegenSpec - f.x))
}

# b / |x-c|^(2*a)
f_unknown_1 <- function(param, gegenFreq, gegenSpec) {
  denom <- abs(gegenFreq - param[3])
  for(i in 1:length(denom)) {
    if(is.na(denom[i])) {
      denom[i] <- 1e10
    }
    if(is.infinite(denom[i])) {
      denom[i] <- 1e10
    }
    if(denom[i] == 0) {
      denom[i] <- 1e10
    }
  }
  
  f.x <- param[2] / ( (denom)^(2*param[1]) )
  sum(abs(gegenSpec - f.x))
}

# From Ferrara paper
ferraraSpectral <- function(param, gegenFreq, gegenSpec, GF, sigmaWN) {
  f.x <- (sigmaWN/(2*pi)) * ( abs( 4 * sin( (gegenFreq + GF) / 2 ) * sin( (gegenFreq - GF) / 2) ) )^(-2*param)
  for(i in 1:length(f.x)) {
    if(is.na(f.x[i])) {
      f.x[i] <- 1e10
    }
    if(is.infinite(f.x[i])) {
      f.x[i] <- 1e10
    }
    if(f.x[i] == 0) {
      f.x[i] <- 1e10
    }
  }
  sum(abs(gegenSpec - f.x))
}

res_f_1 <- optim(c(0.5, 0.01), fn=f_1, gegenFreq=specGegen$freq, gegenSpec=specGegen$spec, GF=estGF)
res_f_1$par
plot(specGegen$freq, specGegen$spec, type='l')
curve((res_f_1$par[2]) / ( abs(x - estGF)^(2*res_f_1$par[1]) ), col='green', lwd='2', add=T)
curve((res_f_1$par[2]) / ( abs(x - estGF)^(2*res_f_1$par[1]) ), col='green', lwd='2')

estD <- res_f_1$par[1]
glue('Estimated d: {estD}')
diffD <- d - estD
glue('Difference in d: {diffD}')

res_f_2 <- optim(c(0.5, 0.01), fn=f_2, gegenFreq=specGegen$freq, gegenSpec=specGegen$spec, GF=estGF, maxSpec=maxSpec)
res_f_2$par
plot(specGegen$freq, specGegen$spec, type='l')
curve((maxSpec*res_f_2$par[2]) / ( abs(x - estGF)^(2*res_f_2$par[1]) ), col='green', lwd='2', add=T)
curve((maxSpec*res_f_2$par[2]) / ( abs(x - estGF)^(2*res_f_2$par[1]) ), col='green', lwd='2')

estD <- res_f_2$par[1]
glue('Estimated d: {estD}')
diffD <- d - estD
glue('Difference in d: {diffD}')

res_f_3 <- optim(c(0.5, 0.01, 0.1), fn=f_3, gegenFreq=specGegen$freq, gegenSpec=specGegen$spec, GF=estGF, maxSpec=maxSpec)
res_f_3$par

plot(specGegen$freq, specGegen$spec, type='l')
curve( ( (maxSpec*res_f_3$par[2]) * abs(x - estGF)^(-2*res_f_3$par[1]) ) + res_f_3$par[3],
       col='green', lwd='2', add=T)

curve( ( (maxSpec*res_f_3$par[2]) * abs(x - estGF)^(-2*res_f_3$par[1]) ) + res_f_3$par[3], col='green', lwd='2')

estD <- res_f_3$par[1]
glue('Estimated d: {estD}')
diffD <- d - estD
glue('Difference in d: {diffD}')

res_unknown_1 <- optim(c(0.5, 0.01, 0.1), fn=f_unknown_1, gegenFreq=specGegen$freq, gegenSpec=specGegen$spec)
res_unknown_1$par

plot(specGegen$freq, specGegen$spec, type='l')
curve( res_unknown_1$par[2] / ( abs(x - res_unknown_1$par[3])^(2*res_unknown_1$par[1]) ),
       col='green', lwd='2', add=T)
curve( res_unknown_1$par[2] / ( abs(x - res_unknown_1$par[3])^(2*res_unknown_1$par[1]) ), col='green', lwd='2')

estD <- res_unknown_1$par[1]
glue('Estimated d: {estD}')
diffD <- d - estD
glue('Difference in d: {diffD}')

estGF <- res_unknown_1$par[3]
glue('Estimated GF: {estGF}')
diffGF <- trueGF - estGF
glue('Difference in GF: {diffGF}')

res_f_ferrara <- optim(0.5, fn=ferraraSpectral, gegenFreq=specGegen$freq, gegenSpec=specGegen$spec, GF=estGF, sigmaWN=1)
res_f_ferrara$par

plot(specGegen$freq, specGegen$spec, type='l')
curve((1/(2*pi)) * abs( 4*sin((x + trueGF)/2)*sin((x - trueGF)/2))^(-2*res_f_ferrara$par), col='red', lwd=2, add=T)

# estD <- res_f_ferrara$minimum
estD <- res_f_ferrara$par
glue('Estimated d: {estD}')
diffD <- d - estD
glue('Difference in d: {diffD}')

