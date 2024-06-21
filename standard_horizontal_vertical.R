library(RandomFields)
library(plot3D)
library(doParallel)
library(glue)
library(dplyr)

registerDoParallel(detectCores())

# Generate a random field, or read one from memory
x <- seq(-30, 30, 0.1)
xLen <- length(x)
y <- x

xRes <- abs(x[1] - x[2])
yRes <- abs(y[1] - y[2])
cbind(xRes, yRes)

model <- RMgauss(var=1, scale=1)
# model <- RMbessel(var=1, nu=1, scale=1)
sim <- RFsimulate(model, x=x, y=y)
zMatrix <- RFspDataFrame2conventional(sim)$data
image2D(x=x, y=y, z=zMatrix)

# Estimate mean and variance
mu.hat <- mean(zMatrix)
mu.hat
sample.var <- 1/(length(zMatrix) - 1) * sum((zMatrix - mu.hat)^2)
sample.var

emp.cov <- RFcov(data=sim, bin=seq(0, 28, by=0.1))
plot(emp.cov, model=model)

# Using lines and then time-series methods to derive the isotropic correlation function.
verticalLines <- list()
horizontalLines <- list()

# Vertical
for(j in 1:ncol(zMatrix)) {
    column <- zMatrix[, j]
    verticalLines[[j]] <- column
}

# Horizontal
for(i in 1:nrow(zMatrix)) {
    row <- zMatrix[i, ]
    horizontalLines[[i]] <- row
}

# Average verticals without FT
verticalAcf <- list()
horizontalAcf <- list()
for(i in 1:length(x)) {
    verticalAcf[[i]] <- acf(verticalLines[[i]], type='correlation', plot=F, lag.max=length(verticalLines[[i]]))$acf
    horizontalAcf[[i]] <- acf(horizontalLines[[i]], type='correlation', plot=F, lag.max=length(horizontalLines[[i]]))$acf
}

verticalAcfMatrix <- matrix(unlist(verticalAcf), nrow=length(verticalAcf), ncol=length(verticalAcf[[1]]), byrow=T)
horizontalAcfMatrix <- matrix(unlist(horizontalAcf), nrow=length(horizontalAcf), ncol=length(horizontalAcf[[1]]), byrow=T)

estimatedAcf <- colMeans(rbind(colMeans(verticalAcfMatrix), colMeans(horizontalAcfMatrix)))