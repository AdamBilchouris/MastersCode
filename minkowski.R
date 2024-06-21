# Running a number of iterations of the overall tests to check for normality
library(RandomFields)
library(plot3D)
library(doParallel) 
library(glue)
library(nortest)
library(moments)

registerDoParallel(detectCores())

x <- seq(-10, 10, 0.1)
xLen <- length(x)
y <- x
yLen <- xLen

getMinkowski <- function(m, a, k) {
    whichGreaterA <- ifelse((m^k) > a, 1, 0)
    return(sum(whichGreaterA == 1))
}

getMinkowskiRF <- function(model, x, y, a, k) {
    sim <- RFsimulate(model, x=x,y=y)
    zMatrix <- RFspDataFrame2conventional(sim)$data
    return(getMinkowski(zMatrix, a, k))
}

runSims <- function(iters, model, x, y, a, k) {
    minkEst <- foreach(i=1:iters, .combine=c, .packages=c('RandomFields'),
    .export=c('getMinkowski', 'getMinkowskiRF')) %dopar% {
        getMinkowskiRF(model, x, y, a, k)
    }
    llTest <- lillie.test(minkEst)
    
    return( rbind('lilliefors'=llTest$p.value) )
}

iters <- 4
resDf <- data.frame(matrix(nrow=iters, ncol=1))
colnames(resDf) <- c('lilliefors')

model <- RMgauss(scale=1)
# model <- RMbessel(nu=1, scale=1)
# model <- RMcauchy(gamma=1, scale=1)

t1 <- Sys.time()

for(i in 1:iters) {
    resDf[i, ] <- runSims(5000, model, x, y, 0.01, 1)
}

View(resDf)
t2 <- Sys.time()
print(t2 - t1)

resDf2 <- data.frame(matrix(nrow=iters, ncol=1))
colnames(resDf2) <- c('lilliefors')
t1 <- Sys.time()
for(i in 1:iters) {
    resDf2[i, ] <- runSims(5000, model, x, y, 0.01, 2)
}
t2 <- Sys.time()
print(t2 - t1)
View(resDf2)

resDf3 <- data.frame(matrix(nrow=iters, ncol=1))
colnames(resDf3) <- c('lilliefors')
t1 <- Sys.time()
for(i in 1:iters) {
    resDf3[i, ] <- runSims(5000, model, x, y, 0.01, 3)
}
t2 <- Sys.time()
print(t2 - t1)
View(resDf3)

resDf4 <- data.frame(matrix(nrow=iters, ncol=1))
colnames(resDf4) <- c('lilliefors')
t1 <- Sys.time()
for(i in 1:iters) {
    resDf4[i, ] <- runSims(5000, model, x, y, 0.01, 4)
}
t2 <- Sys.time()
print(t2 - t1)
View(resDf4)

hist(minkEst, freq=F)
curve(dnorm(x, mean=mean(minkEst), sd=sd(minkEst)), lwd=2, lty=2, col='red', add=T)