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

thetas <- seq(0, pi, length.out=63)
epsilon <- 2
kep <- floor(epsilon/xRes)
acf_len <- rep(0,2*kep+1)
start_vec <- c(0,0)
end_vec <- c(1,1)
set_of_acf_round <- c()
set_of_acf_points <- c()

startT <- Sys.time()
for (t in 1:length(thetas)) {
  a <- getPointsInTube(start_vec, end_vec, angle=thetas[t], epsilon, x, y,  e.type='y', parallel.line=F, y.shift=NA, just.line=T)
  set_of_acf_dir <- c()
  for (i in (-kep):kep) {
    a$baseLine1 <- a$baseLine

    ix <- 0
    iy <- i

    a$baseLine1$x <- a$baseLine1$x + xRes*ix
    a$baseLine1$y <- a$baseLine1$y + yRes*iy
    a$baseLine1 <- a$baseLine1[a$baseLine1$x >=-10 & a$baseLine1$x <= 10 & a$baseLine1$y >=-10 & a$baseLine1$y <= 10,]

    # lines(a$baseLine1, col=i+21)

    acf_dir <-  acf(zMatrix[cbind( (a$baseLine1$x-x[1]) / xRes + 1, (a$baseLine1$y-y[1]) / yRes + 1)], type='correlation', plot=F, lag.max=length(a$line_p1$x))
    acf_len[i + kep + 1] <- max(acf_dir$lag) + 1
    set_of_acf_dir <- c(set_of_acf_dir, list(acf_dir$acf))
  }

  shortest_length <- min(acf_len)

  # Create a matrix with each row being one vector, cut to the shortest length
  my_matrix <- matrix(nrow = length(set_of_acf_dir), ncol = shortest_length)
  for (j in seq_along(set_of_acf_dir)) {
    my_matrix[j, ] <- head(set_of_acf_dir[[j]], shortest_length)
  }

  my_df <- as.data.frame(my_matrix)
  mean_covar <- colMeans(my_df)

  set_of_acf_round <- c(set_of_acf_round, list(mean_covar))
  dirRes <- (1: shortest_length)* sqrt((a$baseLine$x[1]-a$baseLine$x[2])^2 + (a$baseLine$y[1]-a$baseLine$y[2])^2)

  dirRes <- round(dirRes * 10) / 10
  set_of_acf_points <-  c(set_of_acf_points, list(dirRes))
}

unique_distances <- sort(unique(unlist(set_of_acf_points)))
acf_fin <- rep(0,length(unique_distances))
k <- 1
for (num in unique_distances){
  indices_lst <- sapply(set_of_acf_points, function(x) which(x == num))
  indices_vec <- unlist(indices_lst)
  indices_list <- rep(seq_along(set_of_acf_points), sapply(indices_lst, length))

  sum1 <- 0
  for (i in 1: length(indices_list)){sum1 <- sum1 + set_of_acf_round [[indices_list[i]]][indices_vec[i]]}
  acf_fin[k] <- sum1/length(indices_list)
  k <- k+1
}

endT <- Sys.time()
print(endT - startT)

plot(unique_distances, acf_fin)

