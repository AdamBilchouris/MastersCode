library(readr)
library(tidyverse)
library(dplyr)
library(caret)
library(GA)
library(glue)
library(memoise)

MODEL_DIR <- 'models_txt_iqr'
RESULTS_DIR <- 'results_iqr'
IMAGES_DIR <- 'images_iqr'

iters <- 5

columns <-  c('full', 'aic', 'rmse', 'mae', 'p90', 'p95', 'p50', 'p95_90_50', 'percErr')
maeDf <- data.frame(matrix(nrow=0, ncol=length(columns)))
colnames(maeDf) <- columns

rmseDf <- data.frame(matrix(nrow=0, ncol=length(columns)))
colnames(rmseDf) <- columns

percErrDf <- data.frame(matrix(nrow=0, ncol=length(columns)))
colnames(percErrDf) <- columns


# Dmytro's code
df <- readRDS("Houses_75%_NAs_IQR.RDS")
df <- df[complete.cases(df), ]
df <- df %>%  select(-sewage)
df_90_complete_cases <- df 
##### full model
# We need the vector below the create  the quadratic effect of continuous variables for the model
nb_vars <- names(df_90_complete_cases[,])[apply(df_90_complete_cases,2,function(x) {!all(x %in% 0:1)})]
nb_vars <- nb_vars[nb_vars!="sold_price"]
nb_vars

n2 <- paste0("I(", nb_vars, "^2)") # dropping the dep. var - sold price and factor var.
n2

fml <- as.formula(paste("sold_price ~.+.^2+", paste(n2, collapse = "+"))) # .^2 - 2nd degree of interaction effect
# fml <- as.formula(paste("sold_price ~.+", paste(n2, collapse = "+"))) # .^2 - 2nd degree of interaction effect
fml


# FITNESS FUNCTIONS


for(it in 1:iters)
{
  train_rows <- sample(rownames(df_90_complete_cases), size=0.7*nrow(df_90_complete_cases))
  train <- df_90_complete_cases[rownames(df_90_complete_cases) %in% train_rows, ]
  test <- df_90_complete_cases[!rownames(df_90_complete_cases) %in% train_rows, ]
  
  model_1 <- lm(fml, data=train)
  summary(model_1)
  
  model_1_pred <- predict(model_1, newdata=test)
  MAE_model_1 <- mean(abs(test$sold_price - model_1_pred))
  RMSE_model_1 <- sqrt( mean((test$sold_price - model_1_pred)^2) )
  percErr_model_1 <- mean( abs((test$sold_price - model_1_pred) / test$sold_price) )
  
  maeDf[it, 'full'] <- MAE_model_1
  rmseDf[it, 'full'] <- RMSE_model_1
  percErrDf[it, 'full'] <- percErr_model_1
  
  
  # GA
  x <- model.matrix(model_1)[, -1]
  
  y <- train$sold_price
  
    fitness0 <- function(string) {
    inc <- which(string == 1)
    X <- cbind(1, x[, inc])
    mod <- .lm.fit(X, y)
    class(mod) <- "lm"
    
    -AIC(mod)
  }
  
  mFitness0 <- memoise(fitness0)
  
  parallelStart <- Sys.time()
  GA0 <- ga(
    "binary",
    fitness = mFitness0,
    nBits = ncol(x),
    names = colnames(x),
    monitor = plot,
    parallel = T
  )
  parallelEnd <- Sys.time()
  print(parallelEnd - parallelStart)
  
  c <- GA0
  indx <- which(c@solution[1, ]== 1)
  print('AIC DONE')
  
  f0 <- names(indx)
  f20Str <- paste('sold_price ~ ', paste(f0, collapse='+'))
  
  model_aic_best <- lm(f20Str, data=train)
  
  fileStr <- glue('{MODEL_DIR}/selectedGA_aic_{it}.txt')
  write(f20Str, file=fileStr)
  
  best_aic_pred <- predict(model_aic_best, newdata=test)
  MAE_aic <- mean(abs(test$sold_price - best_aic_pred))
  RMSE_aic <- sqrt( mean((test$sold_price - best_aic_pred)^2) )
  percErr_aic <- mean( abs((test$sold_price - best_aic_pred) / test$sold_price) )
  
  maeDf[it, 'aic'] <- MAE_aic
  rmseDf[it, 'aic'] <- RMSE_aic
  percErrDf[it, 'aic'] <- percErr_aic  
  
  
  fitness1 <- function(string) {
    inc <- which(string == 1)
    X <- cbind(1, x[, inc])
    mod <- .lm.fit(X, y)
    class(mod) <- "lm"
    
    mse <- mean((mod$residuals)^2)
    rmse <- sqrt(mse)
    -rmse
  }
  
  mFitness1 <- memoise(fitness1)
  
  parallelStart <- Sys.time()
  GA1 <- ga(
    "binary",
    fitness = mFitness1,
    nBits = ncol(x),
    names = colnames(x),
    monitor = plot,
    parallel = T
  )
  parallelEnd <- Sys.time()
  print(parallelEnd - parallelStart)
  
  c <- GA1
  indx <- which(c@solution[1, ]== 1)
  print('RMSE DONE')
  
  f1 <- names(indx)
  f2Str <- paste('sold_price ~ ', paste(f1, collapse='+'))
  
  model_rmse_best <- lm(f2Str, data=train)
  
  fileStr <- glue('{MODEL_DIR}/selectedGA_rmse_{it}.txt')
  write(f2Str, file=fileStr)
  
  best_rmse_pred <- predict(model_rmse_best, newdata=test)
  MAE_rmse <- mean(abs(test$sold_price - best_rmse_pred))
  RMSE_rmse <- sqrt( mean((test$sold_price - best_rmse_pred)^2) )
  percErr_rmse <- mean( abs((test$sold_price - best_rmse_pred) / test$sold_price) )
  
  maeDf[it, 'rmse'] <- MAE_rmse
  rmseDf[it, 'rmse'] <- RMSE_rmse
  percErrDf[it, 'rmse'] <- percErr_rmse
  
  #=======
  fitness2 <- function(string) {
    inc <- which(string == 1)
    # incc <<- inc
    X <- cbind(1, x[, inc])
    mod <- .lm.fit(X, y)
    class(mod) <- "lm"
    mae <- mean(abs(mod$residuals))
    -mae
  }
  
  mFitness2 <- memoise(fitness2)
  
  parallelStart <- Sys.time()
  GA2 <- ga(
    "binary",
    fitness = mFitness2,
    nBits = ncol(x),
    names = colnames(x),
    monitor = plot,
    parallel = T
  )
  parallelEnd <- Sys.time()
  print(parallelEnd - parallelStart)
  
  c2 <- GA2
  indx2 <- which(c2@solution[1, ]== 1)
  
  print('MAE DONE')
  
  f12 <- names(indx2)
  f22Str <- paste('sold_price ~ ', paste(f12, collapse='+'))
  
  model_mae_best <- lm(f22Str, data=train)
  
  fileStr <- glue('{MODEL_DIR}/selectedGA_mae_{it}.txt')
  write(f22Str, file=fileStr)
  
  best_mae_pred <- predict(model_mae_best, newdata=test)
  MAE_mae <- mean(abs(test$sold_price - best_mae_pred))
  RMSE_mae <- sqrt( mean((test$sold_price - best_mae_pred)^2) )
  percErr_mae <- mean( abs((test$sold_price - best_mae_pred) / test$sold_price) )
  
  maeDf[it, 'mae'] <- MAE_mae
  rmseDf[it, 'mae'] <- RMSE_mae
  percErrDf[it, 'mae'] <- percErr_mae
  
  #=======
  fitness3 <- function(string) {
    inc <- which(string == 1)
    X <- cbind(1, x[, inc])
    mod <- lm.fit(X, y)
    class(mod) <- "lm"
    
    # Extract variables and make it into a new model
    lmStr <- paste(names(mod$coefficients), collapse='+')
    lmStr2 <- substr(lmStr, 2, nchar(lmStr))
    lmStr3 <- paste0('sold_price ~ ', lmStr2)
    lmMod <- lm(lmStr3, data=train)
    
    isIn90Percent <- 0
    tryCatch(
      expr={
        # Compute prediciton interval for this model
        pred90 <- predict(lmMod, newdata=train, interval='prediction', level=0.90)
        
        for(i in 1:nrow(train)) {
          # Check if the true value is within the 90% prediction interval generated from the model
          if(train[i, 'sold_price'] >= pred90[i, 'lwr'] && train[i, 'sold_price'] <= pred90[i, 'upr']) {
            isIn90Percent <- isIn90Percent + 1
          }
        }
        
        return(isIn90Percent)
      },
      error = function(e) {
        message(e)
        # Return 0 as it contributes nothing to the overall fitness of the population
        return(0)
      }
    )
  }
  
  mFitness3 <- memoise(fitness3)
  
  parallelStart <- Sys.time()
  GA3 <- ga(
    "binary",
    fitness = mFitness3,
    nBits = ncol(x),
    names = colnames(x),
    monitor = plot,
    parallel=T
  )
  parallelEnd <- Sys.time()
  print(parallelEnd - parallelStart)
  
  c3 <- GA3
  indx3 <- which(c3@solution[1, ]== 1)
  
  print('90% DONE')
  
  f13 <- names(indx3)
  f23Str <- paste('sold_price ~ ', paste(f13, collapse='+'))
  
  model_p90_best <- lm(f23Str, data=train)
  
  fileStr <- glue('{MODEL_DIR}/selectedGA_p90_{it}.txt')
  write(f23Str, file=fileStr)
  
  best_p90_pred <- predict(model_p90_best, newdata=test)
  MAE_p90 <- mean(abs(test$sold_price - best_p90_pred))
  RMSE_p90 <- sqrt( mean((test$sold_price - best_p90_pred)^2) )
  percErr_p90 <- mean( abs((test$sold_price - best_p90_pred) / test$sold_price) )
  
  maeDf[it, 'p90'] <- MAE_p90
  rmseDf[it, 'p90'] <- RMSE_p90
  percErrDf[it, 'p90'] <- percErr_p90
  
  #=======
  fitness4 <- function(string) {
    inc <- which(string == 1)
    X <- cbind(1, x[, inc])
    mod <- lm.fit(X, y)
    class(mod) <- "lm"
    
    # Extract variables and make it into a new model
    lmStr <- paste(names(mod$coefficients), collapse='+')
    lmStr2 <- substr(lmStr, 2, nchar(lmStr))
    lmStr3 <- paste0('sold_price ~ ', lmStr2)
    lmMod <- lm(lmStr3, data=train)
    
    isIn95Percent <- 0
    tryCatch(
      expr={
        # Compute prediciton interval for this model
        pred95 <- predict(lmMod, newdata=train, interval='prediction', level=0.95)
        
        for(i in 1:nrow(train)) {
          # Check if the true value is within the 95% prediction interval generated from the model
          if(train[i, 'sold_price'] >= pred95[i, 'lwr'] && train[i, 'sold_price'] <= pred95[i, 'upr']) {
            isIn95Percent <- isIn95Percent + 1
          }
        }
        
        return(isIn95Percent)
      },
      error = function(e) {
        message(e)
        # Return 0 as it contributes nothing to the overall fitness of the population
        return(0)
      }
    )
  }
  
  mFitness4 <- memoise(fitness4)
  
  parallelStart <- Sys.time()
  GA4 <- ga(
    "binary",
    fitness = mFitness4,
    nBits = ncol(x),
    names = colnames(x),
    monitor = plot,
    parallel = T
  )
  parallelEnd <- Sys.time()
  print(parallelEnd - parallelStart)
  
  c4 <- GA4
  indx4 <- which(c4@solution[1, ]== 1)
  
  print('95% DONE')
  
  f14 <- names(indx4)
  f24Str <- paste('sold_price ~ ', paste(f14, collapse='+'))
  
  model_p95_best <- lm(f24Str, data=train)
  
  fileStr <- glue('{MODEL_DIR}/selectedGA_p95_{it}.txt')
  write(f24Str, file=fileStr)
  
  best_p95_pred <- predict(model_p95_best, newdata=test)
  MAE_p95 <- mean(abs(test$sold_price - best_p95_pred))
  RMSE_p95 <- sqrt( mean((test$sold_price - best_p95_pred)^2) )
  percErr_p95 <- mean( abs((test$sold_price - best_p95_pred) / test$sold_price) )
  
  maeDf[it, 'p95'] <- MAE_p95
  rmseDf[it, 'p95'] <- RMSE_p95
  percErrDf[it, 'p95'] <- percErr_p95
  
  #=======
  fitness5 <- function(string) {
    inc <- which(string == 1)
    X <- cbind(1, x[, inc])
    mod <- lm.fit(X, y)
    class(mod) <- "lm"
    
    # Extract variables and make it into a new model
    lmStr <- paste(names(mod$coefficients), collapse='+')
    lmStr2 <- substr(lmStr, 2, nchar(lmStr))
    lmStr3 <- paste0('sold_price ~ ', lmStr2)
    lmMod <- lm(lmStr3, data=train)
    
    isIn50Percent <- 0
    tryCatch(
      expr={
        # Compute prediciton interval for this model
        pred50 <- predict(lmMod, newdata=train, interval='prediction', level=0.50)
        
        for(i in 1:nrow(train)) {
          # Check if the true value is within the 95% prediction interval generated from the model
          if(train[i, 'sold_price'] >= pred50[i, 'lwr'] && train[i, 'sold_price'] <= pred50[i, 'upr']) {
            isIn50Percent <- isIn50Percent + 1
          }
        }
        
        return(isIn50Percent)
      },
      error = function(e) {
        message(e)
        # Return 0 as it contributes nothing to the overall fitness of the population
        return(0)
      }
    )
  }
  
  mFitness5 <- memoise(fitness5)
  
  parallelStart <- Sys.time()
  GA5 <- ga(
    "binary",
    fitness = fitness5,
    nBits = ncol(x),
    names = colnames(x),
    monitor = plot,
    parallel = T
  )
  parallelEnd <- Sys.time()
  print(parallelEnd - parallelStart)
  
  c5 <- GA5
  indx5 <- which(c5@solution[1, ]== 1)
  print('50% DONE')
  
  f15 <- names(indx5)
  f25Str <- paste('sold_price ~ ', paste(f15, collapse='+'))
  
  model_p50_best <- lm(f25Str, data=train)
  summary(model_p50_best)
  
  fileStr <- glue('{MODEL_DIR}/selectedGA_p50_{it}.txt')
  write(f25Str, file=fileStr)
  
  best_p50_pred <- predict(model_p50_best, newdata=test)
  MAE_p50 <- mean(abs(test$sold_price - best_p50_pred))
  RMSE_p50 <- sqrt( mean((test$sold_price - best_p50_pred)^2) )
  percErr_p50 <- mean( abs((test$sold_price - best_p50_pred) / test$sold_price) )
  
  maeDf[it, 'p50'] <- MAE_p50
  rmseDf[it, 'p50'] <- RMSE_p50
  percErrDf[it, 'p50'] <- percErr_p50
  
  #=======
  sum95_90_50 <- 50 + 90 + 95
  wght <- c(95/sum95_90_50, 90/sum95_90_50, 50/sum95_90_50)
  
  fitness6 <- function(string) {
    inc <- which(string == 1)
    X <- cbind(1, x[, inc])
    mod <- lm.fit(X, y)
    class(mod) <- "lm"
    
    # Extract variables and make it into a new model
    lmStr <- paste(names(mod$coefficients), collapse='+')
    lmStr2 <- substr(lmStr, 2, nchar(lmStr))
    lmStr3 <- paste0('sold_price ~ ', lmStr2)
    lmMod <- lm(lmStr3, data=train)
    
    
    isIn50Percent <- 0
    isIn90Percent <- 0
    isIn95Percent <- 0
    
    tryCatch(
      expr={
        # Compute prediciton interval for this model
        pred95 <- predict(lmMod, newdata=train, interval='prediction', level=0.95)
        pred90 <- predict(lmMod, newdata=train, interval='prediction', level=0.90)
        pred50 <- predict(lmMod, newdata=train, interval='prediction', level=0.50)
        
        for(i in 1:nrow(train)) {
          if(train[i, 'sold_price'] >= pred95[i, 'lwr'] && train[i, 'sold_price'] <= pred95[i, 'upr']) {
            isIn95Percent <- isIn95Percent + 1
          }
          if(train[i, 'sold_price'] >= pred90[i, 'lwr'] && train[i, 'sold_price'] <= pred90[i, 'upr']) {
            isIn90Percent <- isIn90Percent + 1
          }
          if(train[i, 'sold_price'] >= pred50[i, 'lwr'] && train[i, 'sold_price'] <= pred50[i, 'upr']) {
            isIn50Percent <- isIn50Percent + 1
          }
        }
        
        # Weighted sum of these
        wmean <- (wght[1] * isIn95Percent) + (wght[2] * isIn90Percent) + (wght[3] * isIn50Percent)
        return(wmean)
      },
      error = function(e) {
        message(e)
        # Return 0 as it contributes nothing to the overall fitness of the population
        return(0)
      }
    )
    
  }
  
  mFitness6 <- memoise(fitness6)
  
  GA6 <- ga(
    "binary",
    fitness = mFitness6,
    nBits = ncol(x),
    names = colnames(x),
    monitor = plot,
    parallel = T
  )
  parallelEnd <- Sys.time()
  print(parallelEnd - parallelStart)
  
  c6 <- GA6
  indx6 <- which(c6@solution[1, ]== 1)
  
  print('95_90_50% DONE')
  
  f16 <- names(indx6)
  f26Str <- paste('sold_price ~ ', paste(f16, collapse='+'))
  
  model_p95_90_50_best <- lm(f26Str, data=train)
  summary(model_p95_90_50_best)
  
  fileStr <- glue('{MODEL_DIR}/selectedGA_p95_90_50__{it}.txt')
  write(f26Str, file=fileStr)
  
  best_p95_90_50_pred <- predict(model_p95_90_50_best, newdata=test)
  MAE_p95_90_50 <- mean(abs(test$sold_price - best_p95_90_50_pred))
  RMSE_p95_90_50 <- sqrt( mean((test$sold_price - best_p95_90_50_pred)^2) )
  percErr_p95_90_50 <- mean( abs((test$sold_price - best_p95_90_50_pred) / test$sold_price) )
  
  maeDf[it, 'p95_90_50'] <- MAE_p95_90_50
  rmseDf[it, 'p95_90_50'] <- RMSE_p95_90_50
  percErrDf[it, 'p95_90_50'] <- percErr_p95_90_50
  
  
  #=======
  fitness7 <- function(string) {
    inc <- which(string == 1)
    X <- cbind(1, x[, inc])
    mod <- .lm.fit(X, y)
    class(mod) <- "lm"
    
    percErr <- mean(abs(mod$residuals / train$sold_price))
    -percErr
  }
  
  mFitness7 <- memoise(fitness7)
  
  parallelStart <- Sys.time()
  GA7 <- ga(
    "binary",
    fitness = mFitness7,
    nBits = ncol(x),
    names = colnames(x),
    monitor = plot,
    parallel = T
  )
  
  c7 <- GA7
  indx7 <- which(c7@solution[1, ]== 1)
  
  print('percErr DONE')
  
  f17 <- names(indx7)
  f27Str <- paste('sold_price ~ ', paste(f17, collapse='+'))
  
  model_percErr_best <- lm(f27Str, data=train)
  summary(model_percErr_best)
  
  fileStr <- glue('{MODEL_DIR}/selectedGA_percErr_{it}.txt')
  write(f27Str, file=fileStr)
  
  best_percErr_pred <- predict(model_percErr_best, newdata=test)
  MAE_percErr <- mean(abs(test$sold_price - best_percErr_pred))
  RMSE_percErr <- sqrt( mean((test$sold_price - best_percErr_pred)^2) )
  percErr_percErr <- mean( abs((test$sold_price - best_percErr_pred) / test$sold_price) )
  
  maeDf[it, 'percErr'] <- MAE_percErr
  rmseDf[it, 'percErr'] <- RMSE_percErr
  percErrDf[it, 'percErr'] <- percErr_percErr
  
  maeFile <- glue('{RESULTS_DIR}/maeDf_{it}.csv')
  rmseFile <- glue('{RESULTS_DIR}/rmseDf_{it}.csv')
  percErrFile <- glue('{RESULTS_DIR}/percErrDf_{it}.csv')
  
  write.csv(maeDf, file=maeFile, row.names=F)
  write.csv(rmseDf, file=rmseFile, row.names=F)
  write.csv(percErrDf, file=percErrFile, row.names=F)
  
  imageFile <- glue('{IMAGES_DIR}/image_{it}.RData')
  save.image(file=imageFile)
}

maeFile <- glue('{RESULTS_DIR}/maeDf_all.csv')
rmseFile <- glue('{RESULTS_DIR}/rmseDf_all.csv')
percErrFile <- glue('{RESULTS_DIR}/percErrDf_all.csv')

write.csv(maeDf, file=maeFile, row.names=F)
write.csv(rmseDf, file=rmseFile, row.names=F)
write.csv(percErrDf, file=percErrFile, row.names=F)
