Singleblock_algo_Param_init <- function(y, x, u, K, tau){
  # proportions
  pi <- runif(K)
  pi <- pi / sum(pi)
  # parameters of clustering
  alpha <- list()
  see <- sample(1:nrow(x), K)
  for (j in 1:ncol(x)){
    if (is.factor(x[,j])){
      m <- nlevels(x[,j])
      alpha[[j]] <- matrix(runif(m * K), K, m, dimnames = list(paste0("Class.", 1:K), paste0("Level.", 1:m)))
      for (k in 1:K){
        if (! is.na(x[see[k],j])){
          alpha[[j]][k,x[see[k],j]] <- alpha[[j]][k,x[see[k],j]] + 1
          alpha[[j]][k,] <- alpha[[j]][k,] / sum(alpha[[j]][k,])
        }
      }
    }else{
      alpha[[j]] <- matrix(NA, K, 2, dimnames = list(paste0("Class.", 1:K), c("mean", "variance")))

      for (k in 1:K){
        if (!is.na(x[see[k],j])){
          alpha[[j]][k,1] <- x[see[[k]],j]
        }else{
          alpha[[j]][k,1] <- mean(x[,j], na.rm = TRUE)
        }
      }
      alpha[[j]][,2] <- var(x[,j], na.rm = TRUE)
    }
  }
  # parameters of prediction
  dv <- ifelse(!is.null(u), ncol(u), 0) + K
  beta <- runif(dv, min = -1, max = 1)
  list(pi = pi,
       alpha = alpha,
       beta = beta,
       tunereg = 1,
       tau = tau)
}

Singleblock_algo_Param_PDF <- function(y, x, u, model.reg, param){
  K <- length(param$pi)
  out <- matrix(log(param$pi), length(y), K, byrow = TRUE)
  for (k in 1:K){
    # parameters of clustering
    for (j in 1:ncol(x)){
      obs <- which(!is.na(x[,j]))
      if (is.factor(x[,j])){
        out[obs,k] <- out[obs,k] + log(param$alpha[[j]][k,x[obs,j]])
      }else{
        out[obs,k] <- out[obs,k] + dnorm(x[obs,j], param$alpha[[j]][k,1], sqrt(param$alpha[[j]][k,2]), log = TRUE)
      }
    }
    # parameters of regression
    if (model.reg == "mean"){
      if (!is.null(u)){
        out[,k] <- out[,k] + dnorm(y, u%*%param$beta[1:ncol(u)] + param$beta[k + ncol(u)], sd = sqrt(param$tunereg), log = TRUE)
      }else{
        out[,k] <- out[,k] + dnorm(y, param$beta[k], sd =  sqrt(param$tunereg), log = TRUE)
      }
    }else if (model.reg == "quantile"){
      if (!is.null(u)){
        out[,k] <- out[,k] + log(dALD(y = y, mu = u%*%param$beta[1:ncol(u)] + param$beta[k + ncol(u)], sigma =  param$tunereg, p = param$tau))
      }else{
        out[,k] <- out[,k] + log(dALD(y = y, mu = param$beta[k], sigma = param$tunereg, p = param$tau))
      }
    }
  }
  out
}

Singleblock_algo_Param_Mstep <- function(y, x, u, model.reg, tik, param){
  # proportions
  K <- ncol(tik)
  n <- length(y)
  param$pi <- colSums(tik) / n
  # parameters of clustering
  for (k in 1:K){
    # parameters of clustering
    for (j in 1:ncol(x)){
      if (is.factor(x[,j])){
        for (h in 1:nlevels(x[,j])){
          param$alpha[[j]][k,h] <- sum(tik[which(x[,j]==h), k])
        }
        param$alpha[[j]][k,] <- param$alpha[[j]][k,] / sum(param$alpha[[j]][k,])
      }else{
        param$alpha[[j]][k, 1] <- weighted.mean(x = x[,j], w = tik[,k], na.rm = TRUE)
        param$alpha[[j]][k, 2] <- weighted.mean(x = (x[,j]-param$alpha[[j]][k,1])**2, w = tik[,k], na.rm = TRUE)
      }
    }
  }
  # parameters of regression
  covariates <- u
  if (K>1){
    for (k in 2:K) covariates <- rbind(covariates, u)
  }
  covariates <- cbind.data.frame(covariates, factor(rep(1:K, each = n)))
  covariates <- cbind.data.frame(rep(y, K), covariates)
  colnames(covariates)[1] <- "Y"

  if (model.reg == "mean"){
    if (K>1){
      tmp <- glm(Y~. + 0, weights = as.numeric(tik), data = covariates)
      param$tunereg <- sum((tmp$residuals**2) * as.numeric(tik)) /n
      param$beta <- coefficients(tmp)
    }else{
      tmp <- glm(Y~., weights = as.numeric(tik), data = covariates[,-ncol(covariates)])
      param$tunereg <- sum((tmp$residuals**2) * as.numeric(tik)) /n
      param$beta <- c(coefficients(tmp)[-1], coefficients(tmp)[1])
    }
  }else if (model.reg == "quantile"){
    if (K>1){
      tmp <- rq(Y~. + 0, weights = as.numeric(tik), data = covariates, tau = param$tau)
      resid <- tmp$residuals
      param$tunereg <- - sum( (resid * ((resid<0) - param$tau))  * as.numeric(tik)) /n
      param$beta <- coefficients(tmp)
    }else{
      tmp <- rq(Y~. , weights = as.numeric(tik), data = covariates[,-ncol(covariates)], tau = param$tau)
      resid <- tmp$residuals
      param$tunereg <- - sum( (resid * ((resid<0) - param$tau))  * as.numeric(tik)) /n
      param$beta <- c(coefficients(tmp)[-1], coefficients(tmp)[1])
    }
  }
  if (is.null(u)){
    names(param$beta) <- paste0("coeff.Class", 1:K)
  }else{
    names(param$beta) <- c(paste0("coeff.U",1:ncol(u)), paste0("coeff.Class", 1:K))
  }
  param
}

Singleblock_algo_Param_onealgo <- function(y, x, u, model.reg, K, tau, param, tol){
  logPDFwithZ <- Singleblock_algo_Param_PDF(y, x, u, model.reg, param)
  logPDF <- logRowSums(logPDFwithZ)
  loglike <- sum(logPDF)
  repeat{
    tik <- exp(sweep(logPDFwithZ, 1, logPDF, "-"))
    param <- Singleblock_algo_Param_Mstep(y, x, u, model.reg, tik, param)
    logPDFwithZ <- Singleblock_algo_Param_PDF(y, x, u, model.reg, param)
    logPDF <- logRowSums(logPDFwithZ)
    prec <- loglike
    loglike <- sum(logPDF)
    if (is.nan(loglike)){
      loglike <- prec <- -Inf
      break
    }
    if ( (loglike - prec) < tol) break
  }
  list(param = param, loglike = loglike, tik = tik,
       zhat = apply(tik, 1, which.max))
}


Singleblock_algo_Param <- function(y, x, u, model.reg, nbinit, K, tau, nbCPU, tol) {
  all.init <- replicate(nbinit, Singleblock_algo_Param_init(y, x, u, K, tau), simplify = FALSE)
  all.res <- mclapply(all.init, function(param) Singleblock_algo_Param_onealgo(y, x, u, model.reg, K, tau, param, tol), mc.cores = nbCPU)
  all.res <- all.res[[which.max(sapply(all.res, function(u) u$loglike))]]
  all.res
}
