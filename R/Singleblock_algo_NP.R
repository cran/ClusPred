Singleblock_algo_NP_init <- function(y, x, u, K, tau, band){
  # proportions
  pi <- runif(K)
  pi <- pi / sum(pi)
  # parameters of clustering
  n <- nrow(x)
  weights <- matrix(0.5/n, n, K)
  see <- sample(1:n, K)
  for (k in 1:K) weights[see[k],k] <- 0.5 + 0.5/n
  # parameters of prediction
  dv <- ifelse(!is.null(u), ncol(u), 0) + K
  beta <- runif(dv, min = -1, max = 1)
  list(pi = pi,
       weights = weights,
       band = band,
       beta = beta,
       tunereg = 1,
       tau = tau)
}

Singleblock_algo_NP_smoothPDF <- function(y, x, u, model.reg, param){
  K <- length(param$pi)
  all.weights <- unlist(param$weights)
  out <- matrix(log(param$pi), length(y), K, byrow = TRUE)
  n <- length(y)
  for (k in 1:K){
    # parameters of clustering
    for (j in 1:ncol(x)){
      who <- which(!is.na(x[,j]))
      if (is.factor(x[,j])){
        alpha <- rep(0, nlevels(x[,j]))
        for (h in 1:nlevels(x[,j])){
          alpha[h] <- sum(param$weights[which(x[,j]==h), k])
        }
        alpha <- alpha / sum(alpha)
        out[who,k] <- out[who,k] + log(alpha[x[who,j]])
      }else{
        se <- seq(min(na.omit(x[,j])) - 5 * param$band, max(na.omit(x[,j]))+ 5 * param$band, length.out = 1000)
        size <- (max(se) - min(se)) / length(se)
        out[who,k] <- out[who,k] + as.numeric(obj3Cpp(se, xj=x[who,j], weightsR=param$weights[who,k],band=param$band)) * size
      }
    }
  }
  epsilon <- list()
  if (!is.null(u)){
    epsilon <-  lapply(1:K, function(k) y -  u%*%param$beta[1:ncol(u)] - param$beta[k + ncol(u)])
  }else{
    epsilon <-  lapply(1:K, function(k) y -  param$beta[k])
  }
  all.residuals <- unlist(epsilon)
  se <- seq(min(all.residuals) - 5 * param$band, max(all.residuals)+ 5 * param$band, length.out = 1000)
  size <- (max(se) - min(se)) / length(se)
  out <- out + matrix(as.numeric(obj3Cpp(se, xj=all.residuals, weightsR=as.numeric(param$weights),band=param$band)) * size, ncol = K, nrow=length(y))
  out
}


Singleblock_algo_NP_Mstep <- function(y, x, u, model.reg, tik, param){
  # proportions
  K <- ncol(tik)
  n <- length(y)
  param$pi <- colSums(tik) / n
  covariates <- u
  if (K>1){
    for (k in 2:K) covariates <- rbind(covariates, u)
  }
  if (model.reg %in% c("mean", "quantile")){
    covariates <- cbind.data.frame(covariates, factor(rep(1:K, each = n)))
    covariates <- cbind.data.frame(rep(y, K), covariates)
    colnames(covariates)[1] <- "Y"
  }else{
    for (k in 1:K){
      loc <- rep(0, n*K)
      loc[(n * (k-1)) + 1:n] <- 1
      covariates <- cbind(covariates, loc)
    }
  }
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
  }else if (model.reg == "expectile"){
    param <- regexpectile(rep(y, K), covariates, as.numeric(tik), param)
  }else if (model.reg == "huber"){
    param <- reghuber(rep(y, K), covariates, as.numeric(tik), param)
  }else if (model.reg == "logcosh"){
    param <- reglogcosh(rep(y, K), covariates, as.numeric(tik), param)
  }
  if (is.null(u)){
    names(param$beta) <- paste0("coeff.Class", 1:K)
  }else{
    names(param$beta) <- c(paste0("coeff.U",1:ncol(u)), paste0("coeff.Class", 1:K))
  }
  param
}

Singleblock_algo_NP_onealgo <- function(y, x, u, model.reg, K, tau, param, tol){
  logSmoothPDFwithZ <- Singleblock_algo_NP_smoothPDF(y, x, u, model.reg, param)
  if (any((logSmoothPDFwithZ)=="NaN")) logSmoothPDFwithZ <-computeSmoothlogPDFwithZTwoSteps(scale(x), param)
  logSmoothPDF <- logRowSums(logSmoothPDFwithZ)
  logSmoothlike <- sum(logSmoothPDF)
  repeat{
    param$weights <- tik <- exp(sweep(logSmoothPDFwithZ, 1, logSmoothPDF, "-"))
    param <- Singleblock_algo_NP_Mstep(y, x, u, model.reg, tik, param)
    logSmoothPDFwithZ <- Singleblock_algo_NP_smoothPDF(y, x, u, model.reg, param)
    logSmoothPDF <- logRowSums(logSmoothPDFwithZ)
    prec <- logSmoothlike
    logSmoothlike <- sum(logSmoothPDF)
    if (is.nan(logSmoothlike)){
      logSmoothlike <- prec <- -Inf
      break
    }
    if ( (logSmoothlike - prec) < tol){
      break
    }
  }
  list(param = param, logSmoothlike = logSmoothlike, tik = tik,
        zhat = apply(tik, 1, which.max))
}

Singleblock_algo_NP <- function(y, x, u, model.reg, nbinit, K, tau, nbCPU, tol, band) {
  all.init <- replicate(nbinit, Singleblock_algo_NP_init(y, x, u, K, tau, band), simplify = FALSE)
  all.res <- mclapply(all.init, function(param) Singleblock_algo_NP_onealgo(y, x, u, model.reg, K, tau, param, tol), mc.cores = nbCPU)
  all.res <- all.res[[which.max(sapply(all.res, function(u) u$logSmoothlike))]]
  all.res
}
