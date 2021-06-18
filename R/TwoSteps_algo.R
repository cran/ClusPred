oneEstimEMTwoSteps <- function(x, K, tau, param, tol){
  logPDFwithZ <- computelogPDFwithZTwoSteps(x, param)
  logPDF <- logRowSums(logPDFwithZ)
  loglike <- sum(logPDF)
  repeat{
    tik <- exp(sweep(logPDFwithZ, 1, logPDF, "-"))
    param <- MstepTwoSteps(x, tik, param)
    logPDFwithZ <- computelogPDFwithZTwoSteps(x, param)
    logPDF <- logRowSums(logPDFwithZ)
    prec <- loglike
    loglike <- sum(logPDF)
    if (is.nan(loglike)){
      loglike <- prec <- -Inf
      break
    }
    if ( (loglike - prec) < tol) break
  }
  list(param = param, loglike = loglike, tik = tik, zhat = apply(tik, 1, which.max))
}

oneEstimEMNPTwoSteps <- function(x, K, tau, param, tol){
  logSmoothPDFwithZ <- computeSmoothlogPDFwithZTwoSteps(x, param)
  logSmoothPDF <- logRowSums(logSmoothPDFwithZ)
  logSmoothlike <- sum(logSmoothPDF)
  repeat{
    param$weights <- tik <- exp(sweep(logSmoothPDFwithZ, 1, logSmoothPDF, "-"))
    param <- MstepNPTwoSteps(x, tik, param)
    logSmoothPDFwithZ <- computeSmoothlogPDFwithZTwoSteps(x, param)
    logSmoothPDF <- logRowSums(logSmoothPDFwithZ)
    prec <- logSmoothlike
    logSmoothlike <- sum(logSmoothPDF)
    if (is.nan(logSmoothlike)){
      logSmoothlike <- prec <- -Inf
      break
    }
    if ( (logSmoothlike - prec) < tol) break
  }
  list(param = param, logSmoothlike = logSmoothlike, tik = tik, zhat = apply(tik, 1, which.max))
}

TwostepsAlgo <- function(y, x, u, np, model.reg, nbinit, K, tau, nbCPU, tol, band){
  results <- list()
  # Step 1: clustering
  if (np){
    all.res <- mclapply(replicate(nbinit, Singleblock_algo_NP_init(y, x, u, K, tau, band), simplify = FALSE),
                        function(param) oneEstimEMNPTwoSteps(x, K, tau, param, tol), mc.cores = nbCPU)
    results <- all.res[[which.max(sapply(all.res, function(u) u$logSmoothlike))]]
  }else{
    all.res <- mclapply(replicate(nbinit, Singleblock_algo_Param_init(y, x, u, K, tau), simplify = FALSE),
                        function(param) oneEstimEMTwoSteps(x, K, tau, param, tol), mc.cores = nbCPU)
    results <- all.res[[which.max(sapply(all.res, function(u) u$loglike))]]
  }
  # Step 2: prediction
  n <- length(y)
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
      tmp <- glm(Y~. + 0, weights = as.numeric(results$tik), data = covariates)
      results$param$tunereg <- sum((tmp$residuals**2) * as.numeric(results$tik)) /n
      results$param$beta <- coefficients(tmp)
    }else{
      tmp <- glm(Y~., weights = as.numeric(results$tik), data = covariates[,-ncol(covariates)])
      results$param$tunereg <- sum((tmp$residuals**2) * as.numeric(results$tik)) /n
      results$param$beta <- c(coefficients(tmp)[-1], coefficients(tmp)[1])
    }
  }else if (model.reg == "quantile"){
    if (K>1){
      tmp <- rq(Y~. + 0, weights = as.numeric(results$tik), data = covariates, tau = tau)
      resid <- tmp$residuals
      results$param$tunereg <- - sum( (resid * ((resid<0) - results$param$tau))  * as.numeric(results$tik)) /n
      results$param$beta <- coefficients(tmp)
    }else{
      tmp <- rq(Y~. , weights = as.numeric(results$tik), data = covariates[,-ncol(covariates)], tau = tau)
      resid <- tmp$residuals
      results$param$tunereg <- - sum( (resid * ((resid<0) - results$param$tau))  * as.numeric(results$tik)) /n
      results$param$beta <- c(coefficients(tmp)[-1], coefficients(tmp)[1])
    }
  }else if (model.reg == "expectile"){
    results$param <- regexpectile(rep(y, K), covariates, as.numeric(results$tik), list(beta=rep(0, ncol(covariates)), tau=tau))
  }else if (model.reg == "huber"){
    results$param <- reghuber(rep(y, K), covariates, as.numeric(results$tik), list(beta=rep(0, ncol(covariates)), tau=tau))
  }else if (model.reg == "logcosh"){
    results$param <- reglogcosh(rep(y, K), covariates, as.numeric(results$tik),list(beta=rep(0, ncol(covariates)), tau=tau))
  }
  if (is.null(u)){
    names(results$param$beta) <- paste0("coeff.Class", 1:K)
  }else{
    names(results$param$beta) <- c(paste0("coeff.U",1:ncol(u)), paste0("coeff.Class", 1:K))
  }
  return(results)
}

