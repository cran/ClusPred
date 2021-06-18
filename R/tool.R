logRowSums <- function(u){
  tmp <- apply(u, 1, max)
  tmp + log(rowSums(exp(sweep(u, 1, tmp, "-"))))
}

obj2CppExport <- function(uR, xij, xj, weightsR, band)
  obj2Cpp(uR, xij, xj, weightsR, band)


############# Expectile
rhoexpectile <- function(u, tau) ( (1-tau) * (u <= 0) + tau * (u>0)) * u

objexpectile <- function(beta, y, x, weights, tau){
  v <- y - x %*% beta
  L <- abs(tau - (v<=0)) * v * v
  mean(L * weights)
}

gradexpectile <-function(beta, y, x, weights, tau){
  tmp <- rhoexpectile(y - x%*%beta, tau)
  -colMeans(sweep(x, 1, tmp*weights, "*"))
}

regexpectile <- function(y, x, weights, param){
  tmp <-  optim(param$beta, objexpectile, gr = gradexpectile, method = "BFGS", y=y, x=x, tau=param$tau, weights=weights)
  param$beta <- tmp$par
  param$tunereg <- tmp$value
  param
}

############# Huber
rhohuber <- function(u, tau) u * (abs(u) < tau) + sign(u) * tau * (abs(u) >= tau)

objhuber <- function(beta, y, x, weights, tau){
  v <- y - x %*% beta
  L <- (v**2)/2 * (abs(v)<tau) + (tau * abs(v) - tau * tau /2) * (abs(v)>=tau)
  mean(L * weights)
}

gradhuber <-function(beta, y, x, weights, tau){
  tmp <- rhohuber(y - x%*%beta, tau)
  -colMeans(sweep(x, 1, tmp*weights, "*"))
}

reghuber <- function(y, x, weights, param){
  tmp <- optim(param$beta, objhuber, gr = gradhuber, method = "BFGS", y=y, x=x, tau=param$tau, weights=weights)
  param$beta <- tmp$par
  param$tunereg <- tmp$value
  param
}

############# Logcosh
rhologcosh <- function(u) tanh(u)

objlogcosh <- function(beta, y, x, weights)
  mean(log(cosh(y - x %*% beta)) * weights)

gradlogcosh <-function(beta, y, x, weights){
  tmp <- rhologcosh(y - x %*% beta)
  -colMeans(sweep(x, 1, tmp*weights, "*"))
}

reglogcosh <- function(y, x, weights,param){
  tmp <- optim(param$beta, objlogcosh, gr = gradlogcosh, method = "BFGS", y=y, x=x, weights=weights)
  param$beta <- tmp$par
  param$tunereg <- tmp$value
  param
}

## checks
checkinputs <- function(y, x, u, K, model.reg, tau, simultaneous, np, nbinit, nbCPU, tol, band, seed){
  if (!is.numeric(y)) stop("y must be numeric")
  if ((!is.matrix(x))&(!is.data.frame(x))) stop("x must be data.frame or matrix")
  if (length(y) != nrow(x)) stop("length(y) must be equal to nrow(x)")
  if ((!is.null(u))&(!is.matrix(u))&(!is.data.frame(u))) stop("u must be null or matrix or data.frame")
  if (length(K)!=1) stop("The number of clusters must be length one")
  if ((K<1)|(K!=ceiling(K))) stop("The number of clusters must be a positive integer of length one")
  if (!is.logical(simultaneous)) stop("simultaneous must be a logical of length one")
  if (!is.logical(np)) stop("np must be a logical of length one")
  if (np){
    if (model.reg %in% c("mean", "quantile", "logcosh", "huber", "expectile") == FALSE)
      stop("model.reg must be equal to mean, quantile, logcosh, huber or expectile for the non-parametric approach")
  }else{
    if (model.reg %in% c("mean", "quantile") == FALSE)
      stop("model.reg must be equal to mean or quantile the approaches that are not parametric approach")
  }
}
