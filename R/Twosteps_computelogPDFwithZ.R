computelogPDFwithZTwoSteps <- function(x, param){
  K <- length(param$pi)
  out <- matrix(log(param$pi), nrow(x), K, byrow = TRUE)
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
  }
  out
}

computeSmoothlogPDFwithZTwoSteps <- function(x, param){
  K <- length(param$pi)
  n <- nrow(x)
  out <- matrix(log(param$pi), n, K, byrow = TRUE)
  for (k in 1:K){
    for (j in 1:ncol(x)){
      se <- seq(min(na.omit(x[,j])) - 5 * param$band, max(na.omit(x[,j]))+ 5 * param$band, length.out = 1000)
      size <- (max(se) - min(se)) / length(se)
      out[,k] <- out[,k] + as.numeric(obj3Cpp(se, xj=x[,j], weightsR=param$weights[,k],band=param$band)) * size

    }
  }
  out
}
