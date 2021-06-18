##' Prediction (clustering and target variable)
##'
##' @description
##' Prediction for new observations
##'
##' @param x covariates used for clustering
##' @param u covariates of the regression (can be null)
##' @param result results provided by function cluspred
##' @param np boolean indicating whether nonparametric estimation is used (TRUE) or not (FALSE)
##'
##' @return predictboth returns a list containing the predicted cluster membership (zhat) and the predicted value of the target variable (yhat).
##' @examples
##' require(ClusPred)
##' # data loading
##' data(simdata)
##'
##' # mean regression with two latent groups in parametric framework and two covariates
##' res <- cluspred(simdata$y, simdata$x, simdata$u, K=2,
##' np=FALSE, nbCPU = 1, nbinit = 10)
##' # coefficient of the regression
##' res$param$beta
##' # proportions of the latent groups
##' res$param$pi
##' # posterior probability of the group memberships
##' head(res$tik)
##' # partition
##' res$zhat
##' # loglikelihood
##' res$loglike
##' # prediction (for possible new observations)
##' pred <- predictboth(simdata$x, simdata$u, res, np = FALSE)
##' # predicted cluster membreships
##' pred$zhat
##' # predicted value of the target variable
##' pred$yhat
##'
##' \donttest{
##' # median regression with two latent groups in nonparametric framework and two covariates
##' res <- cluspred(simdata$y, simdata$x, simdata$u, K=2,
##'  model.reg = "quantile", tau = 0.5, nbinit = 10)
##' # coefficient of the regression
##' res$param$beta
##' # proportions of the latent groups
##' res$param$pi
##' # posterior probability of the group memberships
##' head(res$tik)
##' # partition
##' res$zhat
##' # smoothed loglikelihood
##' res$logSmoothlike
##' # prediction (for possible new observations)
##' pred <- predictboth(simdata$x, simdata$u, res, np = TRUE)
##' # predicted cluster membreships
##' pred$zhat
##' # predicted value of the target variable
##' pred$yhat
##' }
##'
##'
##' @export
predictboth <- function(x, u=NULL, result, np=FALSE){
  K <- length(result$param$pi)
  out <- matrix(log(result$param$pi), nrow(x), K, byrow = TRUE)
  if(!np){
    for (k in 1:K){
      # parameters of clustering
      for (j in 1:ncol(x)){
        obs <- which(!is.na(x[,j]))
        if (is.factor(x[,j])){
          out[obs,k] <- out[obs,k] + log(result$param$alpha[[j]][k,x[obs,j]])
        }else{
          out[obs,k] <- out[obs,k] + dnorm(x[obs,j], result$param$alpha[[j]][k,1], sqrt(result$param$alpha[[j]][k,2]), log = TRUE)
        }
      }
    }
  }else{
    all.weights <- unlist(result$param$weights)
    for (k in 1:K){
      # parameters of clustering
      for (j in 1:ncol(x)){
        who <- which(!is.na(x[,j]))
        if (is.factor(x[,j])){
          alpha <- rep(0, nlevels(x[,j]))
          for (h in 1:nlevels(x[,j])){
            alpha[h] <- sum(result$param$weights[which(x[,j]==h), k])
          }
          alpha <- alpha / sum(alpha)
          out[who,k] <- out[who,k] + log(alpha[x[who,j]])
        }else{
          se <- seq(min(na.omit(x[,j])) - 5 * result$param$band, max(na.omit(x[,j]))+ 5 * result$param$band, length.out = 1000)
          size <- (max(se) - min(se)) / length(se)
          out[,k] <- out[,k] + as.numeric(obj3Cpp(se, xj=x[who,j], weightsR=result$param$weights[who,k], band=result$param$band)) * size
        }
      }
    }
  }
  ztest <- unlist(apply(out, 1, which.max))
  yhat <- result$param$beta[ztest]
  if (is.null(u)==FALSE){
    u <- as.matrix(u)
    yhat <- as.matrix(u) %*% result$param$beta[1:ncol(u)] + result$param$beta[ncol(u) + ztest]
  }
  list(zhat = ztest, yhat = as.numeric(yhat))
}

