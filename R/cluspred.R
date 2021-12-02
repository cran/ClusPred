
##' ClusPred.
##'
##'
##'
##' \tabular{ll}{
##'   Package: \tab ClusPred\cr
##'   Type: \tab Package\cr
##'   Version: \tab 1.1.0\cr
##'   Date: \tab 2021-12-01\cr
##'   License: \tab GPL-3\cr
##'   LazyLoad: \tab yes\cr
##' }
##'
##' @description
##' Parameter estimation of regression models with fixed group effects, when the group variable is missing while group-related variables are available.
##' @references Simultaneous semi-parametric estimation of clustering and regression, Matthieu Marbac and Mohammed Sedki and Christophe Biernacki and Vincent Vandewalle (2020) <arXiv:2012.14159>.
##' @name ClusPred-package
##' @aliases ClusPred
##' @rdname ClusPred-package
##' @docType package
##' @keywords package
##' @import parallel
##' @import ALDqr
##' @import ald
##' @import quantreg
##' @import Rcpp
##' @import VGAM
##' @importFrom stats dnorm glm integrate optim runif var weighted.mean na.omit
##' @useDynLib  ClusPred
NULL

##' Simulated data
##'
##'
##'
##'
##' @description simulated data used for the pacakge examples.
##' @name simdata
##' @docType data
##' @keywords datasets
##'
##' @examples
##' data(simdata)
NULL

##' Function used for clustering and fitting the regression model
##'
##' @description
##' Estimation of the group-variable Z based on covariates X and estimation of the parameters of the regression of Y on (U, Z)
##'
##' @param y numeric vector of the traget variable (must be numerical)
##' @param x matrix used for clustering (can contain numerical and factors)
##' @param u matrix of the covariates used for regression (can contain numerical and factors)
##' @param K number of clusters
##' @param model.reg indicates the type of the loss ("mean", "quantile", "expectile", "logcosh", "huber"). Only the losses "mean" and "quantile" are implemented if simultaneous=FALSE or np=FALSE
##' @param tau specifies the level for the loss (quantile, expectile or huber)
##' @param simultaneous oolean indicating whether the clustering and the regression are performed simultaneously (TRUE) or not (FALSE)
##' @param np boolean indicating whether nonparameteric model is used (TRUE) or not (FALSE)
##' @param nbinit number of random initializations
##' @param nbCPU number of CPU only used for linux
##' @param tol to specify the stopping rule
##' @param band bandwidth selection
##' @param seed value of the seed (used for drawing the starting points)
##'
##' @return cluspred returns a list containing the model parameters (param), the posterior probabilities of cluster memberships (tik), the partition (zhat) and the (smoothed) loglikelihood)
##' @references Simultaneous semi-parametric estimation of clustering and regression, Matthieu Marbac and Mohammed Sedki and Christophe Biernacki and Vincent Vandewalle (2020) <arXiv:2012.14159>.
##'
##' @examples
##' require(ClusPred)
##' # data loading
##' data(simdata)
##'
##' # mean regression with two latent groups in parametric framework and two covariates
##' res <- cluspred(simdata$y, simdata$x, simdata$u, K=2,
##'  np=FALSE, nbCPU = 1, nbinit = 10)
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
##' model.reg = "quantile", tau = 0.5, nbinit = 10)
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
cluspred <- function(y, x, u = NULL,
                     K = 2, model.reg = "mean",  tau = 0.5,
                     simultaneous=TRUE, np=TRUE,
                     nbinit = 20, nbCPU = 1, tol = 0.01,
                     band = (length(y)**(-1/5)), seed=134) {
  checkinputs(y, x, u, K, model.reg, tau, simultaneous, np, nbinit, nbCPU, tol, band, seed)
  set.seed(seed)
  if (!np){
    if (!(model.reg%in%c("mean", "quantile"))) stop("model.red must be mean or quantile for the parametric approach")
  }
  if (simultaneous){
    if (np){
      return(Singleblock_algo_NP(y, x, u, model.reg, nbinit, K, tau, nbCPU, tol, band))
    }else{
      return(Singleblock_algo_Param(y, x, u, model.reg, nbinit, K, tau, nbCPU, tol))
    }
  }else{
    return(TwostepsAlgo(y, x, u, np, model.reg, nbinit, K, tau,nbCPU, tol, band))
  }
}
