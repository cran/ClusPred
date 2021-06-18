MstepTwoSteps <- function(x, tik, param){
  # proportions
  K <- ncol(tik)
  n <- nrow(x)
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
  param
}

MstepNPTwoSteps <- function(x, tik, param){
  # proportions
  K <- ncol(tik)
  n <- nrow(x)
  param$pi <- colSums(tik) / n
  # parameters of clustering
  param
}
