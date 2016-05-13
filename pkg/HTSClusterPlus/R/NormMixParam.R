#' Calculate the mean and variance parameters for a normal mixture model
#' 
#' Corresponds to pK_Lk_Ck model
#'
#' @param y_profiles y (\emph{n} x \emph{q}) matrix of observed profiles for \emph{n}
#' observations and \emph{q} variables
#' @param probaPost Matrix containing the conditional
#' probabilities of belonging to each cluster for all observations
#'
#' @return
#' \item{pi }{ ...}
#' \item{mu }{ ...}
#' \item{Sigma }{ ...}
#' 
#' @export
#'

##   il faut supprimer les boucles for, à améliorer
NormMixParam <- function(y_profiles, probaPost){
  pi <- apply(probaPost,2,sum)/nrow(y_profiles)
  mu <- matrix(0,nrow=ncol(probaPost),ncol=ncol(y_profiles))
  Sigma <- array(0,dim=c(ncol(y_profiles),ncol(y_profiles),ncol(probaPost)))
  for (k in 1:ncol(probaPost)){
    mu[k,]<-apply(probaPost[,k] * y_profiles,2,sum) / sum(probaPost[,k])
    for (i in 1:nrow(y_profiles)){        # a ameliorer ici
      Sigma[,,k]<-Sigma[,,k] + (probaPost[i,k]*(t(as.matrix(y_profiles[i,]-mu[k,])) %*% 
                                                  as.matrix(y_profiles[i,]-mu[k,])))
    }
    Sigma[,,k]<-Sigma[,,k] / sum(probaPost[,k])
  }
  param<-list(pi=pi, mu=mu, Sigma=Sigma)
  return(param)           
}
