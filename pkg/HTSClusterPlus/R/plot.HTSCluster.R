#' @importFrom graphics legend
#############################################
##    graphe comparaison ICL sur donnees transf sur PP
plot.HTSCluster <- function(resarcsin, reslogit, reslogMedianRef, y_profiles) {
  n=dim(y_profiles)[1]
  p=dim(y_profiles)[2]
  qarcsin<-(n*p*log(2)) + (0.5*sum(sum(log(y_profiles*(1-y_profiles)))))
  qlogit<-(n*p*log(log(2))) + (sum(sum(log(y_profiles*(1-y_profiles)))))
  qlogmedianref<-(n*p*log(log(2))) + sum(sum(log(y_profiles)))
  
  plot(resarcsin$nbClust,resarcsin$ICLvalue + (2*qarcsin),col="red",type="l",
       xlab="nb cluster",ylab="ICL sur les pi")
  points(reslogit$nbClust,reslogit$ICLvalue + (2*qlogit),col="blue",type="l")
  points(reslogMedianRef$nbClust, reslogMedianRef$ICLvalue + (2*qlogmedianref),
         col="magenta",type="l")
  legend("topright",legend=c("arcsin","logit","logMedianRef"),
         col=c("red","blue","magenta"),lty=1)
}