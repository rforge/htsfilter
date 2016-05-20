#' @importFrom graphics matplot
#' @importFrom grDevices heat.colors
#' 
# plot.NormMixClus_K <- function(x, y_profiles, K, threshold=0.8, conds=NULL, 
#                                average_conditions = TRUE, 
#                                graphs=c("profiles_K", "boxplots_K", 
#                                         "profiles", "boxplots"), ...) {
# 
#   labels <- apply(x$ICL.results$probaPost, 1, which.max)
#   proba <- apply(x$ICL.results$probaPost, 1, max)
#   
#   if("profiles_K" %in% graphs) {
#     I <- which(labels==K)
#     I1 <- which(labels==K & proba > threshold)
#     par(mfrow=c(1,2))
#     matplot(t(y_profiles[I1,]),type="l",col="black",main=paste("Cluster ",K,sep=""),ylab="")
#     matplot(t(y_profiles[I,]),type="l",col="black",main=paste("Cluster ",K,sep=""),ylab="") 
#   }
#   
#   if("boxplots_K" %in% graphs) {
#     I <- which(labels==K)
#     colo <- NULL
#     if(is.null(conds)==FALSE){    
#       pal<-heat.colors(length(unique(conds)))
#       a <- tabulate(conds)
#       for(aa in 1:length(a)){
#         colo <- c(colo,rep(pal[aa],a[aa]))
#       }
#     }
#     boxplot(y_profiles[I,], ylim=c(0,max(y_profiles)), 
#             main=paste("Cluster ",cluster," (",length(I),")",sep=""), col=colo)
#     points(apply(y_profiles[I,],2,mean), type="l", col="red")
#   }
#   
#   if("profiles" %in% graphs) {
#     
#     
#   }
#   
#   if("boxplots" %in% graphs) {
# 
#     orderCluster <- seq(1, max(labels))
#     y_profiles_plot <- y_profiles
#     if (average_conditions == TRUE & is.null(conds)==FALSE) {     
#       a <- unique(conds)
#       A <- apply(y_profiles[,which(conds==a[1])],1,mean)
#       if(length(a)>1) {
#         for(k in 2:length(a)) {
#           A <- cbind(A, apply(y_profiles[,which(conds==a[k])],1,mean))
#         }
#       }
#       y_profiles_plot <- A
#       conds<-seq(1,length(unique(conds)))
#     }
#     v <- floor(max(labels)/6)
#     if (v>0){
#       for (u in 1:v){
#         op <- par(mfrow=c(2,3),mar=1.8*rep(1,4))  
#         for (k in ((u-1)*6 +1):(u*6)){
#           if (is.null(conds)==T){
#             BoxplotProfilPP(Res,PP,cluster=ordercluster[k])  
#          }else{BoxplotProfilPP(Res,PP,cluster=ordercluster[k],conds=conds)}
#         }    
#       }
#     }
#     if (max(label) > (6*v)){
#       op<-par(mfrow=c(3,3),mar=1.8*rep(1,4))  
#       for (k in (v*6 +1):max(label)){
#         if (is.null(conds)==T){
#           BoxplotProfilPP(Res,PP,cluster=ordercluster[k])  
#         }else{BoxplotProfilPP(Res,PP,cluster=ordercluster[k],conds=conds)}
#       }   
#     }   
#     
#   }
# }
# 
# PPmoy<-function(PP, conds){
#   a<-unique(conds)
#   A<-apply(PP[,which(conds==a[1])],1,mean)
#   if (length(a)>1){
#     for (k in 2:length(a))
#       A<-cbind(A,apply(PP[,which(conds==a[k])],1,mean))
#   }
#   return(A)
# }

# graphprofilsPP<-function(Res,PP,ordercluster=NULL,conds=NULL,group=F){
#   label<-apply(Res$ICL.results$probaPost,1,which.max)
#   proba<-apply(Res$ICL.results$probaPost,1,max)
#   if (is.null(ordercluster)==T){
#     ordercluster<-seq(1,max(label))
#   }
#   if (group==T & is.null(conds)==F){     #il faut rajouter des sécurités sur conds, group=T nécessite conds= un bon vecteur, ....
#     PP<-PPmoy(PP,conds)
#     conds<-seq(1,length(unique(conds)))
#   }
#   v<-floor(max(label)/6)
#   if (v>0){
#     for (u in 1:v){
#       op<-par(mfrow=c(2,3),mar=1.8*rep(1,4))  
#       for (k in ((u-1)*6 +1):(u*6)){
#         if (is.null(conds)==T){
#           BoxplotProfilPP(Res,PP,cluster=ordercluster[k])  
#         }else{BoxplotProfilPP(Res,PP,cluster=ordercluster[k],conds=conds)}
#       }    
#     }
#   }
#   if (max(label) > (6*v)){
#     op<-par(mfrow=c(3,3),mar=1.8*rep(1,4))  
#     for (k in (v*6 +1):max(label)){
#       if (is.null(conds)==T){
#         BoxplotProfilPP(Res,PP,cluster=ordercluster[k])  
#       }else{BoxplotProfilPP(Res,PP,cluster=ordercluster[k],conds=conds)}
#     }   
#   }   
# }