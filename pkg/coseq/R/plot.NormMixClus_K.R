#' Plot NormMixClus_K object
#' 
#' Blah blah blah
#' 
#' @param x An object of class \code{"NormMixClus_K"}
#' @param y_profiles y (\emph{n} x \emph{q}) matrix of observed profiles for \emph{n}
#' observations and \emph{q} variables to be used for graphing
#' @param K If desired, the specific cluster number to use for plotting. If \code{NULL},
#' all clusters will be visualized
#' @param threshold Threshold used for maximum conditional probability; only observations
#' with maximum conditional probability greater than this threshold are visualized 
#' @param conds Condition labels, if desired
#' @param average_over_conds If \code{TRUE}, average values of \code{y_profiles} within
#' each condition identified by \code{conds}
#' @param graphs Graphs to be produced, one (or more) of the following: 
#' \code{"profiles_K"}, ...
#' @param ...  Additional optional plotting arguments
#'
#' @importFrom graphics matplot
#' @importFrom grDevices heat.colors
#' @importFrom scales alpha
#' @import ggplot2
 
### TODO: average over conditions option
### TODO: save files 
plot.NormMixClus_K <- function(x, y_profiles, K=NULL, threshold=0.8, conds=NULL,
                                average_over_conds=TRUE,
                                graphs=c("profiles", "boxplots"), ...) {
  
  
  labels <- apply(x$probaPost, 1, which.max)
  proba <- apply(x$probaPost, 1, max)
  
  if(length(conds) > 0) {
    conds <- as.factor(conds);
    conds_vec <- rep(conds, each=nrow(y_profiles))
  }
  if(length(conds) == 0) {
    conds_vec <- rep(NA, nrow(y_profiles)*ncol(y_profiles))
  }
  
  ## Parse ellipsis function
  arg.user <- list(...)
  if(is.null(arg.user$alpha)) arg.user$alpha<-0.3;
  
  
  if(average_over_conds == FALSE) {
    pl_data <- data.frame(ID=ifelse(rep(length(rownames(y_profiles))==0, nrow(y_profiles)), 
                                    rep(1:nrow(y_profiles), times=ncol(y_profiles)),
                                    rownames(y_profiles)),
                          y_prof=as.vector(y_profiles), 
                          col_num=rep(1:ncol(y_profiles), each=nrow(y_profiles)),
                          col_nam=rep(colnames(y_profiles), each=nrow(y_profiles)),
                          conds=conds_vec,
                          labels=rep(labels, times=ncol(y_profiles)),
                          proba=rep(proba, times=ncol(y_profiles)))
  }
  
  if(average_over_conds == TRUE) {
    if(length(conds) == 0) stop("Conds argument needed when average_over_conds == TRUE")
    y_profiles_c <- t(rowsum(t(y_profiles), conds))
    conds_vec <- factor(rep(colnames(y_profiles_c), each=nrow(y_profiles_c)))
    pl_data <- data.frame(ID=ifelse(rep(length(rownames(y_profiles_c))==0, nrow(y_profiles_c)), 
                                    rep(1:nrow(y_profiles_c), times=ncol(y_profiles_c)),
                                    rownames(y_profiles_c)),
                          y_prof=as.vector(y_profiles_c), 
                          col_num=rep(1:ncol(y_profiles_c), each=nrow(y_profiles_c)),
                          col_nam=rep(colnames(y_profiles_c), each=nrow(y_profiles_c)),
                          conds=conds_vec,
                          labels=rep(labels, times=ncol(y_profiles_c)),
                          proba=rep(proba, times=ncol(y_profiles_c)))
  }
  
  
  if("profiles" %in% graphs) {
    if(average_over_conds == FALSE) {
      ## For one specific value of K
      if(is.null(K) == FALSE) {
        pl_data_tmp <- pl_data[which(pl_data$labels == K),]
        g1 <- ggplot(pl_data_tmp[which(pl_data_tmp$proba > threshold),]) +
          geom_line(colour=alpha("black", arg.user$alpha), 
                    aes_string(x="col_num", y="y_prof", group="ID")) +
          geom_line(data=pl_data_tmp[which(pl_data_tmp$proba < threshold),],
                    colour=alpha("red", arg.user$alpha), 
                    aes_string(x="col_num", y="y_prof", group="ID")) +
          theme_bw() + ggtitle(paste("Cluster", K)) +
          scale_y_continuous(name="y") + scale_x_continuous(name="Sample number")
        
        print(g1)
      }
      ## For all values of K
      if(is.null(K) == TRUE) {
        g2 <- ggplot(pl_data[which(pl_data$proba > threshold),]) +
          geom_line(colour=alpha("black", arg.user$alpha), 
                    aes_string(x="col_num", y="y_prof", group="ID")) +
          geom_line(data=pl_data[which(pl_data$proba < threshold),],
                    colour=alpha("red", arg.user$alpha), 
                    aes_string(x="col_num", y="y_prof", group="ID")) +
          theme_bw() +
          scale_y_continuous(name="y") + scale_x_continuous(name="Sample number") +
          facet_wrap(~labels)
        
        print(g2)
      }
    }
    if(average_over_conds == TRUE) {
      ## For one specific value of K
      if(is.null(K) == FALSE) {
        pl_data_tmp <- pl_data[which(pl_data$labels == K),]
        g1b <- ggplot(pl_data_tmp[which(pl_data_tmp$proba > threshold),]) +
          geom_line(colour=alpha("black", arg.user$alpha), 
                    aes_string(x="conds", y="y_prof", group="ID")) +
          geom_line(data=pl_data_tmp[which(pl_data_tmp$proba < threshold),],
                    colour=alpha("red", arg.user$alpha), 
                    aes_string(x="conds", y="y_prof", group="ID")) +
          theme_bw() + ggtitle(paste("Cluster", K)) +
          scale_y_continuous(name="Average y") + scale_x_discrete(name="Conditions")
        
        print(g1b)
      }
      ## For all values of K
      if(is.null(K) == TRUE) {
        g2b <- ggplot(pl_data[which(pl_data$proba > threshold),]) +
          geom_line(colour=alpha("black", arg.user$alpha), 
                    aes_string(x="conds", y="y_prof", group="ID")) +
          geom_line(data=pl_data[which(pl_data$proba < threshold),],
                    colour=alpha("red", arg.user$alpha), 
                    aes_string(x="conds", y="y_prof", group="ID")) +
          theme_bw() +
          scale_y_continuous(name="Average y") + scale_x_discrete(name="Conditions") +
          facet_wrap(~labels)
        
        print(g2b)
      }
    }
  
  }
  
  
  if("boxplots" %in% graphs) {
    ## For one specific value of K
    pl_data_tmp <- pl_data[which(pl_data$labels == K),]
    pl_data_tmp$col_num <- factor(pl_data_tmp$col_num)
    pl_data_tmp$conds <- factor(pl_data_tmp$conds)
    if(average_over_conds == FALSE) {
      if(is.null(K) == FALSE) {
        if(length(conds)==0) {
          g3 <- ggplot(pl_data_tmp, 
                       aes_string(x="col_num", y="y_prof")) +
            geom_boxplot() +
            stat_summary(fun.y=mean, geom="line", aes(group=1), colour="red")  + 
            stat_summary(fun.y=mean, geom="point", colour="red") +
            ggtitle(paste("Cluster", K)) +
            scale_y_continuous(name="y") + scale_x_discrete(name="Sample number")
          print(g3)
        }
        if(length(conds)>0) {
          g4 <- ggplot(pl_data_tmp, 
                       aes_string(x="col_num", y="y_prof")) +
            geom_boxplot(aes_string(fill="conds")) +
            stat_summary(fun.y=mean, geom="line", aes(group=1), colour="red")  + 
            stat_summary(fun.y=mean, geom="point", colour="red") +
            ggtitle(paste("Cluster", K)) +
            scale_y_continuous(name="y") + scale_x_discrete(name="Sample number") +
            scale_fill_discrete(name="Conditions")
          print(g4)
        }
      }
      if(is.null(K) == TRUE) {
        if(length(conds)==0) {
          g5 <- ggplot(pl_data, aes_string(x="col_num", y="y_prof")) +
            geom_boxplot() +
            stat_summary(fun.y=mean, geom="line", aes(group=1), colour="red")  + 
            stat_summary(fun.y=mean, geom="point", colour="red") +
            scale_y_continuous(name="y") + scale_x_discrete(name="Sample number") +
            facet_wrap(~labels)
          print(g5)
        }
        if(length(conds)>0) {
          g6 <- ggplot(pl_data, aes_string(x="col_num", y="y_prof")) +
            geom_boxplot(aes_string(fill="conds")) +
            stat_summary(fun.y=mean, geom="line", aes(group=1), colour="red")  + 
            stat_summary(fun.y=mean, geom="point", colour="red") +
            facet_wrap(~labels) + 
            scale_y_continuous(name="y") + scale_x_discrete(name="Sample number") +
            scale_fill_discrete(name="Conditions")
          print(g6)
        }
      }
    }
    if(average_over_conds == TRUE) {
      if(is.null(K) == FALSE) {
        g7 <- ggplot(pl_data_tmp, aes_string(x="conds", y="y_prof")) +
          geom_boxplot(aes_string(fill="conds")) +
          stat_summary(fun.y=mean, geom="line", aes(group=1), colour="red")  + 
          stat_summary(fun.y=mean, geom="point", colour="red") +
          ggtitle(paste("Cluster", K)) +
          scale_y_continuous(name="Average y") + scale_x_discrete(name="Conditions") +
          scale_fill_discrete(name="Conditions")
        print(g7)
      }
      
      if(is.null(K) == TRUE) {
        g8 <- ggplot(pl_data, aes_string(x="cibds", y="y_prof")) +
          geom_boxplot(aes_string(fill="conds")) +
          stat_summary(fun.y=mean, geom="line", aes(group=1), colour="red")  + 
          stat_summary(fun.y=mean, geom="point", colour="red") +
          facet_wrap(~labels) + 
          scale_y_continuous(name="Average y") + scale_x_discrete(name="Conditions") +
          scale_fill_discrete(name="Conditions")
        print(g8)
      }
    }
  }
}




# 
# ## boxplot des probapost
# if("boxplots" %in% graphs) {
#   probapost <- x$ICL.results$probaPost
#   label <- apply(probapost,1,which.max)
#   A <- boxplot(apply(probapost,1,max)~label, plot=F)
#   if(order==TRUE){
#     J <- sort.int(A$stat[3,],index.return=T,decreasing=T)$ix   #ordre par rapport à la mediane
#   } else{
#     J<-seq(1,ncol(probapost),1)
#   }  
#   boxplot(A$stat[,J],axes=F,outline=T)
#   axis(side=1, at=seq(1,ncol(probapost)), labels=J, cex.lab=0.5)
#   axis(side=2)
# }
# 
# ## graphe barplot de la proportion de probapost >seuil et < seuil par classe
# if("barplots" %in% graphs) {
#   A<-matrix(0, nrow=2, ncol=ncol(probapost))
#   for (k in 1:ncol(probapost)){
#     I<-which(label==k)
#     A[,k]<-c(sum(probapost[I,k]>threshold), sum(probapost[I,k]<=threshold))
#   }
#   barplot(A[,J],names.arg=J)
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
 # 
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
 #    }
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