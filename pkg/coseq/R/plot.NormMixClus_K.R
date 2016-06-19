#' Plot NormMixClus_K object
#' 
#' Plot a NormMixClus_K object.
#' 
#' @param x An object of class \code{"NormMixClus_K"}
#' @param y_profiles y (\emph{n} x \emph{q}) matrix of observed profiles for \emph{n}
#' observations and \emph{q} variables to be used for graphing
#' @param K If desired, the specific cluster number(s) to use for plotting. If \code{NULL},
#' all clusters will be visualized
#' @param threshold Threshold used for maximum conditional probability; only observations
#' with maximum conditional probability greater than this threshold are visualized 
#' @param conds Condition labels, if desired
#' @param average_over_conds If \code{TRUE}, average values of \code{y_profiles} within
#' each condition identified by \code{conds}
#' @param graphs Graphs to be produced, one (or more) of the following: 
#' \code{"profiles"} (line plots of profiles in each cluster), \code{"boxplots"} 
#' (boxplots of profiles in each cluster), \code{"probapost_boxplots"} (boxplots of
#' maximum conditional probabilities per cluster), \code{"probapost_barplots"} 
#' (number of observations with a maximum conditional probability greater than 
#' \code{threshold} per cluster), \code{"probapost_histogram"} (histogram of maximum
#' conditional probabilities over all clusters) ...
#' @param order If \code{TRUE}, order clusters in \code{probapost_boxplot} by median and
#' \code{probapost_barplot} by number of observations with maximum conditional probability
#' greater than \code{threshold}
#' @param ...  Additional optional plotting arguments
#'
#' @importFrom graphics matplot boxplot
#' @importFrom grDevices heat.colors
#' @importFrom scales alpha
#' @import ggplot2
 
plot.NormMixClus_K <- function(x, y_profiles, K=NULL, threshold=0.8, conds=NULL,
                                average_over_conds=FALSE,
                                graphs=c("profiles", "boxplots",
                                         "probapost_boxplots",
                                         "probapost_barplots",
                                         "probapost_histogram"), 
                               order = FALSE, ...) {
  
  if(is.null(x$probaPost) == TRUE) stop("Selected model is empty.")
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
  
  
  #####################################################
  ## SET UP PLOTTING DATA.FRAME
  #####################################################
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
  
  #####################################################
  ## PROFILE PLOTS
  #####################################################
  if("profiles" %in% graphs) {
    if(average_over_conds == FALSE) {
      ## For one specific value of K
      if(is.null(K) == FALSE & length(K) == 1) {
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
      ## For a subset of values of K
      if(is.null(K) == FALSE & length(K) > 1) {
        pl_data_tmp <- pl_data[which(pl_data$labels %in% K),]
        g2bb <- ggplot(pl_data_tmp[which(pl_data_tmp$proba > threshold),]) +
          geom_line(colour=alpha("black", arg.user$alpha), 
                    aes_string(x="col_num", y="y_prof", group="ID")) +
          geom_line(data=pl_data_tmp[which(pl_data_tmp$proba < threshold),],
                    colour=alpha("red", arg.user$alpha), 
                    aes_string(x="col_num", y="y_prof", group="ID")) +
          theme_bw() +
          scale_y_continuous(name="y") + scale_x_continuous(name="Sample number") +
          facet_wrap(~labels)
        print(g2bb)
      }
    }
    if(average_over_conds == TRUE) {
      ## For one specific value of K
      if(is.null(K) == FALSE & length(K) == 1) {
        pl_data_tmp <- pl_data[which(pl_data$labels == K),]
        g1b <- ggplot(pl_data_tmp[which(pl_data_tmp$proba > threshold),]) +
          geom_line(colour=alpha("black", arg.user$alpha), 
                    aes_string(x="conds", y="y_prof", group="ID")) +
          geom_line(data=pl_data_tmp[which(pl_data_tmp$proba < threshold),],
                    colour=alpha("red", arg.user$alpha), 
                    aes_string(x="conds", y="y_prof", group="ID")) +
          theme_bw() + ggtitle(paste("Cluster", K)) +
          scale_y_continuous(name="Average y") + 
          scale_x_discrete(name="Conditions")
        
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
      ## For a subset of values of K
      if(is.null(K) == FALSE & length(K) > 1) {
        pl_data_tmp <- pl_data[which(pl_data$labels %in% K),]
        g2bb <- ggplot(pl_data_tmp[which(pl_data_tmp$proba > threshold),]) +
          geom_line(colour=alpha("black", arg.user$alpha), 
                    aes_string(x="conds", y="y_prof", group="ID")) +
          geom_line(data=pl_data_tmp[which(pl_data_tmp$proba < threshold),],
                    colour=alpha("red", arg.user$alpha), 
                    aes_string(x="conds", y="y_prof", group="ID")) +
          theme_bw() +
          scale_y_continuous(name="Average y") + scale_x_discrete(name="Conditions") +
          facet_wrap(~labels)
        print(g2bb)
      }
    }
  }
  
  #####################################################
  ## PROFILE BOXPLOTS
  #####################################################
  if("boxplots" %in% graphs) {
    ## For one specific value of K
    pl_data_tmp <- pl_data[which(pl_data$labels %in% K),]
    pl_data_tmp$col_num <- factor(pl_data_tmp$col_num)
    pl_data_tmp$conds <- factor(pl_data_tmp$conds)
    pl_data$col_num <- factor(pl_data$col_num)
    
    if(average_over_conds == FALSE) {
      ## Single value of K
      if(is.null(K) == FALSE & length(K) == 1) {
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
      ## All values of K
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
      ## A subset of values of K
      if(is.null(K) == FALSE & length(K) > 1) {
        if(length(conds)==0) {
          g5b <- ggplot(pl_data_tmp, aes_string(x="col_num", y="y_prof")) +
            geom_boxplot() +
            stat_summary(fun.y=mean, geom="line", aes(group=1), colour="red")  + 
            stat_summary(fun.y=mean, geom="point", colour="red") +
            scale_y_continuous(name="y") + scale_x_discrete(name="Sample number") +
            facet_wrap(~labels)
          print(g5b)
        }
        if(length(conds)>0) {
          g6b <- ggplot(pl_data_tmp, aes_string(x="col_num", y="y_prof")) +
            geom_boxplot(aes_string(fill="conds")) +
            stat_summary(fun.y=mean, geom="line", aes(group=1), colour="red")  + 
            stat_summary(fun.y=mean, geom="point", colour="red") +
            facet_wrap(~labels) + 
            scale_y_continuous(name="y") + scale_x_discrete(name="Sample number") +
            scale_fill_discrete(name="Conditions")
          print(g6b)
        }
      }
    }
    
    if(average_over_conds == TRUE) {
      ## Single value of K
      if(is.null(K) == FALSE & length(K) == 1) {
        g7 <- ggplot(pl_data_tmp, aes_string(x="conds", y="y_prof")) +
          geom_boxplot(aes_string(fill="conds")) +
          stat_summary(fun.y=mean, geom="line", aes(group=1), colour="red")  + 
          stat_summary(fun.y=mean, geom="point", colour="red") +
          ggtitle(paste("Cluster", K)) +
          scale_y_continuous(name="Average y") + scale_x_discrete(name="Conditions") +
          scale_fill_discrete(name="Conditions")
        print(g7)
      }
      ## All values of K
      if(is.null(K) == TRUE) {
        g8 <- ggplot(pl_data, aes_string(x="conds", y="y_prof")) +
          geom_boxplot(aes_string(fill="conds")) +
          stat_summary(fun.y=mean, geom="line", aes(group=1), colour="red")  + 
          stat_summary(fun.y=mean, geom="point", colour="red") +
          facet_wrap(~labels) + 
          scale_y_continuous(name="Average y") + scale_x_discrete(name="Conditions") +
          scale_fill_discrete(name="Conditions")
        print(g8)
      }
      ## A subset of values of K
      if(is.null(K) == FALSE & length(K) > 1) {
        g9 <- ggplot(pl_data_tmp, aes_string(x="conds", y="y_prof")) +
          geom_boxplot(aes_string(fill="conds")) +
          stat_summary(fun.y=mean, geom="line", aes(group=1), colour="red")  + 
          stat_summary(fun.y=mean, geom="point", colour="red") +
          facet_wrap(~labels) + 
          scale_y_continuous(name="Average y") + scale_x_discrete(name="Conditions") +
          scale_fill_discrete(name="Conditions")
        print(g9)
      }
    }
  }
  
  #####################################################
  ## PROBAPOST BOXPLOTS
  #####################################################
  if("probapost_boxplots" %in% graphs) {
    pl_data <- data.frame(ID=ifelse(rep(length(rownames(y_profiles))==0, nrow(y_profiles)), 
                                    rep(1:nrow(y_profiles), times=ncol(y_profiles)),
                                    rownames(y_profiles)),
                          y_prof=as.vector(y_profiles), 
                          col_num=rep(1:ncol(y_profiles), each=nrow(y_profiles)),
                          col_nam=rep(colnames(y_profiles), each=nrow(y_profiles)),
                          conds=conds_vec,
                          labels=rep(labels, times=ncol(y_profiles)),
                          proba=rep(proba, times=ncol(y_profiles)))
    pl_data$labels <- factor(pl_data$labels)
    ## Single value of K or subset of values
    if(is.null(K) == FALSE) {
      pl_data_tmp <- pl_data[which(pl_data$labels %in% K),]
      gg <- ggplot(pl_data_tmp, aes_string(x="labels", y="proba")) +
          geom_boxplot() +  scale_x_discrete(name="Cluster") +
          scale_y_continuous(name="Max conditional probability")
      print(gg)
    }
    ## All K
    if(is.null(K) == TRUE) {
      if(order == TRUE) {
        A <- boxplot(pl_data$proba~pl_data$label, plot=F)
        J <- sort.int(A$stat[3,],index.return=T,decreasing=T)$ix 
        pl_data$labels <- factor(pl_data$labels, levels=J)
      }
      gg <- ggplot(pl_data, aes_string(x="labels", y="proba")) +
        geom_boxplot() +  scale_x_discrete(name="Cluster") +
        scale_y_continuous(name="Max conditional probability")
      print(gg)
    }
  }

  #####################################################
  ## PROBAPOST BARPLOTS
  #####################################################
  if("probapost_barplots" %in% graphs) {
    pl_data <- data.frame(ID=ifelse(rep(length(rownames(y_profiles))==0, nrow(y_profiles)),
                                    rep(1:nrow(y_profiles), times=ncol(y_profiles)),
                                    rownames(y_profiles)),
                          y_prof=as.vector(y_profiles),
                          col_num=rep(1:ncol(y_profiles), each=nrow(y_profiles)),
                          col_nam=rep(colnames(y_profiles), each=nrow(y_profiles)),
                          conds=conds_vec,
                          labels=rep(labels, times=ncol(y_profiles)),
                          proba=rep(proba, times=ncol(y_profiles)))
    pl_data <- pl_data[which(pl_data$col_num==1),]
    pl_data$goodproba <- factor(ifelse(pl_data$proba > threshold,
                                paste(">", threshold), paste("<", threshold)),
                                levels=c(paste(">", threshold),
                                         paste("<", threshold)))
    pl_data$labels <- factor(pl_data$labels)
    ## Single value of K or subset of values
    if(is.null(K) == FALSE) {
      pl_data_tmp <- pl_data[which(pl_data$labels %in% K),]
      gg <- ggplot(pl_data_tmp, aes_string(x="labels", fill="goodproba")) +
        geom_bar() +
        scale_fill_brewer(direction=-1, palette="Accent",
                          name="Max\nconditional\nprobability") +
        scale_x_discrete(name="Cluster") +
        scale_y_continuous(name="Number of observations")
      print(gg)
    }
    ## All K
    if(is.null(K) == TRUE) {
      if(order == TRUE) {
        pl_data$labels <- factor(pl_data$labels, 
                levels=names(sort(table(pl_data$labels[which(
                  pl_data$goodproba == paste(">", threshold))]), 
                  decreasing=TRUE)))
      }
      gg <- ggplot(pl_data, aes_string(x="labels", fill="goodproba")) +
        geom_bar() +
        scale_fill_brewer(direction=-1, palette="Accent",
                          name="Max\nconditional\nprobability") +
        scale_x_discrete(name="Cluster") +
        scale_y_continuous(name="Number of observations")
      print(gg)
    
    }
  }
  
  #####################################################
  ## PROBAPOST HISTOGRAM
  #####################################################
  if("probapost_histogram" %in% graphs) {
    pl_data <- data.frame(ID=ifelse(rep(length(rownames(y_profiles))==0, nrow(y_profiles)), 
                                    rep(1:nrow(y_profiles), times=ncol(y_profiles)),
                                    rownames(y_profiles)),
                          y_prof=as.vector(y_profiles), 
                          col_num=rep(1:ncol(y_profiles), each=nrow(y_profiles)),
                          col_nam=rep(colnames(y_profiles), each=nrow(y_profiles)),
                          conds=conds_vec,
                          labels=rep(labels, times=ncol(y_profiles)),
                          proba=rep(proba, times=ncol(y_profiles)))
    pl_data_tmp <- pl_data[which(pl_data$col_num == 1),]
    gg <- ggplot(pl_data_tmp, aes_string(x="proba")) +
      geom_histogram(binwidth = 0.01) +
      scale_x_continuous(name = "Maximum conditional probability") + theme_bw()
    print(gg)
  }

}


