syntheticData <- function(H0number, plot = FALSE, plot.name = NA) {
  dat <- data(initialDataset, package = "HTSDiff")
  dat <- get(dat)
                                        #  dat <- initialDataset
  fix <- grep(".DE", dat[,1])
  modify <- dat[-fix, 1]
  synth <- dat[,1:5] ## BF/F data
  ## Choose H0number genes among those we can modify
  if(H0number < 1 & H0number >= 0) {
    subset <- sample(modify, H0number*length(modify))
  }
  if(H0number >= 1) {
    subset <- sample(modify, H0number)
  }
  tochange <- is.element(synth[,1], subset)
  synth[tochange, -1] <- dat[tochange, 4:7]
  synth[,1] <- as.character(synth[,1])
  synth[tochange, 1] <- paste(synth[tochange, 1], ".Hnull", sep="")
  if(plot == TRUE) {
    if(is.na(plot.name) == FALSE) pdf(plot.name);
    x<-apply(synth[,2:3],1,mean)
    y<-apply(synth[,4:5],1,mean)
    
    plot(log2(x),log2(y),pch=46, xlab = expression(paste(log[2](x))), ylab = expression(paste(log[2](y))))   
    points(log2(x)[grep(".Hnull",synth[,1])],log2(y)[grep(".Hnull",synth[,1])],pch=20,col="cyan",cex=0.75)
    points(log2(x)[grep("DE",synth[,1])],log2(y)[grep("DE",synth[,1])],pch=20,col=2)
    points(log2(x)[grep("NDE",synth[,1])],log2(y)[grep("NDE",synth[,1])],pch=20,col="blue")
    legend("topleft", bty="n", col=c("cyan", 2, "blue"), pch = c(20,20,20), c(expression(paste(H[0])), "DE", "NDE"),
           cex=1.5, pt.cex=c(0.75, 1, 1))
    if(is.na(plot.name) == FALSE) dev.off();
  }
  rownames(synth) <- synth[,1]
  return(synth[,-1])
}
