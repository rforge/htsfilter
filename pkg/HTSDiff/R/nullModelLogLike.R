nullModelLogLike <- function(counts, conds, norm="DESeq")
{
	counts <- as.matrix(counts)
	conds <- as.vector(conds)
	if(length(unique(conds)) != 2) {
		stop("The number of unique conditions must be 2.")
	}
   	if (norm == "TC") 
            s <- colSums(counts)/sum(counts)
      if(norm == "UQ") 
            s <- apply(counts, 2, quantile, 0.75)/sum(apply(counts, 2, 
                quantile, 0.75))
      if(norm == "Med") 
            s <- apply(counts, 2, median)/sum(apply(counts, 2, median))
      if(norm == "DESeq") {
            loggeomeans <- rowMeans(log(counts))
            s <- apply(counts, 2, function(x) exp(median((log(x) - 
                loggeomeans)[is.finite(loggeomeans)])))
            s <- s/sum(s)
      }
      if(norm == "TMM") {
            f <- calcNormFactors(as.matrix(counts), method = "TMM")
            s <- colSums(counts) * f/sum(colSums(counts) * f)
      }

	n <- nrow(counts)
	mean.calc <- PoisMixMean(y=counts, g=1, conds=conds, s=s, 
		lambda=matrix(1, ncol=1, nrow=2))
      LL.tmp <- logLikePoisMix(y=counts, mean.calc, pi=1)
      LL <- LL.tmp$ll
      t <- probaPost(y=counts, g=1, conds, pi=1, s=s, lambda=matrix(1,ncol=1,nrow=2))
      map <- unlist(apply(t, 1, function(x) which(x == max(x, 
            na.rm = TRUE))[1]))
      z <- matrix(1, nrow = n, ncol = 1)
 	np <- n   ## Just the w's to be estimated    
   
      BIC <- LL - (np/2) * log(n)
	return(list(logLike=LL, BIC=BIC, ICL=BIC))
}