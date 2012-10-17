.perConditionSimilarityIndex <-
function(data.norm.perCondition, log.s){

	countsNormLog <- log(data.norm.perCondition)
	countsNormLog_binary <- countsNormLog > log.s
	sum <- 0
	nbindiv <- ncol(countsNormLog)
	calc <- mapply(function(i,j) {
		x <- countsNormLog_binary[,i]
		y <- countsNormLog_binary[,j]
		a <- length(which(x == 1 & y == 1))
		b <- length(which(x == 1 & y == 0))
		c <- length(which(x == 0 & y == 1))
		d <- length(which(x == 0 & y == 0))
		## Sanity check
		if(sum(c(a,b,c,d)) != length(x) & sum(c(a,b,c,d)) != length(y)) {
			stop("Something's wrong.")
		}
		a / (a+b+c)
		}, 
		rep(1:(nbindiv), times = c((nbindiv-1):0)),
		unlist(lapply(2:nbindiv, function(x) seq(x,nbindiv))) )
	return(sum(calc))
}


.normalizeData <-
function(data, conds, normalization) {
	if(normalization == "TMM") {
		N <- colSums(data)
		f <- calcNormFactors(data,method="TMM")
		TMM <- N*f / mean(N*f)
		norm.factor <- TMM
		data.norm <- scale(data, center=FALSE, scale=TMM)
	}
	if(normalization == "DESeq") {
		deseq <- estimateSizeFactorsForMatrix(data, locfunc = median)
		norm.factor <- deseq
		data.norm <- scale(data, center=FALSE, scale=deseq)
	}
	if(normalization == "none") {
		data.norm <- data
		norm.factor <- NA
	}
	return(list(data.norm = data.norm, norm.factor = norm.factor))
}


.HTSFilterBackground <- 
function(data, conds, s.min, s.max, s.len, 
	loess.span, normalization,  plot, plot.name) {

	## Sanity checks: replicated data?
	if(min(table(conds) == 1)) {
		stop("All conditions must have at least TWO replicates.")
	}
	if(s.min <= 0) stop("s.min must be > 0.");
	if(s.min > s.max) stop(paste(sQuote(s.min), "must be <", sQuote(s.max)));
	if(dim(data)[2] != length(conds)) {
		stop(paste(sQuote(conds), " must be of length ", dim(data)[2], 
		" (the number of columns of", sQuote(data), ")", sep = ""))
	}

	## Normalization (calculated on full dataset)
	norm <- .normalizeData(data = data, conds = conds, 
		normalization = normalization)
	data.norm <- norm$data.norm
	norm.factor <- norm$norm.factor

	## Calculate index for each threshold value
	s.test <- seq(log(s.min), log(s.max), length = s.len)
	index <- rowSums(matrix(mapply(function(condition, s) {
		.perConditionSimilarityIndex(data.norm.perCondition = data.norm[,which(conds == unique(conds)[condition])], 
		log.s = s)}, 
		rep(1:length(unique(conds)), times = length(s.test)), 
		rep(s.test, each = length(unique(conds)))), ncol = length(unique(conds)), 
		byrow = TRUE))

	## Choosing s.optimal: Fit a loess curve to data
	index[which(is.nan(index) == TRUE)] <- NA
	p <- predict(loess(index ~ exp(s.test), span = loess.span), 
		seq(s.min, s.max, by = 0.001))
	s.optimal <- seq(s.min, s.max, by = 0.001)[which(p == max(p, na.rm = TRUE))]
	if(plot == TRUE) {
		if(is.na(plot.name) != TRUE) pdf(plot.name, width = 6, height = 6);
		plot(exp(s.test), index, log = "x", xlab = "Threshold", 
			ylab = "Similarity index")
		lines(seq(s.min, s.max, by = 0.001), p, col = "blue", lwd = 2)
		points(s.optimal, max(p, na.rm = TRUE), col = "red", pch = "X", cex= 1.5)
		abline(v = s.optimal, lty = 2, col= "red")
		legend("bottomleft", pch = "X", col="red", cex = 1.25, paste("s =", s.optimal), bty = "n")
		if(is.na(plot.name) != TRUE) dev.off();
	}

	## Calculate which genes are present in at least one sample with s.optimal
	data.binary <- matrix(NA, nrow = length(conds), ncol = dim(data)[1])
	data.binary <- data.norm > s.optimal
	on <- ifelse(rowSums(data.binary) == 0, 0, 1)

	data.filter <- data[which(on == 1),]
	colnames(data.binary) <- colnames(data)
	rownames(data.binary) <- rownames(data)

	filter.results <- list(on = on, s.optimal = s.optimal, data.binary = data.binary,
		index.values = data.frame(threshold = exp(s.test), index = index), 
		data.filter = data.filter, norm.factor = norm.factor)

	return(filter.results)
}