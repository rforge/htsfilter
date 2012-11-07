## matrix
setMethod(
	f= "HTSFilter",
	signature = signature(x="matrix"),
	definition = function(x, conds, s.min=1, s.max=200, s.len=100, 
		loess.span=0.3, normalization=c("TMM", "DESeq", "none"), 
		plot=TRUE, plot.name=NA) 

	{
		normalization <- match.arg(normalization)
		data <- x

		## Run filter
		filter <- .HTSFilterBackground(data=data, conds=conds, s.min=s.min,
			s.max=s.max, s.len=s.len, loess.span=loess.span,
			normalization=normalization, plot=plot, plot.name=plot.name)

		## Return various results
		filter.results <- list(filteredData =  filter$data.filter,
			on = filter$on, s = filter$s.optimal,
			indexValues = filter$index.values, normFactor = filter$norm.factor,
			removedData = data[which(filter$on == 0),])

		return(filter.results)
	}
)

## data.frame
setMethod(
	f= "HTSFilter",
	signature = signature(x="data.frame"),
	definition = function(x, conds, s.min=1, s.max=200, s.len=100, 
		loess.span=0.3, normalization=c("TMM", "DESeq", "none"), 
		plot=TRUE, plot.name=NA) 

	{
		normalization <- match.arg(normalization)
		data <- as.matrix(x)

		## Run filter
		filter <- .HTSFilterBackground(data=data, conds=conds, s.min=s.min,
			s.max=s.max, s.len=s.len, loess.span=loess.span,
			normalization=normalization, plot=plot, plot.name=plot.name)

		## Return various results
		filter.results <- list(filteredData =  filter$data.filter,
			on = filter$on, s = filter$s.optimal,
			indexValues = filter$index.values, normFactor = filter$norm.factor,
			removedData = data[which(filter$on == 0),])

		return(filter.results)
	}
)

## DGElist
setMethod(
	f="HTSFilter",
	signature = signature(x="DGEList"),
	definition = function(x, s.min=1, s.max=200, s.len=100,
		loess.span=0.3, normalization=c("TMM","DESeq","none"),
		plot=TRUE, plot.name=NA)
	{

		normalization <- match.arg(normalization)
		data <- x$counts
		conds <- x$samples$group
		

		## Run filter
		filter <- .HTSFilterBackground(data=data, conds=conds, s.min=s.min,
			s.max=s.max, s.len=s.len, loess.span=loess.span,
			normalization=normalization, plot=plot, plot.name=plot.name)

		## Create a new DGEList
		filteredData <- DGEList(counts = filter$data.filter,
			lib.size=x$lib.size, norm.factors=x$norm.factors,
			group=conds)

		## Return various results
		filter.results <- list(filteredData = filteredData,
			on = filter$on, s = filter$s.optimal,
			indexValues = filter$index.values, normFactor = filter$norm.factor,
			removedData = data[which(filter$on == 0),])

		return(filter.results)
	}
)


## CountDataSet
## setMethod(
## 	f="HTSFilter",
## 	signature = signature(x="CountDataSet"),
##	definition = function(x, conds=NA, s.min=1, s.max=200, s.len=100,
## 		loess.span=0.3, normalization=c("TMM","DESeq","none"),
## 		plot=TRUE, plot.name=NA)
## 	{
## 
## 		normalization <- match.arg(normalization)
## 		data <- counts(x)
## 		## TO DO: WHAT TO DO IF MORE THAN ONE FACTOR
## 		if(is.na(conds)[1] == TRUE) conds <- pData(x);
## 
## 		## Run filter
## 		filter <- .HTSFilterBackground(data=data, conds=conds, s.min=s.min,
## 			s.max=s.max, s.len=s.len, loess.span=loess.span,
## 			normalization=normalization, plot=plot, plot.name=plot.name)
## 
## 		## Assign filtered data to new count data set
## 		filteredData <- newCountDataSet(filter$data.filter, conditions = conds)		
## 		pData(filteredData) <- pData(x)	
## 
## 		## Return various results
## 		filter.results <- list(filteredData = filteredData,
## 			on = filter$on, s = filter$s.optimal,
## 			indexValues = filter$index.values, normFactor = filter$norm.factor,
## 			removedData = data[which(filter$on == 0),])
## 	
## 		return(filter.results)
## 	}
## )


## SeqExpressionSet
## setMethod(
## 	f="HTSFilter",
## 	signature = signature(x="SeqExpressionSet"),
## 	definition = function(x, conds=NA, s.min=1, s.max=200, s.len=100,
## 		loess.span=0.3, normalization=c("TMM","DESeq","none"),
## 		plot=TRUE, plot.name=NA)
## 	{
##		normalization <- match.arg(normalization)
## 		data <- exprs(x)
## 		## TO DO: WHAT TO DO IF MORE THAN ONE FACTOR
## 		if(is.na(conds)[1] == TRUE) conds <- pData(x);
## 
## 		## Run filter
## 		filter <- .HTSFilterBackground(data=data, conds=conds, s.min=s.min,
## 			s.max=s.max, s.len=s.len, loess.span=loess.span,
## 			normalization=normalization, plot=plot, plot.name=plot.name)
## 
## 		## Assign filtered data to matrix
## 		filteredData <- newSeqExpressionSet(filter$data.filter,
## 			 phenoData = phenoData(x))
## 
## 		## Return various results
## 		filter.results <- list(filteredData = filteredData,
##			on = filter$on, s = filter$s.optimal,
## 			indexValues = filter$index.values, normFactor = filter$norm.factor,
## 			removedData = data[which(filter$on == 0),])
## 
## 		return(filter.results)
## 	}
## )



