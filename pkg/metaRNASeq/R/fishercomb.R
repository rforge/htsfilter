fishercomb <-
function(indpval,BHth = 0.05) 
{
	listres = vector("list", 4)
	logpval=do.call(cbind,lapply(indpval, log))
	statc=apply(logpval,1, FUN=function(x) -2*sum(x,na.rm=TRUE))
	notNA=apply(logpval,1,FUN=function(x) sum(!(is.na(x))))	
	rpvalc = 1 - pchisq(statc, df=(2*notNA))
	res = which(p.adjust(rpvalc, method = "BH") <= BHth)
	listres[[1]] = res
    	listres[[2]] = statc
	listres[[3]] = rpvalc
   	listres[[4]] = p.adjust(rpvalc, method = "BH")
    	names(listres) = c("DEindices", "TestStatistic", "rawpval", "adjpval")
    	return(listres)
}
