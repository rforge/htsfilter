invnorm <-
function(pvalonesided, nrep, BHth = 0.05) 
{
    listres = vector("list", 4)
    nbstudies = length(pvalonesided)
    nbreptot = sum(nrep)
    weight = sqrt(nrep/nbreptot)
    qnormpval= do.call(cbind,lapply(pvalonesided, FUN=function(x) qnorm(1-x)))
    statc=apply(qnormpval,1, FUN=function(x) sum(weight*x,na.rm=TRUE))
    rpvalc = 1 - pnorm(statc)
    res = which(p.adjust(rpvalc, method = "BH") <= BHth)
    listres[[1]] = res
    listres[[2]] = statc
    listres[[3]] = rpvalc
    listres[[4]] = p.adjust(rpvalc, method = "BH")
    names(listres) = c("DEindices", "TestStatistic", "rawpval", "adjpval")
    return(listres)
}
