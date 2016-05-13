set.seed(12345)

## Simulate data as shown in Rau et al. (2011)
## Library size setting "A", low cluster separation
## n = 200 observations

simulate <- PoisMixSim(n = 200, libsize = "A", separation = "low")
y <- simulate$y
conds <- simulate$conditions
w <- rowSums(y)               ## Estimate of w
r <- table(conds)             ## Number of replicates per condition
d <- length(unique(conds))    ## Number of conditions
s <- colSums(y) / sum(y)      ## TC estimate of lib size
s.dot <- rep(NA, d)           ## Summing lib size within conditions
for(j in 1:d) s.dot[j] <- sum(s[which(conds == unique(conds)[j])]);

## Initial guess for pi and lambda
g.true <- 4
pi.guess <- simulate$pi
## Recalibrate so that (s.dot * lambda.guess) = 1
lambda.sim <- simulate$lambda
lambda.guess <- matrix(NA, nrow = d, ncol = g.true)
for(k in 1:g.true) {
    tmp <- lambda.sim[,k]/sum(lambda.sim[,k])
    lambda.guess[,k] <- tmp/s.dot
}

## Run the PMM model for K = 4
## with EM algorithm and "TC" library size parameter
run <- PoisMixClus_K(y, K = 4, norm = "TC", conds = conds)
pi.est <- run$pi
lambda.est <- run$lambda

## Mean values for each of the parameter sets
mean.guess <- PoisMixMean(y, 4, conds, s, lambda.guess)
mean.est <- PoisMixMean(y, 4, conds, s, lambda.est)

## Difference in log likelihoods
LL.guess <- logLikePoisMix(y, mean.guess, pi.guess)
LL.est <- logLikePoisMix(y, mean.est, pi.est)
## Difference in log likelihoods       
LL.diff <- logLikePoisMixDiff(y, mean.guess, pi.guess, mean.est, pi.est)
