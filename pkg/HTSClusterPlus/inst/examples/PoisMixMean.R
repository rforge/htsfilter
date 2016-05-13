set.seed(12345)

## Simulate data as shown in Rau et al. (2011)
## Library size setting "A", high cluster separation
## n = 200 observations

simulate <- PoisMixSim(n = 200, libsize = "A", separation = "high")
y <- simulate$y
conds <- simulate$conditions
s <- colSums(y) / sum(y) 	## TC estimate of lib size

## Run the PMM-II model for K = 3
## "TC" library size estimate, EM algorithm

run <- PoisMixClus_K(y, K = 3, norm = "TC", conds = conds)
pi.est <- run$pi
lambda.est <- run$lambda

## Calculate the per-cluster mean for each observation
means <- PoisMixMean(y, K = 3, conds, s, lambda.est)

