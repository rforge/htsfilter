
set.seed(12345)

## Simulate data as shown in Rau et al. (2011)
## Library size setting "A", high cluster separation
## n = 2000 observations

simulate <- PoisMixSim(n = 200, libsize = "A", separation = "high")
y <- simulate$y
conds <- simulate$conditions

## Run the PMM model for K = 3
## "TC" library size estimate, EM algorithm
run <- PoisMixClus_K(y, K=3, conds=conds, norm="TC")

## Estimates of pi and lambda for the selected model
pi.est <- run$pi
lambda.est <- run$lambda

