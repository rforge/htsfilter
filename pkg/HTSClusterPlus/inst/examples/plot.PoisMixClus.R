set.seed(12345)

## Simulate data as shown in Rau et al. (2011)
## Library size setting "A", high cluster separation
## n = 2000 observations
simulate <- PoisMixSim(n = 200, libsize = "A", separation = "high")
y <- simulate$y
conds <- simulate$conditions

## Run the PMM-II model for g = 3
## "TC" library size estimate, EM algorithm
run <- PoisMixClus(y, K=c(2:3), 
                   norm = "TC", conds = conds)

## Visualization of results (not run):
## plot(run)
