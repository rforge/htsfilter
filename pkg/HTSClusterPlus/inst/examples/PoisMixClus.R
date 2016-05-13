## Simulate data as shown in Rau et al. (2011)
## Library size setting "A", high cluster separation
## n = 200 observations

simulate <- PoisMixSim(n = 200, libsize = "A", separation = "high")
y <- simulate$y
conds <- simulate$conditions

## Run the PMM model for K = 3
## "TC" library size estimate, EM algorithm
run <- PoisMixClus_K(y, K = 3, conds = conds, norm = "TC") 
summary(run)
## Estimates of pi and lambda for the selected model
pi.est <- run$pi
lambda.est <- run$lambda
norm <- run$norm

## Calculate the conditional probability of belonging to each cluster
proba <- probaPost(y, K = 3, conds = conds, pi = pi.est, s = norm,
                   lambda = lambda.est)

## Run the PMM model for K = 3 and 4
run <- PoisMixClus(y, K = 3:4, conds = conds, norm="TC")
summary(run)
