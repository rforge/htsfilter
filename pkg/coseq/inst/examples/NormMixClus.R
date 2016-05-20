## Simulate data as shown in Rau et al. (2011)
## Library size setting "A", high cluster separation
## n = 200 observations

y <- matrix(rnorm(300*8), nrow=300, ncol=8)
conds <- rep(c("A","B","C","D"), each=2)

## Run the Normal mixture model for K = 3 and 4
run <- coseq(y, K = 3:4, norm="none")
