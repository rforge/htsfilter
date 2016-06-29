## Simulate toy data, n = 300 observations
countmat <- matrix(rnorm(300*4), nrow=300, ncol=4)
conds <- rep(c("A","B","C","D"), each=2)

## Run the Normal mixture model for K = 2,3
run <- NormMixClus(y=countmat, K=2:3, iter=5)

