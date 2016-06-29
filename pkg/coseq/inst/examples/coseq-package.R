## Simulate toy data, n = 300 observations
countmat <- matrix(rnorm(300*8), nrow=300, ncol=8)
conds <- rep(c("A","B","C","D"), each=2)

## Run the Normal mixture model for K = 2,3,4
run <- coseq(y=countmat, K=3:4, norm="none")
