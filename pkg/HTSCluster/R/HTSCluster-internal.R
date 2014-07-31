.myloopfxn <- function(k, lambda, w.mat, s.mat, r, n, cols) {
	lambda.mat <- matrix(rep(rep(lambda[,k], times = r), each = n),
		nrow = n, ncol = cols)
	return(w.mat * s.mat * lambda.mat)
}
.myfxn <- function(var1, var2) {
	rowSums(ifelse(var1 == 0, var2, var1 * log(var1/var2) + var2 - var1))
}