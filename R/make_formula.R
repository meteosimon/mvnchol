#' Formula Generator
#'
#' @param formula, e.g. O | C | E | A | N ~ x | y | z
#' @param character. Type of Cholesky decomposition.
#' @export
make_formula <- function(formula, type = "chol") {
	FORM <- Formula::as.Formula(formula)
	l <- length(FORM)
	k <- l[1]
	seq_k <- seq_len(k)

	if (l[2] == 1) {
		FORM <- update(FORM, . ~ . | 1 | 1)
	if (l[2] == 2) {
		FORM <- update(FORM, . ~ . | . | 1)
	}
	
	fam <- mvnchol_bamlss(k = k, type = type)
	nms <- fam$names

	rval <- c(
		lapply(seq_k, formula, x = FORM, rhs = 1),
		rep(list(formula(FORM, lhs = FALSE, rhs = 2)), k),
		rep(list(formula(FORM, lhs = FALSE, rhs = 3)), k*(k-1)/2)
	)
	names(rval) <- nms
	rval
}


