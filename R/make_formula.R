#' Formula Generator
#'
#' Generate a formula for a MVN Cholesky model.
#'
#' This is a helper to generate a formula for a bamlss model with \code{k}-dimensional
#' multi-variate normal distribution and Cholesky decomposed variance-covariance matrix.
#' It is helpful if one formula should be used for means, another for all diagonal
#' entries of the Cholesky factor, and a third one for all lower triangular entries
#' of the Cholesky factor.
#' The left hand side has \code{k} elements separated by \code{|}.
#' The right hand side has one to three elements separated by \code{|} specifying
#' the formulas used for all means, diagonal entries of the Cholesky factor and
#' lower triangular entries of the Cholesky factor, respectively.
#'
#' @param formula formula. 
#' @param type character. Type of Cholesky decomposition.
#' @examples
#' f <- O | C | E | A | N ~ s(x1) + s(x2) | s(y) | z 
#' f2 <- make_formula(f)
#' f2
#' @export
make_formula <- function(formula, type = "basic") {
	FORM <- Formula::as.Formula(formula)
	l <- length(FORM)
	k <- l[1]
	seq_k <- seq_len(k)

	if (l[2] == 1) {
		FORM <- stats::update(FORM, . ~ . | 1 | 1)
	}
	if (l[2] == 2) {
		FORM <- stats::update(FORM, . ~ . | . | 1)
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


