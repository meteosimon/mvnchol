#' Cholesky MVN
#'
#' BAMLSS Families for MVN with Cholesky Parameterization
#'
#' BAMLSS families that models a multivariate Normal (Gaussian)
#' distribution by (modified) Cholesky decomposition of the covariance
#' matrix.
#' 
#'
#' @param k integer. The dimension of the multivariate distribution.
#' @param type character. Choose standard Cholesky "chol" or modified
#'        Cholesky "modified".
#' @param ... not used.
#' @return a bamlss family.
#' @export
mvnchol_bamlss <- function(k, type = c("chol", "modified"), ...) {
	type <- match.arg(type)
	switch(type,
		chol = mvn_chol(k = k, ...),
		modified = mvn_modchol(k = k, ...)
	)
}



