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
#' @param type character. Choose \code{"basic"} Cholesky decomposition or \code{"modified"}
#'        Cholesky decomposition.
#' @param ... not used.
#' @return a bamlss family.
#' @export
mvnchol_bamlss <- function(k, type = c("basic", "modified", "chol"), ...) {
	type <- match.arg(type)
	switch(type,
		basic = mvn_chol(k = k, ...),
		chol = mvn_chol(k = k, ...),
		modified = mvn_modchol(k = k, ...)
	)
}



