#' Cholesky MVN (disttree)
#'
#' disttree Families for MVN with Cholesky Parameterization
#'
#' NOTE: These functions are under development!! 
#' disttree families that models a multivariate Normal (Gaussian)
#' distribution by (modified) Cholesky decomposition of the covariance
#' matrix.
#'
#' @param k integer. The dimension of the multivariate distribution.
#' @param type character. Choose \code{"basic"} Cholesky decomposition or \code{"modified"}
#'        Cholesky decomposition. (For back compatibility \code{"chol"} is identical to \code{"basic"}.)
#' @param ... not used.
#' @return a bamlss family.
#' @export
dist_mvnchol <- function(k, type = c("basic", "modified", "chol"), ...) {
	type <- match.arg(type)
	switch(type,
		basic    = dist_mvn_chol(k = k, ...),
		chol     = dist_mvn_chol(k = k, ...),
		modified = dist_mvn_modchol(k = k, ...)
	)
}

