#' MVN Simulator
#'
#' Simulate from MVN
#'
#' This simulator generates synthtic data from a multivaraite normal
#' with zero mean and unit standard deviation.
#' Various parameterizations can be selected (such as "AR1", "AR2",
#' "AD1", "AD2"). 
#'
#' @param n integer. Sample size.
#' @param k integer. Dimension of mvn distribution.
#' @param type character. Type of parameterization (see details).
#' @param rho numeric. Vector of parameters.
#' @param ... not used.
#' @importFrom stats toeplitz
#' @importFrom mvtnorm rmvnorm
#' @export
sim_mvn <- function(n, k = 3, type = "AR1", rho = 0.9, ...) {

  cv <- switch(type,
    "AR1" = {
      stopifnot(length(rho) == 1)
      stats::toeplitz(rho^(0:(k-1)))
#    },
#    "AR2" = {
#      stopifnot(length(rho) == 2)
#      toeplitz(rho^(0:(k-1)))
#    },
#    "AD1" = {
#      stopifnot(length(rho) == 2)
#      toeplitz(rho^(0:(k-1)))
#    },
#    "AD2" = {
#      stopifnot(length(rho) == 2)
#      toeplitz(rho^(0:(k-1)))
    }
  )

  mvtnorm::rmvnorm(n, sigma = cv)
}

