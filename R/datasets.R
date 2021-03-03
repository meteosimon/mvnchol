#' Temperature data.
#'
#' Temperature Data for Innsbruck Airport
#'
#' Numerical weather predictions (NWP) and observations of
#' 2 meter temperature at Innsbruck Airport.
#' The observations from the SYNOP station 11120 cover 5 years from
#  2015-01-01 to 2019-31-12.
#' The NWP data are derived from GEFS reforecasts (Hamill et al. 2013).
#' The data contain following variables:
#'
#' \itemize{
#'   \item \code{init}: Time of initialization of the NWP model.
#'   \item \code{obs_*}: Observations for lead time \code{*}.
#'   \item \code{mean_ens_*}: NWP ensemble mean for lead time \code{*}.
#'   \item \code{logsd_ens_*}: NWP logarithm of ensemble standard deviation for lead time \code{*}.
#'   \item \code{yday}: Yearday.
#' }
#'
#' @references
#'   Hamill TM, Bates GT, Whitaker JS, Murray DR, Fiorino M, Galarneau Jr TJ,
#'     Zhu Y, Lapenta W (2013). NOAA's Second-Generation Global Medium-Range
#'     Ensemble Reforecast Data Set. \emph{Bulletin of the American Meteorological
#'     Society}, 94(10), 1553-1565.
#'
#' @usage data(TempIbk)
#'
#' @examples
#' \dontrun{
#' ## packages and data
#' library("bamlss")
#' library("mvnchol")
#' data("TempIbk", package = "mvnchol")
#' 
#' ## five lead times
#' lead <- seq(192, 216, by = 6)
#' 
#' ## set up formulas
#' f <- c(
#'   ## mu equations
#'   sprintf('obs_%s ~ s(yday, bs = "cc") + s(yday, bs = "cc", by = mean_ens_%s)', lead, lead),
#' 
#'   ## lambda diag equations
#'   sprintf('lamdiag%s ~ s(yday, bs = "cc") + s(yday, bs = "cc", by = logsd_ens_%s)', 1:5, lead),
#' 
#'   ## lambda off-diag equations
#'   sprintf('lambda%s ~ s(yday, bs = "cc")', apply(combn(1:5, 2), 2, paste, collapse = ""))
#' )
#' f <- lapply(f, as.formula)
#' 
#' ## multivariate normal family with basic Cholesky parameterization
#' fam <- mvnchol_bamlss(k = 5, type = "basic")
#' 
#' ## fit boosting model
#' b <- bamlss(f, family = fam, data = TempIbk, optimizer = opt_boost, maxit = 1000, sampler = FALSE)
#' 
#' ## predict sample case
#' nd <- subset(TempIbk, init == "2015-10-10")
#' fit <- predict(b, newdata = nd, type = "parameter", mstop = 5000)
#' 
#' ## plot correlation matrix for GEFS initialization 2015-10-10
#' image(lead, lead, fam$correlation(fit)[[1]][5:1, ],
#' 	 col = hcl.colors(12, "Blues 3", rev = TRUE), axes = FALSE,
#' 	 xlab = "lead time in hours", ylab = "lead time in hours",
#' 	 main = "Correlation matrix fitted for 2015-10-10")
#' axis(1, lead)
#' axis(2, lead)
#' box()
#' 
#' ## plot means and standard deviations
#' stdev <- sqrt(diag(fam$covariance(fit)[[1]]))
#' means <- fam$means(fit)[[1]]
#' lower <- means - stdev
#' upper <- means + stdev
#' 
#' plot(lead, means, type = 'b', cex = 2, lwd = 1, lty = 2, axes = FALSE,
#' 	 ylim = c(min(lower), max(upper)),
#' 	 ylab = expression("Temperature in " * degree * "C"),
#' 	 xlab = "lead time in hours",
#' 	 main = "Means +/- one st. dev. for 2015-10-10")
#' segments(lead, y0 = lower, y1 = upper)
#' axis(1, lead)
#' axis(2)
#' box()
#' }
#'
"TempIbk"

#' Reference data.
#'
#' Simulated data to test the implementation of the bamlss families.
#'
#' @usage data(simdata)
#'
#' @examples
#' \dontrun{
#' ## Reproducing code
#' set.seed(111)
#' n <- 2000
#'
#' ## build orthogonal rotation matrix
#' thetax <- pi/4
#' thetay <- pi/4
#' thetaz <- pi/4
#' Rx <- matrix( c(1,0,0, 0,cos(thetax),sin(thetax), 0,-sin(thetax),cos(thetax) ), 3, 3 )
#' Ry <- matrix( c(cos(thetay),0,-sin(thetay), 0,1,0,  sin(thetay),0,cos(thetay) ), 3, 3 )
#' Rz <- matrix( c(cos(thetaz),sin(thetaz),0, -sin(thetaz),cos(thetaz),0, 0,0,1 ), 3, 3 )
#' R <- Rx %*% Ry %*% Rz
#' 
#' ## non-linear functions
#' f1 <- function(x) (sin(pi * x))^2
#' f2 <- function(x) (cos(pi * x))^2
#'
#' ## random derivitates  
#' x <- runif(n)
#' 
#' ## eigenvalues
#' val1 <- f1(x) 
#' val2 <- f2(x)
#' val3 <- rep(0, n)
#' 
#' ## initialize vectors for parameter lists
#' p12 <- NULL
#' p13 <- NULL
#' p23 <- NULL
#' sig <- matrix(0, n, 3)
#' 
#' lamdiag <- matrix(0, n, 3)
#' lambda <- matrix(0, n, 3)
#' 
#' y <- matrix(0, n, 3)
#' log_dens_ref <- rep(0, n)
#' 
#' tau <- .1  ## offset on diagonal
#' l <- 0     ## count occasions with invertible cv 
#' dens1 <- NULL
#' for ( ii in seq(n) ) {
#'     mu <- rep(0, 3)
#'   
#'     val <- diag( c(val1[ii], val2[ii], val3[ii]) ) + diag(tau, 3)
#'     ## compute covariance matrix from rotation matrix and eigenvalues
#'     cv <- R %*% val %*% t(R)
#'   
#'     ## compute parameters for parameter list
#'     sig[ii,] <- sqrt(diag(cv))
#'     p12[ii] <- cv[1,2]
#'     p13[ii] <- cv[1,3]
#'     p23[ii] <- cv[2,3]
#'   
#'     ## compute paramters for Cholesky family
#'     chol_cv <- solve(chol(cv))  # lambdas come from L^-1 not L
#'     lamdiag[ii,] <- diag(chol_cv)
#'     lambda[ii,] <- chol_cv[upper.tri(chol_cv)]
#'   
#'     ## Check if cv is invertible 
#'     if ( !is.matrix(try(chol(cv))) ) l <- l + 1
#'   
#'     y[ii,] <- mvtnorm::rmvnorm(1, mu, cv)
#'   
#'     log_dens_ref[ii] <- mvtnorm::dmvnorm(y[ii,], mu, cv, log = TRUE)
#' }
#' print(l)
#' 
#' ## Data
#' d <- as.data.frame(y)
#' names(d) <- paste0("y", 1:3)
#' d$x <- x
#' 
#' ## make parameter list for mvn chol family
#' par <- list()
#' par[["mu1"]] <- rep(0,n)
#' par[["mu2"]] <- rep(0,n)
#' par[["mu3"]] <- rep(0,n)
#' par[["lamdiag1"]] <- lamdiag[,1]
#' par[["lamdiag2"]] <- lamdiag[,2]
#' par[["lamdiag3"]] <- lamdiag[,3]
#' par[["lambda12"]] <- lambda[,1]
#' par[["lambda13"]] <- lambda[,2]
#' par[["lambda23"]] <- lambda[,3]
#' 
#' simdata <- list(
#'     d   = d,
#'     par = par,
#'     y   = y
#' )
#' 
#' save(simdata, file = "simdata.rda")
#' ## End of simulation
#' }
"simdata"

