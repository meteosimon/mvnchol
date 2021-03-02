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
#' \begin{itemize}
#'   \item \code{init}: Time of initialization of the NWP model.
#'   \item \code{obs_*}: Observations for lead time \code{*}.
#'   \item \code{mean_ens_*}: NWP ensemble mean for lead time \code{*}.
#'   \item \code{logsd_ens_*}: NWP logarithm of ensemble standard deviation for lead time \code{*}.
#'   \item \code{yday}: Yearday.
#' \end{itemize}
#'
#' @usage data(TempIbk)
#'
#' @examples
#' \dontrun{
#' library(bamlss)
#' library(mvnchol)
#' load("TempIbk.rda")
#' 
#' mean_names <- names(TempIbk)[grep("mean_", names(TempIbk))]
#' logsd_names <- names(TempIbk)[grep("logsd_", names(TempIbk))]
#' obs_names <- names(TempIbk)[grep("obs", names(TempIbk))]
#' 
#' f_mus <- paste0(obs_names, " ~ s(yday, bs = 'cc') + s(yday, bs = 'cc', by = ", mean_names, ")")
#' f_mus[[1]]
#' f_lamdiags <- paste0("lamdiag", seq_along(logsd_names),
#' 	" ~ s(yday, bs = 'cc') + s(yday, bs = 'cc', by = ", logsd_names, ")")
#' f_lamdiags[[1]]
#' 
#' paste_lambda <- function(x) paste0("lambda", x[1], x[2], " ~ s(yday, bs = 'cc')")
#' f_lambdas <- combn(seq_along(logsd_names), 2, paste_lambda, simplify = FALSE)
#' f_lambdas[[1]]
#' 
#' f <- lapply(c(f_mus, f_lamdiags, f_lambdas), as.formula)
#' 
#' fam <- mvnchol_bamlss(k = length(logsd_names), type = "basic")
#' 
#' b <- bamlss(f, family = fam, data = TempIbk,
#' 	sampler = FALSE, optimizer = opt_boost, maxit = 1000)
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

