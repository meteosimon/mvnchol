## Simulate test case

## build orthogonal rotation matrix
thetax <- pi/4
thetay <- pi/4
thetaz <- pi/4
Rx <- matrix( c(1,0,0, 0,cos(thetax),sin(thetax), 0,-sin(thetax),cos(thetax) ), 3, 3 )
Ry <- matrix( c(cos(thetay),0,-sin(thetay), 0,1,0,  sin(thetay),0,cos(thetay) ), 3, 3 )
Rz <- matrix( c(cos(thetaz),sin(thetaz),0, -sin(thetaz),cos(thetaz),0, 0,0,1 ), 3, 3 )
R <- Rx %*% Ry %*% Rz


f0 <- function(x) qnorm(x)
f1 <- function(x) (sin(pi * x))^2
f2 <- function(x) (cos(pi * x))^2

n <- 2000
set.seed(111)
x0 <- runif(n); x1 <- runif(n); x2 <- runif(n)

## initialize vectors for parameter lists
p12 <- NULL
p13 <- NULL
p23 <- NULL
sig <- matrix(0, n, 3)

lamdiag <- matrix(0, n, 3)
lambda <- matrix(0, n, 3)

y <- matrix(0, n, 3)
## eigenvalues
val1 <- f1(x0)
val2 <- f2(x0)
val3 <- rep(0, n)

log_dens_ref <- rep(0, n)

tau <- .1
l <- 0
dens1 <- NULL
for ( ii in seq(n) ) {

  mu <- rep(0, 3)
  
  val <- diag( c(val1[ii], val2[ii], val3[ii]) ) + diag(tau, 3)
  ## compute covariance matrix from rotation matrix and eigenvalues
  cv <- R %*% val %*% t(R)

  ## compute parameters for parameter list
  sig[ii,] <- sqrt(diag(cv))
  p12[ii] <- cv[1,2]
  p13[ii] <- cv[1,3]
  p23[ii] <- cv[2,3]

  ## compute paramters for Cholesky family
  chol_cv <- chol(cv)
  lamdiag[ii,] <- diag(chol_cv)
  lambda[ii,] <- chol_cv[upper.tri(chol_cv)]

  ## Check if cv is invertible 
  if ( !is.matrix(try(chol(cv))) ) l <- l + 1

  y[ii,] <- mvtnorm::rmvnorm(1, mu, cv)

  log_dens_ref[ii] <- mvtnorm::dmvnorm(y[ii,], mu, cv, log = TRUE)
}

## make parameter list for mvn chol family
par <- list()
par[["mu1"]] <- rep(0,n)
par[["mu2"]] <- rep(0,n)
par[["mu3"]] <- rep(0,n)
par[["lamdiag1"]] <- lamdiag[,1]
par[["lamdiag2"]] <- lamdiag[,2]
par[["lamdiag3"]] <- lamdiag[,3]
par[["lambda12"]] <- lambda[,1]
par[["lambda13"]] <- lambda[,2]
par[["lambda23"]] <- lambda[,3]

source("../R/prototype_chol.R")

dens_fun <- mvn_chol(k = 3)$d
log_dens <- dens_fun(y, par, log = TRUE)

par(mfrow = c(1, 2))
plot(log_dens, log_dens_ref)

## make parameter list for mvnorm family
par2 <- list()
par2[["mu1"]] <- rep(0,n)
par2[["mu2"]] <- rep(0,n)
par2[["mu3"]] <- rep(0,n)
par2[["sigma1"]] <- sig[,1]
par2[["sigma2"]] <- sig[,2]
par2[["sigma3"]] <- sig[,3]
par2[["rho12"]] <- p12
par2[["rho13"]] <- p13
par2[["rho23"]] <- p23

dens_fun2 <- bamlss::mvnorm_bamlss(k = 3)$d
log_dens2 <- dens_fun2(y, par2, log = TRUE)

plot(log_dens2, log_dens_ref)

