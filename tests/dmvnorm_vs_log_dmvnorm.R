library("bamlss")
library("mvtnorm")
## Model

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

p12 <- NULL
p13 <- NULL
p23 <- NULL
sig <- matrix(0, n, 3)

y <- matrix(0, n, 3)
## eigenvalues
val1 <- f1(x0)
val2 <- f2(x0)
val3 <- rep(0,n)

tau <- .1
l <- 0
dens1 <- NULL
for ( ii in seq(n) ) {

  mu <- rep(0,3)
  
  val <- diag( c(val1[ii], val2[ii], val3[ii]) ) + diag(tau, 3)
  ## compute covariance matrix from rotation matrix and eigenvalues
  cv <- R %*% val %*% t(R)
  sig[ii,] <- sqrt(diag(cv))
  p12[ii] <- cv[1,2]
  p13[ii] <- cv[1,3]
  p23[ii] <- cv[2,3]

  ## Check if cv is invertible 
  if ( !is.matrix(try(chol(cv))) ) l <- l + 1

  y[ii,] <- rmvn(1, mu, cv)

  ## compare densities
  dens1[ii] <- dmvnorm(x=y[ii,], mean = mu, sigma = cv, log = TRUE)
}

par <- list()
## make parameter list
par[["mu1"]] <- rep(0,n)
par[["mu2"]] <- rep(0,n)
par[["mu3"]] <- rep(0,n)
par[["sigma1"]] <- sig[,1]
par[["sigma2"]] <- sig[,2]
par[["sigma3"]] <- sig[,3]
par[["rho12"]] <- p12
par[["rho13"]] <- p13
par[["rho23"]] <- p23


#X11()
par(mfrow=c(2,2))
## compare densities
dens2 <- bamlss:::log_dmvnorm(y=y, par=par)
plot(exp(dens1), exp(dens2), xlab="mvtnorm", ylab="bamlss")
abline(c(0,1))
plot(dens1, dens2, xlab="mvtnorm", ylab="bamlss")
abline(c(0,1))

## compare mu_score
for ( jj in seq(3) ) {
#  mu_scoreR <- bamlss:::mu_score_mvnormR(y=y, par=par, j=jj)
#  mu_scoreC <- bamlss:::mu_score_mvnorm(y=y, par=par, j=jj)
#  plot(mu_scoreR, mu_scoreC, main=jj)
#  abline(c(0,1))
}

## compare sigma_score
for ( jj in seq(3) ) {
#  sigma_scoreR <- bamlss:::sigma_score_mvnormR(y=y, par=par, j=jj)
#  sigma_scoreC <- bamlss:::sigma_score_mvnorm(y=y, par=par, j=jj)
#  plot(sigma_scoreR, sigma_scoreC, main=jj)
#  abline(c(0,1))
}

## compare rho_score
for ( ii in seq(3) ) {
 for ( jj in seq(3) ) {
  if ( ii < jj) {
#   rho_scoreR <- bamlss:::rho_score_mvnormR(y=y, par=par, i=ii, j=jj)
#   rho_scoreC <- bamlss:::rho_score_mvnorm(y=y, par=par, i=ii, j=jj)
#   plot(rho_scoreR, rho_scoreC, main=paste("i=",ii,"-- j=",jj))
#   abline(c(0,1))
  }
 }
}


