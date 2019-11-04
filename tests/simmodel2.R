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
  chol_cv <- solve(chol(cv))  # lambdas come from L^-1 not L
  lamdiag[ii,] <- diag(chol_cv)
  lambda[ii,] <- chol_cv[upper.tri(chol_cv)]

  ## Check if cv is invertible 
  if ( !is.matrix(try(chol(cv))) ) l <- l + 1

  y[ii,] <- mvtnorm::rmvnorm(1, mu, cv)

  log_dens_ref[ii] <- mvtnorm::dmvnorm(y[ii,], mu, cv, log = TRUE)
}

## Data
d <- as.data.frame(y)
names(d) <- paste0("y", 0:2)
d$x0 <- x0

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

plot(log_dens, log_dens_ref)

## parse family w/o scores through bamlss.frame to
## get a family with numerical scores
library(bamlss)
fam <- mvn_chol(k = 3)
fam$score <- NULL
f <- list(y0 ~ s(x0), y1 ~ s(x0), y2 ~ s(x0),
	  lamdiag1 ~ s(x0), lamdiag2 ~ s(x0), lamdiag3 ~ s(x0),
	  lambda12 ~ s(x0), lambda13 ~ s(x0), lambda23 ~ s(x0))
bf <- bamlss.frame(f, family = fam, data = d)
fam2 <- bf$family

##
fam <- mvn_chol(k = 3)
smu1 <- fam$score$mu1(y, par)
smu2 <- fam2$score$mu1(y, par)
sld1 <- fam$score$lamdiag1(y, par)
sld2 <- fam2$score$lamdiag1(y, par)
sla1 <- fam$score$lambda12(y, par)
sla2 <- fam2$score$lambda12(y, par)
# these are all equal (i.e. numerical and analytical scores match)


## Model
f <- list(y0 ~ s(x0), y1 ~ s(x0), y2 ~ s(x0),
	  lamdiag1 ~ s(x0), lamdiag2 ~ s(x0), lamdiag3 ~ s(x0),
	  lambda12 ~ s(x0), lambda13 ~ s(x0), lambda23 ~ s(x0))

b <- bamlss(f, family = mvn_chol(k = 3), data = d, sampler = FALSE, optimizer = bfit, eps = 0.001)

fit <- predict(b, type = "parameter")

## Effects
x11()
par(mfrow = c(3, 3))

for (i in 1:3) {
  for (j in 1:3) {
    if (i == j) {
      plot(sort(x0), par[[paste0("lamdiag", i)]][order(x0)], type = "l", lwd = 2)
      lines(sort(x0), fit[[paste0("lamdiag", i)]][order(x0)], col = 2, lty = 2, lwd = 2)
    } else {
      plot(sort(x0), par[[paste0("lambda", min(i, j), max(i, j))]][order(x0)], type = "l", lwd = 2)
      lines(sort(x0), fit[[paste0("lambda", min(i, j), max(i, j))]][order(x0)], col = 2, lty = 2, lwd = 2)
    }
  }
}

## Demonstration of L_inv property

corr_data <- mvtnorm::rmvnorm(1000, mu, cv)
L_inv <- t(solve(chol(cv)))
uncorr_data <- t(L_inv %*% t(corr_data))

par(mfrow = c(2,3))
plot(corr_data[,1] ~ corr_data[,2])
plot(corr_data[,1] ~ corr_data[,3])
plot(corr_data[,2] ~ corr_data[,3])
plot(uncorr_data[,1] ~ uncorr_data[,2])
plot(uncorr_data[,1] ~ uncorr_data[,3])
plot(uncorr_data[,2] ~ uncorr_data[,3])



