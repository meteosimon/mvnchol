library("bamlss")
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

## plot input
X11(); par(mfrow = c(2,3))
plot(sig[,1] ~ x0); plot(sig[,2] ~ x0); plot(sig[,3] ~ x0)
plot((p12/ (sig[,1]*sig[,2])) ~ x0)
plot((p13/ (sig[,1]*sig[,3])) ~ x0)
plot((p23/ (sig[,2]*sig[,3])) ~ x0)

dat <- data.frame(y1=y[,1],y2=y[,2],y3=y[,3],x=x0)

###################################################################
library("bamlss")
load("../data/simmodel2.rda")

f <- list(
  y0 ~ 1,
  y1 ~ 1,
  y2 ~ 1,
  sigma1 ~ s(x0),
  sigma2 ~ s(x0),
  sigma3 ~ s(x0),
  rho12 ~ s(x0),
  rho13 ~ s(x0),
  rho23 ~ s(x0)
)

b <- bamlss(f, family = gF("mvnorm", k = 3), data = dat,
            optimizer = boost, sampler = FALSE,
            maxit = 500, nu = 0.1, scale.d = TRUE)

b <- bamlss(f, family = mvnorm_bamlss(k = 3), data = dat,
            optimizer = boost, sampler = FALSE,
            maxit = 500, nu = 0.1, scale.d = TRUE)


X11(); plot(b, ask=FALSE, pages=1)
###################################################################


s1 <- predict(b, model="sigma1", type="parameter")
s2 <- predict(b, model="sigma2", type="parameter")
s3 <- predict(b, model="sigma3", type="parameter")
r12 <- predict(b, model="rho12", type="parameter")
r13 <- predict(b, model="rho13", type="parameter")
r23 <- predict(b, model="rho23", type="parameter")

## plot output
X11(); par(mfrow=c(2,3))
plot(s1~x0); plot(s2~x0); plot(s3~x0)
plot(r12~x0); plot(r13~x0); plot(r23~x0)

###################################################################
## try to get this parallel
library("parallel")

dat$chunk <- rep(1:4, each=500)
models <- list()
for ( j in seq(4) ) {
  models[[j]] <- list("fo"=f, "dtrain"=subset(dat, chunk!=j),
                              "dverif"=subset(dat, chunk==j))
}

blist <- function(obj, ...) {
            environment(obj$fo) <- environment()
            b <- bamlss(formula = obj$fo,
                     family = gF("mvnorm", k = 3),
                     data = obj$dtrain,
                     optimizer = boost, sampler = FALSE,
                     maxit = 100, nu = 0.1, scale.d = TRUE)
            return(list("bobj"=b, "test"=obj$dverif))
}

predict.blist <- function(obj) {
            fit <- predict(obj$bobj, newdata=obj$test)
            return(fit)
}

## Works:
r <- lapply(models, blist)
fit <- lapply(r, predict.blist)
fit.cv <- NULL
for ( j in seq(4) ) {
  fit.cv <- rbind(fit.cv, do.call("cbind", fit[[j]]))
}
## Does not work:
## bamlss as 
res <- mclapply(models, blist, mc.cores=4L,
         mc.silent=TRUE, mc.allow.recursive=FALSE)
###################################################################

