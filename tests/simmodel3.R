library("bamlss")
## Model

## build orthogonal rotation matrix
theta <- pi/4
R1 <- matrix( c(cos(theta),sin(theta),0,0, -sin(theta),cos(theta),0,0, 0,0,1,0, 0,0,0,1 ), 4, 4 )
R2 <- matrix( c(cos(theta),0,sin(theta),0, 0,1,0,0, -sin(theta),0,cos(theta),0, 0,0,0,1 ), 4, 4 )
R3 <- matrix( c(cos(theta),0,0,sin(theta), 0,1,0,0, 0,0,1,0, -sin(theta),0,0,cos(theta) ), 4, 4 )
R <- R1 %*% R2 %*% R3

f0 <- function(x) qnorm(x)
f1 <- function(x) (sin(pi * x))^2
f2 <- function(x) (cos(pi * x))^2

n <- 2000
set.seed(111)
x0 <- runif(n); x1 <- runif(n); x2 <- runif(n)

p12 <- NULL
p13 <- NULL
p14 <- NULL
p23 <- NULL
p24 <- NULL
p34 <- NULL
sig <- matrix(0, n, 4)

y <- matrix(0, n, 4)
## eigenvalues
val1 <- rep(0,n)
val2 <- f1(x0)
val3 <- f2(x1)
val4 <- rep(0,n)

tau <- .1
l <- 0
for ( ii in seq(n) ) {

  mu <- rep(0,4)
  
  val <- diag( c(val1[ii], val2[ii], val3[ii], val4[ii]) ) + diag(tau, 4)
  ## compute covariance matrix from rotation matrix and eigenvalues
  cv <- R %*% val %*% t(R)
  sig[ii,] <- sqrt(diag(cv))
  p12[ii] <- cv[1,2]
  p13[ii] <- cv[1,3]
  p14[ii] <- cv[1,4]
  p23[ii] <- cv[2,3]
  p24[ii] <- cv[2,4]
  p34[ii] <- cv[3,4]

  ## Check if cv is invertible 
  if ( !is.matrix(try(chol(cv))) ) l <- l + 1

  y[ii,] <- rmvn(1, mu, cv)

}

for( ii in 1:2 ) {
 ## plot input
 X11(); par(mfrow = c(3,4))
 if ( ii==1 )
   x <- x0
 if ( ii==2 )
   x <- x1
 plot(sig[,1] ~ x); plot(sig[,2] ~ x); plot(sig[,3] ~ x); plot(sig[,4] ~ x)
 plot((p12/ (sig[,1]*sig[,2])) ~ x)
 plot((p13/ (sig[,1]*sig[,3])) ~ x)
 plot((p14/ (sig[,1]*sig[,4])) ~ x)
 plot((p23/ (sig[,2]*sig[,3])) ~ x)
 plot((p24/ (sig[,2]*sig[,4])) ~ x)
 plot((p34/ (sig[,3]*sig[,4])) ~ x)
}

dat <- data.frame(y0=y[,1],y1=y[,2],y2=y[,3],y3=y[,4],x0=x0,x1=x1)

###################################################################
#library("bamlss")
#load("../data/simmodel3.rda")

f <- list(
  y0 ~ 1,
  y1 ~ 1,
  y2 ~ 1,
  y3 ~ 1,
  sigma1 ~ s(x0) + s(x1),
  sigma2 ~ s(x0) + s(x1),
  sigma3 ~ s(x0) + s(x1),
  sigma4 ~ s(x0) + s(x1),
  rho12 ~ s(x0) + s(x1),
  rho13 ~ s(x0) + s(x1),
  rho14 ~ s(x0) + s(x1),
  rho23 ~ s(x0) + s(x1),
  rho24 ~ s(x0) + s(x1),
  rho34 ~ s(x0) + s(x1)
)

f2 <- list(
  y0 ~ 1,
  y1 ~ 1,
  y2 ~ 1,
  y3 ~ 1,
  sigma1 ~ s(x0) + s(x1),
  sigma2 ~ s(x0) + s(x1),
  sigma3 ~ s(x1),
  sigma4 ~ 1,
  rho12 ~ s(x0) + s(x1),
  rho13 ~ s(x0) + s(x1),
  rho14 ~ 1,
  rho23 ~ s(x0) + s(x1),
  rho24 ~ 1,
  rho34 ~ 1
)

X11()
b <- bamlss(f2, family = gF("mvnorm", k = 4), data = dat,
            optimizer = boost, sampler = FALSE,
            maxit = 500, nu = 0.1, scale.d = TRUE)

X11()
plot(b, pages=1, ask=FALSE)

## I would also like to try a formula like this:
f <- list(
  y0 ~ 1,
  y1 ~ 1,
  y2 ~ 1,
  y3 ~ 1,
  sigma1 ~ id1 + id2,
  sigma2 ~ id1 + id2,
  id1 ~ s(x0),
  id2 ~ s(x1),
  sigma3 ~ s(x1),
  sigma4 ~ 1,
  rho12 ~ s(x0) + s(x1),
  rho13 ~ id3,
  rho14 ~ 1,
  rho23 ~ id3,
  id3 ~ s(x1),
  rho24 ~ 1,
  rho34 ~ 1
)

###################################################################


s1 <- predict(b, model="sigma1", type="parameter")
s2 <- predict(b, model="sigma2", type="parameter")
s3 <- predict(b, model="sigma3", type="parameter")
s4 <- predict(b, model="sigma4", type="parameter")
r12 <- predict(b, model="rho12", type="parameter")
r13 <- predict(b, model="rho13", type="parameter")
r14 <- predict(b, model="rho14", type="parameter")
r23 <- predict(b, model="rho23", type="parameter")
r24 <- predict(b, model="rho24", type="parameter")
r34 <- predict(b, model="rho34", type="parameter")

## plot output
#X11(); par(mfrow=c(3,4))
#plot(s1~x0); plot(s2~x0); plot(s3~x0); plot(s4~x0)
#plot(r12~x0); plot(r13~x0); plot(r14~x0)
#plot(r23~x0); plot(r24~x0); plot(r34~x0)

## plot input vs. fitted values
X11(); par(mfrow=c(3,4))
plot(sig[,1] ~ s1); abline(c(0,1))
plot(sig[,2] ~ s2); abline(c(0,1))
plot(sig[,3] ~ s3); abline(c(0,1)) 
plot(sig[,4] ~ s4); abline(c(0,1))
plot((p12/ (sig[,1]*sig[,2])) ~ r12); abline(c(0,1))
plot((p13/ (sig[,1]*sig[,3])) ~ r13); abline(c(0,1))
plot((p14/ (sig[,1]*sig[,4])) ~ r14); abline(c(0,1))
plot((p23/ (sig[,2]*sig[,3])) ~ r23); abline(c(0,1))
plot((p24/ (sig[,2]*sig[,4])) ~ r24); abline(c(0,1))
plot((p34/ (sig[,3]*sig[,4])) ~ r34); abline(c(0,1))


