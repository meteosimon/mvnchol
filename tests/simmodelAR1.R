library("bamlss")
## Model

## simulate data
k <- 5
r <- .6

f1 <- function(x) (sin(pi * x))^2
f2 <- function(x) (cos(pi * x))^2

n <- 2000
set.seed(111)
x0 <- runif(n); x1 <- runif(n); x2 <- runif(n)

y <- matrix(0, n, k)
## rho
r <- f1(x0) + f2(x1)
r <- f1(x0)
r <- 0.4 + 0.4 * r/max(r)

for ( ii in seq(n) ) {

  mu <- rep(0,k)
  sig <- matrix(0, ncol=k, nrow=k)
  for ( i in seq(k) ) {
    for ( j in seq(k) ) {
      sig[i,j] <- r[ii]^abs(i-j)
    }
  }
  
  y[ii,] <- rmvn(1, mu, sig)

}

dat <- data.frame(y1=y[,1],y2=y[,2],y3=y[,3],y4=y[,4],y5=y[,5],x0=x0,x1=x1)

###################################################################

f <- list(
  y1 ~ 1, y2 ~ 1, y3 ~ 1, y4 ~ 1, y5 ~ 1,
  sigma1 ~ 1, sigma2 ~ 1, sigma3 ~ 1, sigma4 ~ 1, sigma5 ~ 1,
  rho12 ~ s(x0)+s(x1),
  rho13 ~ s(x0)+s(x1),
  rho14 ~ s(x0)+s(x1),
  rho15 ~ s(x0)+s(x1),
  rho23 ~ s(x0)+s(x1),
  rho24 ~ s(x0)+s(x1),
  rho25 ~ s(x0)+s(x1),
  rho34 ~ s(x0)+s(x1),
  rho35 ~ s(x0)+s(x1),
  rho45 ~ s(x0)+s(x1)
)

f <- list(
  y1 ~ 1, y2 ~ 1, y3 ~ 1, y4 ~ 1, y5 ~ 1,
  sigma1 ~ 1, sigma2 ~ 1, sigma3 ~ 1, sigma4 ~ 1, sigma5 ~ 1,
  rho12 ~ s(x0),
  rho13 ~ s(x0),
  rho14 ~ s(x0),
  rho15 ~ s(x0),
  rho23 ~ s(x0),
  rho24 ~ s(x0),
  rho25 ~ s(x0),
  rho34 ~ s(x0),
  rho35 ~ s(x0),
  rho45 ~ s(x0)
)

b <- bamlss(f, family = gF("mvnorm", k = 5), data = dat,
            optimizer = boost, sampler = FALSE,
            maxit = 2000, nu = 0.1, scale.d = TRUE)

fAR <- list(
  y1 ~ 1, y2 ~ 1, y3 ~ 1, y4 ~ 1, y5 ~ 1,
  sigma1 ~ 1, sigma2 ~ 1, sigma3 ~ 1, sigma4 ~ 1, sigma5 ~ 1,
  rho ~ s(x0)
)

bAR <- bamlss(fAR, family = gF("mvnormAR1", k = 5), data = dat,
            optimizer = boost, sampler = FALSE,
            maxit = 100, nu = 0.1, scale.d = TRUE)

X11(); plot(b, ask=FALSE, pages=1)
###################################################################

r12 <- predict(b, model="rho12", type="parameter")
r13 <- predict(b, model="rho13", type="parameter")
r14 <- predict(b, model="rho14", type="parameter")
r15 <- predict(b, model="rho15", type="parameter")
r23 <- predict(b, model="rho23", type="parameter")
r24 <- predict(b, model="rho24", type="parameter")
r25 <- predict(b, model="rho25", type="parameter")
r34 <- predict(b, model="rho34", type="parameter")
r35 <- predict(b, model="rho35", type="parameter")
r45 <- predict(b, model="rho45", type="parameter")
rAR <- predict(bAR, model="rho", type="parameter")

par(mfrow=c(2,2))
plot(r~x0)
points(rAR~x0, col=2)
points(r12~x0, col=3)
points(r23~x0, col=4)
points(r34~x0, col=5)
points(r45~x0, col=6)

plot(r^2~x0)
points(r13~x0, col=2)
points(r24~x0, col=3)
points(r35~x0, col=4)

plot(r^3~x0)
points(r14~x0, col=2)
points(r25~x0, col=3)

plot(r^4~x0)
points(r15~x0, col=2)
