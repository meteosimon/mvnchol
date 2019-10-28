library("bamlss")
## Model.
f7 <- function(x) {
  eta <- .5 + 2 * cos(pi*x)
  eta / sqrt(1 + eta^2)
}

n <- 2000
set.seed(111)
x0 <- runif(n)

rho12 <- f7(x0)
rho13 <- rep(0, n)
rho23 <- rep(0, n)

y <- matrix(0, n, 3)

l <- 0
for(i in 1:n) {
  V <- diag(c(1,1,1))
  V[1, 2] <- rho12[i]
  V[1, 3] <- rho13[i]
  V[2, 3] <- rho23[i]
  V[2, 1] <- V[1, 2]
  V[3, 1] <- V[1, 3]
  V[3, 2] <- V[2, 3]
  mu <- c(0,0,0)

  if ( !is.matrix(try(chol(V))) ) l <- l + 1

  y[i,] <- rmvn(1, mu, V)
}

dat <- data.frame(y0=y[,1],y1=y[,2],y2=y[,3],x0=x0)

f <- list(
  y0 ~ 1,
  y1 ~ 1,
  y2 ~ 1,
  sigma1 ~ 1,
  sigma2 ~ 1,
  sigma3 ~ 1,
  rho12 ~ s(x0),
  rho13 ~ 1,
  rho23 ~ 1
)

b <- bamlss(f, family = gF("mvnorm", k = 3), data = dat, optimizer = boost, sampler = FALSE, maxit = 100, nu = 0.1, scale.d = TRUE)

p <- predict(b, model="rho12", type="parameter")
par(mfrow=c(1,2))
plot(p~x0)
plot(f7(x0)~x0)

