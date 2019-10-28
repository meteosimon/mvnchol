library("bamlss")
load("simmodel2.rda")

f <- list(
  y0 ~ 1,
  y1 ~ 1,
  y2 ~ 1,
  sigma1 ~ 1,
  sigma2 ~ s(x0),
  sigma3 ~ s(x0),
  rho12 ~ s(x0),
  rho13 ~ s(x0),
  rho23 ~ s(x0)
)

b <- bamlss(f, family = gF("mvnorm", k = 3), data = dat,
            optimizer = boost, sampler = FALSE,
            maxit = 500, nu = 0.1, scale.d = TRUE)

X11(); plot(b, ask=FALSE, pages=1)

