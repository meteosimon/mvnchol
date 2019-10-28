library("mvtnorm")
library("bamlss")

n <- 2000
k <- 3
R <- matrix( c(1, .6, .36, .6, 1, .6, .36, .6, 1), 3, 3)
D <- diag(sqrt(c(4,2,1)))
y <- rmvnorm(n, mean=1:3, sigma= D %*% R %*% D)

parMV <- parAR <- list()
parMV$mu1 <- rep(1, n)
parMV$mu2 <- rep(2, n)
parMV$mu3 <- rep(3, n)
parMV$sigma1 <- rep(2, n)
parMV$sigma2 <- rep(sqrt(2), n)
parMV$sigma3 <- rep(1, n)

parAR <- parMV 

parAR$rho <- rep(.6, n)

parMV$rho12 <- rep(.6 ,n)
parMV$rho13 <- rep(.36, n)
parMV$rho23 <- rep(.6, n)

##
densAR <-  bamlss:::log_dmvnormAR1(y=y, par=parAR)
densMV <-  bamlss:::log_dmvnorm(y=y, par=parMV)
dens3 <- dmvnorm(y, mean=1:3, sigma= D %*% R %*% D, log=TRUE)

plot(densMV, densAR)
abline(c(0,1))

##
j <- 2
muAR <-  bamlss:::mu_score_mvnormAR1(y=y, par=parAR, j=j)
muMV <-  bamlss:::mu_score_mvnorm(y=y, par=parMV, j=j)

plot(muMV, muAR)
abline(c(0,1))

##
j <- 1
sigAR <-  bamlss:::sigma_score_mvnormAR1(y=y, par=parAR, j=j)
sigMV <-  bamlss:::sigma_score_mvnorm(y=y, par=parMV, j=j)

plot(sigMV, sigAR)
abline(c(0,1))

##
rC <-  bamlss:::rho_score_mvnormAR1(y=y, par=parAR)
rR <-  bamlss:::rho_score_mvnormAR1_R(y=y, par=parAR)

plot(rR, rC)
abline(c(0,1))

