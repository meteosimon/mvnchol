library(jmcm)            # for cattle data
library(bamlss)
library(bamlssMVN)

# load and reshape cattle data
data("cattle", package = "jmcm")
df <- reshape(cattle, direction = "wide", timevar = "day", v.names = "weight", sep = "_")
rownames(df) <- df$id
df$id <- NULL

matricize <- function(p) {
  k <- length(grep("lamdiag", names(p)))
  n <- length(p[[1]])

  Linvt <- list()
  Tinvt <- list()
  D <- list()
  Siginv <- list()
  Sigma <- list()
  Omega <- list()

  for (i in 1:n) {
    Linvt[[i]] <- matrix(0, nrow = k, ncol = k)
    for (j in 1:k) {
      Linvt[[i]][j, j] <- as.numeric(p[[paste0("lamdiag", j)]][i])
      if (j < k) {
        for (l in (j+1):k) {
          Linvt[[i]][j, l] <- as.numeric(p[[paste0("lambda", j, l)]][i])
        }
      }
    }

  D[[i]] <- matrix(0, nrow = k, ncol = k)
  diag(D[[i]]) <- diag(Linvt[[i]])

  Tinvt[[i]] <- t(apply(Linvt[[i]], 1, function(x) x / diag(D[[i]])))

  Siginv[[i]] <- Linvt[[i]] %*% t(Linvt[[i]])

  Sigma[[i]] <- solve(Siginv[[i]])

  temp <- matrix(0, nrow = k, ncol = k)
  diag(temp) <- 1/sqrt(diag(Sigma[[i]]))

  Omega[[i]] <- temp %*% Sigma[[i]] %*% temp

  }



  return(list("L" = Linvt, "D" = D, "T" = Tinvt, "precision" = Siginv,
              "covariance" = Sigma, "correlation" = Omega))
}


## Helper function to extract offdiagonals of matrix as vector
offdiag <- function(x, k = 0) {x[which((col(x) - row(x)) == k)]}

## Reproduce Fig. 2 in jmcm package JSS paper
## (log innov variances vs. time) and (autoreg coefficients vs. lag)

# Use model estimated without penalization
b1 <- readRDS("cattle_model1.RDS")
p1 <- predict(b1, type = "parameter")

p1_matrices <- matricize(p1)
T_A1 <- p1_matrices$T[[1]]
loginnovars1 <- log(1/diag(p1_matrices$D[[1]])^2)

genautoregs1 <- matrix(NA, nrow = 10, ncol = 10)
for (i in 1:10) {
  genautoregs1[1:length(offdiag(T_A1, i)), i] <- -offdiag(T_A1, i)
}

# Use model with ridge penalization
b2 <- readRDS("cattle_model2.RDS")
p2 <- predict(b2, type = "parameter")

p2_matrices <- matricize(p2)
T_A2 <- p2_matrices$T[[1]]
loginnovars2 <- log(1/diag(p2_matrices$D[[1]])^2)

genautoregs2 <- matrix(NA, nrow = 10, ncol = 10)
for (i in 1:10) {
  genautoregs2[1:length(offdiag(T_A2, i)), i] <- -offdiag(T_A2, i)
}


xx <- matrix(1:10, ncol = 10, nrow = 10, byrow = TRUE)


x11()
par(mfrow = c(2, 2))
plot(loginnovars1, xlab = "Time", ylab = "Log-innovat. var.")
plot(xx, genautoregs1, xlab = "Lag", ylab = "Autoregres. coeffic.")
plot(loginnovars2, xlab = "Time", ylab = "Log-innovat. var.")
plot(xx, genautoregs2, xlab = "Lag", ylab = "Autoregres. coeffic.")


## Plot 95 percent confidence intervals for parameters

# b3 <- readRDS("cattle_model_mcmc.RDS")
b3 <- readRDS("cattle_model_mcmc_burn1000_thin5.RDS")
p95 <- predict(b3, newdata = df[1,], type = "parameter", FUN = c95) 
p95_matrices <- matricize(p95)

T_025 <- p95_matrices$T[[1]]
T_500 <- p95_matrices$T[[2]]
T_975 <- p95_matrices$T[[3]]

loginnovars025 <- log(1/diag(p95_matrices$D[[1]])^2)
loginnovars500 <- log(1/diag(p95_matrices$D[[2]])^2)
loginnovars975 <- log(1/diag(p95_matrices$D[[3]])^2)


genautoregs025 <- matrix(NA, nrow = 10, ncol = 10)
genautoregs500 <- matrix(NA, nrow = 10, ncol = 10)
genautoregs975 <- matrix(NA, nrow = 10, ncol = 10)

for (i in 1:10) {
  genautoregs025[1:length(offdiag(T_025, i)), i] <- -offdiag(T_025, i)
  genautoregs500[1:length(offdiag(T_500, i)), i] <- -offdiag(T_500, i)
  genautoregs975[1:length(offdiag(T_975, i)), i] <- -offdiag(T_975, i)
}


x11()
par(mfrow = c(2, 2))
plot(loginnovars500, xlab = "Time", ylab = "Log-innovat. var.", pch = 19, 
     ylim = c(1.5, 5))
arrows(x0 = 1:11, y0 = loginnovars025, x1 = 1:11, y1 = loginnovars975, 
       length = .05, angle = 90, code = 3)
plot(xx, genautoregs500, xlab = "Lag", ylab = "Autoregres. coeffic.", 
     ylim = c(-0.7, 2), pch = 19)
arrows(x0 = xx, y0 = genautoregs025, x1 = xx, y1 = genautoregs975,
       length = .05, angle = 90, code = 3)
abline(h = 0, lty = 2)
plot(loginnovars2, xlab = "Time", ylab = "Log-innovat. var.", 
     ylim = c(1.5, 5), pch = 19)
plot(xx, genautoregs2, xlab = "Lag", ylab = "Autoregres. coeffic.", 
     ylim = c(-0.7, 2), pch = 19)
abline(h = 0, lty = 2)


## Compare regularized and unregularized covariance and precision matrix
cor_reg <- p2_matrices$correlation[[1]]
cor_un <- p95_matrices$correlation[[2]]
cov_reg <- p2_matrices$covariance[[1]]
cov_un <- p95_matrices$covariance[[2]]
prec_reg <- p2_matrices$precision[[1]]
prec_un <- p95_matrices$precision[[2]]

x11()
par(mfrow = c(2, 2))
# image(cor_un[, nrow(cor_un):1])
# image(cor_reg[, nrow(cor_reg):1])
image(cov_un[, nrow(cov_un):1], xaxt = 'n', yaxt = 'n')
image(cov_reg[, nrow(cov_reg):1], xaxt ='n', yaxt = 'n')
image(prec_un[, nrow(prec_un):1], xaxt = 'n', yaxt = 'n')
image(prec_reg[, nrow(prec_reg):1], xaxt = 'n', yaxt = 'n')

