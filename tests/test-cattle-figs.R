library(jmcm)            # for cattle data
library(bamlss)
library(bamlssMVN)
library(latex2exp)
library(fields)

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

## Adjust jmcm regressogram function to allow for different mfrow
regressogram2 <- function(object, time, op2){
    debug <- 0
    op <- op2
    opt <- object@opt
    lambda <- opt$lambda
    gamma <- opt$gamma
    llmd <- length(lambda)
    lgma <- length(gamma)
    args <- object@args
    dims <- object@devcomp$dims
    m <- args[["m"]]
    Y <- args[["Y"]]
    X <- args[["X"]]
    Z <- args[["Z"]]
    W <- args[["W"]]
    if (length(unique(m)) != 1)
        stop("No regressograms. Unbalanced longitudinal dataset.")
    DataMat <- t(Y[1:m[1]])
    for (i in 2:length(m)) {
        DataMat <- rbind(DataMat, t(Y[(sum(m[1:(i - 1)]) + 1):sum(m[1:i])]))
    }
    dimnames(DataMat) <- NULL
    S <- cov(DataMat)
    R <- cor(DataMat)
    C <- t(chol(S))
    D <- diag(diag(C))
    Tt <- t(forwardsolve(C %*% diag(diag(C)^(-1)), diag(dim(D)[1])))
    Lt <- t(diag(diag(C)^(-1)) %*% C)
    ts <- seq(min(time), max(time), length.out = 100)
    tlag <- 1
    for (i in 2:10) {tlag <- c(tlag, i:1)}
    tslag <- seq(min(tlag), max(tlag), length.out = 100)
    Z.ts <- NULL
    W.tslag <- NULL
    for (i in 0:(llmd - 1)) Z.ts <- cbind(Z.ts, ts^i)
    for (i in 0:(lgma - 1)) W.tslag <- cbind(W.tslag, tslag^i)
    Zlmd <- Z.ts %*% lambda
    Wgma <- W.tslag %*% gamma
    if (dims["MCD"] == 1) {
        plot(time, log(diag(D)^2), xlab = "", ylab = "Log-innovat. var.", 
	     pch = 19, ylim = c(1.5, 5))
        lines(ts, Zlmd)
        phi <- -Tt[upper.tri(Tt, diag = FALSE)]
        plot(tlag, phi, xlab = "", ylab = "Autoregres. coeffic.", 
	     pch = 19, ylim = c(-.7, 2))
        lines(tslag, Wgma)
    }
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


cattleA <- subset(cattle, group == "A")
fit1 <- jmcm(weight | id | I(day / 14 + 1) ~ 1 | 1, data = cattleA, 
	     triple = c(8, 3, 4), cov.method = "mcd")

x11()
#par(mfrow = c(2, 2))
regressogram2(fit1, time = 1:11, par(mfrow = c(3,2), 
				     mar = c(4,4,2,2)))
plot(loginnovars1, xlab = "", ylab = "Log-innovat. var.")
plot(xx, genautoregs1, xlab = "", ylab = "Autoregres. coeffic.")
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


pdf("/home/thomas/Documents/Projects/PhD-Concept/figures/cattle_comp.pdf")
regressogram2(fit1, time = 1:11, par(mfrow = c(3,2),
                                     mar = c(4,4,0.5,2)))
plot(loginnovars500, xlab = "", ylab = "Log-innovqt. var.", pch = 19, 
     ylim = c(1.5, 5))
arrows(x0 = 1:11, y0 = loginnovars025, x1 = 1:11, y1 = loginnovars975, 
       length = .05, angle = 90, code = 3)
plot(xx, genautoregs500, xlab = "", ylab = "Autoregres. coeffic.", 
     ylim = c(-0.7, 2), pch = 19)
arrows(x0 = xx, y0 = genautoregs025, x1 = xx, y1 = genautoregs975,
       length = .05, angle = 90, code = 3)
abline(h = 0, lty = 2)
plot(loginnovars2, xlab = "Time", ylab = "Log-innovat. var.", 
     ylim = c(1.5, 5), pch = 19)
plot(xx, genautoregs2, xlab = "Lag", ylab = "Autoregres. coeffic.", 
     ylim = c(-0.7, 2), pch = 19)
abline(h = 0, lty = 2)
dev.off()

## Compare regularized and unregularized covariance and precision matrix
cor_reg <- p2_matrices$correlation[[1]]
cor_un <- p95_matrices$correlation[[2]]
cov_reg <- p2_matrices$covariance[[1]]
cov_un <- p95_matrices$covariance[[2]]
prec_reg <- p2_matrices$precision[[1]]
prec_un <- p95_matrices$precision[[2]]

# make values close to 0 NA so that they appear white in image plot
prec_reg[abs(prec_reg) < .00001] <- NA
prec_un[abs(prec_un) < .00001] <- NA


coolwarm_hcl <- colorspace::diverging_hcl(11,
  h = c(250, 10), c = 100, l = c(37, 88), power = c(0.7, 1.7))

brbg_hcl <- colorspace::diverging_hcl(101, h = c(180, 50), c = 80, 
				      l = c(20, 95), power = c(0.7, 1.3))

ylrd <- hcl.colors(10, "YlOrRd", rev = TRUE)

mx <- .33
brks <- seq(-mx, mx, length.out = 12)

pdf("/home/thomas/Documents/Projects/PhD-Concept/figures/cattle_sparse.pdf", 
    width = 7, height = 7)
par(mfrow = c(2, 2), mar = c(2,1,3,1.7), oma = c(0,0,0,4))
image(cov_un[, nrow(cov_un):1], xaxt = 'n', yaxt = 'n', 
      main = TeX('$\\Sigma$'), zlim = c(0, 500), 
      col = ylrd, asp = 1)
image(cov_reg[, nrow(cov_reg):1], xaxt ='n', yaxt = 'n',
      main = TeX('$\\Sigma_{reg}$'), zlim = c(0, 500), 
      col = ylrd, asp = 1)
image.plot(cov_reg, zlim = c(0, 500), 
           legend.only = TRUE, col = hcl.colors(10, "YlOrRd", rev = TRUE), 
	   smallplot= c(.95,1,.1,.85), breaks = 50*(0:10), 
	   axis.args = list(at = 50*(0:10)))
image(prec_un[, nrow(prec_un):1], xaxt = 'n', yaxt = 'n',
      main = TeX('$\\Sigma^{-1}$'), zlim = c(-mx, mx), 
      col = coolwarm_hcl, breaks = brks, asp = 1)
image(prec_reg[, nrow(prec_reg):1], xaxt = 'n', yaxt = 'n',
      main = TeX('$\\Sigma^{-1}_{reg}$'), zlim = c(-mx, mx),
      col = coolwarm_hcl, breaks = brks, asp = 1)
image.plot(zlim = c(-mx, mx), 
           legend.only = TRUE, col = coolwarm_hcl,
           breaks = brks, at = brks, axis.args = list(at = brks),
	   smallplot= c(.95,1,.1,.85))
dev.off()
