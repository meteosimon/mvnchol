library(jmcm)            # for cattle data
library(bamlss)
library(bamlssMVN)

# load and reshape cattle data
data("cattle", package = "jmcm")
df <- reshape(cattle, direction = "wide", timevar = "day", v.names = "weight", sep = "_")
rownames(df) <- df$id
df$id <- NULL


f_mus <- paste0(paste0(names(df)[2:12], " ~ 0 + s(group, bs = 're')"))
f_lamdiags <- paste0(paste0("lamdiag", seq_len(11)), " ~ 0 + s(group, bs = 're')")
f_lambdas <- paste0(combn(seq_len(11), 2, 
			  function(x) paste0("lambda", x[1], x[2])), 
		    " ~ 0 + s(group, bs = 're')")


f <- lapply(c(f_mus, f_lamdiags, f_lambdas), FUN = as.formula)

b <- bamlss(f, family = mvn_chol(k = 11), data = df, 
	    sampler = FALSE, criterion = "BIC", optimizer = bfit)

p <- predict(b, type = "parameter")

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
      Linvt[[i]][j, j] <- p[[paste0("lamdiag", j)]][i]
      if (j < k) {
        for (l in (j+1):k) {
          Linvt[[i]][j, l] <- p[[paste0("lambda", j, l)]][i]
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





