library(jmcm)            # for cattle data
library(bamlss)
library(bamlssMVN)

# load and reshape cattle data
data("cattle", package = "jmcm")
df <- reshape(cattle, direction = "wide", timevar = "day", v.names = "weight")
rownames(df) <- df$id
df$id <- NULL


# source("../R/prototype_chol.R")

### C routines not working...
# lambda_score_mvnchol <- lambda_score_mvnchol_ref
# lamdiag_score_mvnchol <- lamdiag_score_mvnchol_ref
# mu_score_mvnchol <- mu_score_mvnchol_ref
# log_dmvnchol <- log_dmvnchol_ref



f_mus <- paste0(paste0(names(df)[2:12], " ~ 0 + group"))
f_lamdiags <- paste0(paste0("lamdiag", seq_len(11)), " ~ 0 + group")
f_lambdas <- paste0(combn(seq_len(11), 2, 
			  function(x) paste0("lambda", x[1], x[2])), 
		    " ~ 0 + group")


f <- lapply(c(f_mus, f_lamdiags, f_lambdas), FUN = as.formula)

b <- bamlss(f, family = mvn_chol(k = 11), data = df, 
	    sampler = FALSE, optimizer = bfit)

p <- predict(b, type = "parameter")

matricize <- function(p) {
  k <- length(grep("lamdiag", names(fit)))
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



