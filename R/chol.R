#' Cholesky MVN
#'
#' BAMLSS Family for MVN with Cholesky Parameterization
#'
#' This is a prototype implementation of a BAMLSS family that models
#' a multivariate Normal (Gaussian) distribution by a Cholesky
#' decomposition of the covariance matrix.
#'
#' @param k integer. The dimension of the multivariate distribution.
#' @param ... not used.
#' @return a bamlss family.
#' @export
#' @useDynLib mvnchol, .registration = TRUE
mvn_chol <- function(k = 2L, ...) {
  # --- set names of distributional parameters ---
  mu <- paste0("mu", seq_len(k))
  lamdiag <- paste0("lamdiag", seq_len(k))
  lambda <- utils::combn(seq_len(k), 2, function(x) paste0("lambda", x[1], x[2]))
  k_lambda <- k * (k-1) / 2    ## number of lambda parameters

  # --- set names of link functions ---
  links <- c(
    rep("identity", k),
    rep("log", k),
    rep("identity", k_lambda)
  )
  names(links) <- c(mu, lamdiag, lambda)

  # --- family list ---
  rval <- list(
    "family" = "mvnchol",
    "names"  = c(mu, lamdiag, lambda),   # set names of dist parameters
    "links"  = links,
    "d" = function(y, par, log = FALSE) {
      d <- log_dmvnchol(y, par)
      if(!log)
        d <- exp(d)
      return(d)
    } 
  )

  # --- add score funcitons ---
  mu_score_calls <- paste0(
    "function(y, par, ...) {mu_score_mvnchol(y, par, j = ", seq_len(k),")}"
  )
  lamdiag_score_calls <- paste0(
    "function(y, par, ...) {lamdiag_score_mvnchol(y, par, j = ", seq_len(k),")}"
  )
  lambda_score_calls <- utils::combn(seq_len(k), 2, function(x) {paste0(
    "function(y, par, ...) {lambda_score_mvnchol(y, par, i = ", x[1],", j = ", x[2],")}"
  )})

  scores <- list()
  for(j in seq_along(mu)) {
    scores[[mu[j]]] <- eval(parse(text = mu_score_calls[j]))
  }
  for(j in seq_along(lamdiag)) {
    scores[[lamdiag[j]]] <- eval(parse(text = lamdiag_score_calls[j]))
  }
  for(j in seq_along(lambda)) {
    scores[[lambda[j]]] <- eval(parse(text = lambda_score_calls[j]))
  }
  rval$score <- scores

  # --- add precision matrix function ---
  rval$precision <- function(par) {

    k <- length(grep("lamdiag", names(par)))
    n <- length(par[[1]])

    Linvt <- list()
    Siginv <- list()

    for (i in 1:n) {
      Linvt[[i]] <- matrix(0, nrow = k, ncol = k)
      for (j in 1:k) {
        Linvt[[i]][j, j] <- par[[paste0("lamdiag", j)]][i]
        if (j < k) {
          for (l in (j+1):k) {
            Linvt[[i]][j, l] <- par[[paste0("lambda", j, l)]][i]
          }
        }
      }
    Siginv[[i]] <- Linvt[[i]] %*% t(Linvt[[i]])
    }

    return(Siginv)
  }

  # --- add covariance matrix function ---
  rval$covariance <- function(par) {

    k <- length(grep("lamdiag", names(par)))
    n <- length(par[[1]])

    Linvt <- list()
    Sig <- list()

    for (i in 1:n) {
      Linvt[[i]] <- matrix(0, nrow = k, ncol = k)
      for (j in 1:k) {
        Linvt[[i]][j, j] <- par[[paste0("lamdiag", j)]][i]
        if (j < k) {
          for (l in (j+1):k) {
            Linvt[[i]][j, l] <- par[[paste0("lambda", j, l)]][i]
          }
        }
      }
    L <- solve(Linvt[[i]])
    Sig[[i]] <- t(L) %*% L
    }

    return(Sig)
  }

  # --- add correlation matrix function ---
  rval$correlation <- function(par) {

    k <- length(grep("lamdiag", names(par)))
    n <- length(par[[1]])

    Linvt <- list()
    Cor <- list()

    for (i in 1:n) {
      Linvt[[i]] <- matrix(0, nrow = k, ncol = k)
      for (j in 1:k) {
        Linvt[[i]][j, j] <- par[[paste0("lamdiag", j)]][i]
        if (j < k) {
          for (l in (j+1):k) {
            Linvt[[i]][j, l] <- par[[paste0("lambda", j, l)]][i]
          }
        }
      }
    L <- solve(Linvt[[i]])
    Sig <- t(L) %*% L
    inv_sds <- matrix(0, nrow = k, ncol = k)
    diag(inv_sds) <- 1 / sqrt(diag(Sig))
    Cor[[i]] <- inv_sds %*% Sig %*% inv_sds 
    }

    return(Cor)
  }



  # --- set class and return ---
  class(rval) <- "family.bamlss"
  rval
}

# #' @param y amatrix n x k.
# #' @param par a list with k+k+k(k-1)/2 elements named like the parameters of the family,
# #'            each element is a numeric vector of length n.
log_dmvnchol_ref <- function(y, par) {
  n <- nrow(y) # number of observations
  k <- ncol(y) # dimension of gaussian distribution
  mu <- paste0("mu", seq_len(k))
  lamdiag <- paste0("lamdiag", seq_len(k))
  lambda <- utils::combn(seq_len(k), 2, function(x) paste0("lambda", x[1], x[2]))
  par_full <- do.call("cbind", par)
  term_1 <- -k / 2 * log(2 * pi)
  term_2 <- apply(log(subset(par_full, select = lamdiag)), 1, sum)     
  term_3 <- vector(length = n) # initialise term_3 vector
  y_til <- as.matrix(y - subset(par_full, select = mu))
  for (ni in 1:n) {
    y_tild <- y_til[ni, ]
    L_inv <- matrix(0, nrow = k, ncol = k) # initialise L^-1 matrix
    for (l in 1:length(lambda)) { # assign off diagonal values
      i <- utils::combn(k, 2)[1, l]
      j <- utils::combn(k, 2)[2, l]
      L_inv[j, i] <- par[[paste0("lambda", i, j)]][ni] # j,i since lower trian
    }  
    for (i in 1:length(lamdiag)) {
      L_inv[i, i] <- par[[paste0("lamdiag", i)]][ni]
    }
    term_3[ni] <- -1 / 2 * norm(L_inv %*% y_tild, type = "2") ^ 2   
  }
  ll <- term_1 + term_2 + term_3
  return(ll)
}

log_dmvnchol_C <- function(y, par) {
  y <- as.matrix(y)
  storage.mode(y) <- "numeric"
  n <- nrow(y)
  k <- ncol(y)
  par <- do.call("cbind", par)
  .Call("log_dmvncholC", y, par, n, k, PACKAGE = "mvnchol")
}

# choose `log_dmvnchol_ref` or `log_dmvnchol_C` for computing the log-density
log_dmvnchol <- log_dmvnchol_C

# #' @param j dimension of parameter
mu_score_mvnchol_ref <- function(y, par, j) {
  n <- nrow(y) # number of observations
  k <- ncol(y) # dimension of gaussian distribution
  mu <- paste0("mu", seq_len(k))
  lamdiag <- paste0("lamdiag", seq_len(k))
  lambda <- utils::combn(seq_len(k), 2, function(x) paste0("lambda", x[1], x[2]))
  par_full <- do.call("cbind", par)
  y_til <- as.matrix(y - subset(par_full, select = mu))
  dl_dmu <- vector(length = n)
  for (ni in 1:n) {
    y_tild <- y_til[ni, ]
    L_inv <- matrix(0, nrow = k, ncol = k) # initialise L^-1 matrix
    temp <- utils::combn(k, 2)
    for (l in 1:length(lambda)) { # assign off diagonal values
      ii <- temp[1, l]
      jj <- temp[2, l]
      L_inv[jj, ii] <- par[[paste0("lambda", ii, jj)]][ni]
    }
    for (ii in 1:length(lamdiag)) {
      L_inv[ii, ii] <- par[[paste0("lamdiag", ii)]][ni]
    }
    Sigma_inv <- t(L_inv) %*% L_inv
    dl_dmu[ni] <- sum(Sigma_inv[j, ] * y_tild)
  }
  return(dl_dmu)
}      	

mu_score_mvnchol_C <- function(y, par, j) {
  y <- as.matrix(y)
  storage.mode(y) <- "numeric"
  n <- nrow(y)
  k <- ncol(y)
  par <- do.call("cbind", par)
  j <- as.integer(j)
  .Call("mu_score_mvncholC", y, par, n, k, j, PACKAGE = "mvnchol")
}

# choose `mu_score_mvnchol_ref` or `mu_score_mvnchol_C` for computing mu-scores
mu_score_mvnchol <- mu_score_mvnchol_C


lamdiag_score_mvnchol_ref <- function(y, par, j) {
  n <- nrow(y) # number of observations
  k <- ncol(y) # dimension of gaussian distribution
  mu <- paste0("mu", seq_len(k))
  lamdiag <- paste0("lamdiag", seq_len(k))
  lambda <- utils::combn(seq_len(k), 2, function(x) paste0("lambda", x[1], x[2]))
  par_full <- do.call("cbind", par)
  y_til <- as.matrix(y - subset(par_full, select = mu))
  dl_dlamdiag <- vector(length = n)
  for (ni in 1:n) {
    y_tild <- y_til[ni, ]
    L_inv <- matrix(0, nrow = k, ncol = k) # initialise L^-1 matrix
    temp <- utils::combn(k, 2)
    for (l in 1:length(lambda)) { # assign off diagonal values
      ii <- temp[1, l]
      jj <- temp[2, l]
      L_inv[jj, ii] <- par[[paste0("lambda", ii, jj)]][ni]
    }
    for (ii in 1:length(lamdiag)) {
      L_inv[ii, ii] <- par[[paste0("lamdiag", ii)]][ni]
    }
    dl_dlamdiag[ni] <- 1 - L_inv[j, j] * y_tild[j] * 
	               sum(y_tild[1:j] * L_inv[j, 1:j])
  }
  return(dl_dlamdiag)
}


lamdiag_score_mvnchol_C <- function(y, par, j) {
  y <- as.matrix(y)
  storage.mode(y) <- "numeric"
  n <- nrow(y)
  k <- ncol(y)
  par <- do.call("cbind", par)
  j <- as.integer(j)
  .Call("lamdiag_score_mvncholC", y, par, n, k, j, PACKAGE = "mvnchol")
}

# choose `mu_score_mvnchol_ref` or `mu_score_mvnchol_C` for computing mu-scores
lamdiag_score_mvnchol <- lamdiag_score_mvnchol_C



# #' @param i dimension of parameter
lambda_score_mvnchol_ref <- function(y, par, i, j) {
  n <- nrow(y) # number of observations
  k <- ncol(y) # dimension of gaussian distribution
  mu <- paste0("mu", seq_len(k))
  lamdiag <- paste0("lamdiag", seq_len(k))
  lambda <- utils::combn(seq_len(k), 2, function(x) paste0("lambda", x[1], x[2]))
  par_full <- do.call("cbind", par)
  y_til <- as.matrix(y - subset(par_full, select = mu))
  dl_dlambda <- vector(length = n)
  for (ni in 1:n) {
    y_tild <- y_til[ni, ]
    L_inv <- matrix(0, nrow = k, ncol = k) # initialise L^-1 matrix
    temp <- utils::combn(k, 2)
    for (l in 1:length(lambda)) { # assign off diagonal values
      ii <- temp[1, l]
      jj <- temp[2, l]
      L_inv[jj, ii] <- par[[paste0("lambda", ii, jj)]][ni]
    }
    for (ii in 1:length(lamdiag)) {
      L_inv[ii, ii] <- par[[paste0("lamdiag", ii)]][ni]
    }
    dl_dlambda[ni] <- - y_tild[i] *
                       sum(y_tild[1:j] * L_inv[j, 1:j])
  }
  return(dl_dlambda) 

}

lambda_score_mvnchol_C <- function(y, par, i, j) {
  y <- as.matrix(y)
  storage.mode(y) <- "numeric"
  n <- nrow(y)
  k <- ncol(y)
  par <- do.call("cbind", par)
  i <- as.integer(i)
  j <- as.integer(j)
  .Call("lambda_score_mvncholC", y, par, n, k, i, j, PACKAGE = "mvnchol")
}

# choose `mu_score_mvnchol_ref` or `mu_score_mvnchol_C` for computing mu-scores
lambda_score_mvnchol <- lambda_score_mvnchol_C




