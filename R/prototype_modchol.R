#' Modified Cholesky MVN
#'
#' BAMLSS Family for MVN with Modified Cholesky Parameterization
#'
#' This is a prototype implementation of a BAMLSS family that models
#' a multivariate Normal (Gaussian) distribution by a modified Cholesky
#' decomposition of the covariance matrix.
#'
#' @param k integer. The dimension of the multivariate distribution.
#' @param ... not used.
#' @return a bamlss family.
#' @export
#' @useDynLib bamlssMVN, .registration = TRUE
mvn_modchol <- function(k = 2L, ...) {
  # --- set names of distributional parameters ---
  mu <- paste0("mu", seq_len(k))
  innov <- paste0("innov", seq_len(k))
  phi <- utils::combn(seq_len(k), 2, function(x) paste0("phi", x[1], x[2]))
  k_phi <- k * (k-1) / 2    ## number of phi parameters

  # --- set names of link functions ---
  links <- c(
    rep("identity", k),
    rep("log", k),
    rep("identity", k_phi)
  )
  names(links) <- c(mu, innov, phi)

  # --- family list ---
  rval <- list(
    "family" = "mvnmodchol",
    "names"  = c(mu, innov, phi),   # set names of dist parameters
    "links"  = links,
    "d" = function(y, par, log = FALSE) {
      d <- log_dmvnmodchol(y, par)
      if(!log)
        d <- exp(d)
      return(d)
    } 
  )

  # --- add score funcitons ---
  mu_score_calls <- paste0(
    "function(y, par, ...) {mu_score_mvnmodchol(y, par, j = ", seq_len(k),")}"
  )
  innov_score_calls <- paste0(
    "function(y, par, ...) {innov_score_mvnmodchol(y, par, j = ", seq_len(k),")}"
  )
  phi_score_calls <- utils::combn(seq_len(k), 2, function(x) {paste0(
    "function(y, par, ...) {phi_score_mvnmodchol(y, par, i = ", x[1],", j = ", x[2],")}"
  )})

  scores <- list()
  for(j in seq_along(mu)) {
    scores[[mu[j]]] <- eval(parse(text = mu_score_calls[j]))
  }
  for(j in seq_along(innov)) {
    scores[[innov[j]]] <- eval(parse(text = innov_score_calls[j]))
  }
  for(j in seq_along(phi)) {
    scores[[phi[j]]] <- eval(parse(text = phi_score_calls[j]))
  }
  rval$score <- scores

  # --- add hess funcitons ---
  mu_hess_calls <- paste0(
    "function(y, par, ...) {mu_hess_mvnmodchol(y, par, j = ", seq_len(k),")}"
  )
  innov_hess_calls <- paste0(
    "function(y, par, ...) {innov_hess_mvnmodchol(y, par, j = ", seq_len(k),")}"
  )
  phi_hess_calls <- utils::combn(seq_len(k), 2, function(x) {paste0(
    "function(y, par, ...) {phi_hess_mvnmodchol(y, par, i = ", x[1],", j = ", x[2],")}"
  )})

  hesses <- list()
  for(j in seq_along(mu)) {
    hesses[[mu[j]]] <- eval(parse(text = mu_hess_calls[j]))
  }
  for(j in seq_along(innov)) {
    hesses[[innov[j]]] <- eval(parse(text = innov_hess_calls[j]))
  }
  for(j in seq_along(phi)) {
    hesses[[phi[j]]] <- eval(parse(text = phi_hess_calls[j]))
  }
  rval$hess <- hesses

  # --- add precision matrix function ---
  rval$precision <- function(par) {

    k <- length(grep("innov", names(par)))
    n <- length(par[[1]])

    Linvt <- list()
    Siginv <- list()

    for (i in 1:n) {
      Linvt[[i]] <- matrix(0, nrow = k, ncol = k)
      for (j in 1:k) {
        Linvt[[i]][j, j] <- par[[paste0("innov", j)]][i]
        if (j < k) {
          for (l in (j+1):k) {
            Linvt[[i]][j, l] <- -par[[paste0("phi", j, l)]][i] *
		    par[[paste0("innov", l)]][i]
          }
        }
      }
    Siginv[[i]] <- Linvt[[i]] %*% t(Linvt[[i]])
    }

    return(Siginv)
  }

  # --- add covariance matrix function ---
  rval$covariance <- function(par) {

    k <- length(grep("innov", names(par)))
    n <- length(par[[1]])

    Linvt <- list()
    Sig <- list()

    for (i in 1:n) {
      Linvt[[i]] <- matrix(0, nrow = k, ncol = k)
      for (j in 1:k) {
        Linvt[[i]][j, j] <- par[[paste0("innov", j)]][i]
        if (j < k) {
          for (l in (j+1):k) {
            Linvt[[i]][j, l] <- -par[[paste0("phi", j, l)]][i] *
                    par[[paste0("innov", l)]][i]          
	  }
        }
      }
    L <- solve(Linvt[[i]])
    Sig[[i]] <- t(L) %*% L
    }

    return(Sig)
  }

  rval$r <- function(par) {
    n <- length(par[[1]])
    k <- length(grep("innov", names(par)))
    Sigs <- rval$covariance(par)
    
    simdat <- matrix(ncol = k, nrow = n)

    for (i in 1:n) { 
      mu_vec <- vector(length = k)
      for (j in 1:k) {
        mu_vec[j] <- par[[paste0("mu", j)]][i]
      }
      simdat[i, ] <- mvtnorm::rmvnorm(1, mean = mu_vec, sigma = Sigs[[i]])
    }
    return(simdat)
  }
  

  # --- set class and return ---
  class(rval) <- "family.bamlss"
  rval
}


## LOG-LIKELIHOOD

# #' @param y amatrix n x k.
# #' @param par a list with k+k+k(k-1)/2 elements named like the parameters of the family,
# #'            each element is a numeric vector of length n.
log_dmvnmodchol_ref <- function(y, par) {
  n <- nrow(y) # number of observations
  k <- ncol(y) # dimension of gaussian distribution
  mu <- paste0("mu", seq_len(k))
  innov <- paste0("innov", seq_len(k))
  phi <- utils::combn(seq_len(k), 2, function(x) paste0("phi", x[1], x[2]))
  par_full <- do.call("cbind", par)
  term_1 <- -k / 2 * log(2 * pi)
  term_2 <- -1 / 2 * apply(log(subset(par_full, select = innov)), 1, sum)     
  term_3 <- vector(length = n) # initialise term_3 vector
  y_til <- as.matrix(y - subset(par_full, select = mu))
  for (ni in 1:n) {
    y_tild <- y_til[ni, ]
    L_inv <- matrix(0, nrow = k, ncol = k) # initialise L^-1 matrix
    for (l in 1:length(phi)) { # assign off diagonal values
      i <- utils::combn(k, 2)[1, l]
      j <- utils::combn(k, 2)[2, l]
      L_inv[j, i] <- -par[[paste0("phi", i, j)]][ni] /
	     sqrt(par[[paste0("innov", j)]][ni]) 
    }  
    for (i in 1:length(innov)) {
      L_inv[i, i] <- 1/sqrt(par[[paste0("innov", i)]][ni])
    }
    term_3[ni] <- -1 / 2 * norm(L_inv %*% y_tild, type = "2") ^ 2   
  }
  ll <- term_1 + term_2 + term_3
  return(ll)
}

log_dmvnmodchol_C <- function(y, par) {
  y <- as.matrix(y)
  storage.mode(y) <- "numeric"
  n <- nrow(y)
  k <- ncol(y)
  par <- do.call("cbind", par)
  .Call("log_dmvnmodcholC", y, par, n, k, PACKAGE = "bamlssMVN")
}

# choose `log_dmvnchol_ref` or `log_dmvnchol_C` for computing the log-density
log_dmvnmodchol <- log_dmvnmodchol_C


## BEGIN WITH SCORES

# #' @param j dimension of parameter
mu_score_mvnmodchol_ref <- function(y, par, j) {
  n <- nrow(y) # number of observations
  k <- ncol(y) # dimension of gaussian distribution
  mu <- paste0("mu", seq_len(k))
  innov <- paste0("innov", seq_len(k))
  phi <- utils::combn(seq_len(k), 2, function(x) paste0("phi", x[1], x[2]))
  par_full <- do.call("cbind", par)
  y_til <- as.matrix(y - subset(par_full, select = mu))
  dl_dmu <- vector(length = n)
  for (ni in 1:n) {
    y_tild <- y_til[ni, ]
    L_inv <- matrix(0, nrow = k, ncol = k) # initialise L^-1 matrix
    temp <- utils::combn(k, 2)
    for (l in 1:length(phi)) { # assign off diagonal values
      ii <- temp[1, l]
      jj <- temp[2, l]
      L_inv[jj, ii] <- -par[[paste0("phi", ii, jj)]][ni] /
	      sqrt(par[[paste0("innov", jj)]][ni])
    }
    for (ii in 1:length(innov)) {
      L_inv[ii, ii] <- 1 / sqrt(par[[paste0("innov", ii)]][ni])
    }
    Sigma_inv <- t(L_inv) %*% L_inv
    dl_dmu[ni] <- sum(Sigma_inv[j, ] * y_tild)
  }
  return(dl_dmu)
}      	

mu_score_mvnmodchol_C <- function(y, par, j) {
  y <- as.matrix(y)
  storage.mode(y) <- "numeric"
  n <- nrow(y)
  k <- ncol(y)
  par <- do.call("cbind", par)
  j <- as.integer(j)
  .Call("mu_score_mvnmodcholC", y, par, n, k, j, PACKAGE = "bamlssMVN")
}

# choose `mu_score_mvnchol_ref` or `mu_score_mvnchol_C` for computing mu-scores
mu_score_mvnmodchol <- mu_score_mvnmodchol_C


innov_score_mvnmodchol_ref <- function(y, par, j) {
  n <- nrow(y) # number of observations
  k <- ncol(y) # dimension of gaussian distribution
  mu <- paste0("mu", seq_len(k))
  innov <- paste0("innov", seq_len(k))
  phi <- utils::combn(seq_len(k), 2, function(x) paste0("phi", x[1], x[2]))
  par_full <- do.call("cbind", par)
  y_til <- as.matrix(y - subset(par_full, select = mu))
  dl_dinnov <- vector(length = n)
  for (ni in 1:n) {
    y_tild <- y_til[ni, ]
    Phimat <- matrix(0, nrow = k, ncol = k) # initialise -T matrix
    diag(Phimat) <- -1
    temp <- utils::combn(k, 2)
    for (l in 1:length(phi)) { # assign off diagonal values
      ii <- temp[1, l]
      jj <- temp[2, l]
      Phimat[jj, ii] <- par[[paste0("phi", ii, jj)]][ni]
    }

    dl_dinnov[ni] <- -1 / 2 + 1 / (2 * par[[paste0("innov", j)]][ni]) * 
	               sum(y_tild[1:j] * Phimat[j, 1:j])^2
  }
  return(dl_dinnov)
}


innov_score_mvnmodchol_C <- function(y, par, j) {
  y <- as.matrix(y)
  storage.mode(y) <- "numeric"
  n <- nrow(y)
  k <- ncol(y)
  par <- do.call("cbind", par)
  j <- as.integer(j)
  .Call("innov_score_mvnmodcholC", y, par, n, k, j, PACKAGE = "bamlssMVN")
}

# choose `mu_score_mvnchol_ref` or `mu_score_mvnchol_C` for computing mu-scores
innov_score_mvnmodchol <- innov_score_mvnmodchol_C



# #' @param i dimension of parameter
phi_score_mvnmodchol_ref <- function(y, par, i, j) {
  n <- nrow(y) # number of observations
  k <- ncol(y) # dimension of gaussian distribution
  mu <- paste0("mu", seq_len(k))
  innov <- paste0("innov", seq_len(k))
  phi <- utils::combn(seq_len(k), 2, function(x) paste0("phi", x[1], x[2]))
  par_full <- do.call("cbind", par)
  y_til <- as.matrix(y - subset(par_full, select = mu))
  dl_dphi <- vector(length = n)
  for (ni in 1:n) {
    y_tild <- y_til[ni, ]
    Tmat <- matrix(0, nrow = k, ncol = k) # initialise T matrix
    diag(Tmat) <- 1
    temp <- utils::combn(k, 2)
    for (l in 1:length(phi)) { # assign off diagonal values
      ii <- temp[1, l]
      jj <- temp[2, l]
      Tmat[jj, ii] <- -par[[paste0("phi", ii, jj)]][ni]
    }

    dl_dphi[ni] <- y_tild[i] / par[[paste0("innov", j)]][ni] *
                       sum(y_tild[1:j] * Tmat[j, 1:j])
  }
  return(dl_dphi) 

}

phi_score_mvnmodchol_C <- function(y, par, i, j) {
  y <- as.matrix(y)
  storage.mode(y) <- "numeric"
  n <- nrow(y)
  k <- ncol(y)
  par <- do.call("cbind", par)
  i <- as.integer(i)
  j <- as.integer(j)
  .Call("phi_score_mvnmodcholC", y, par, n, k, i, j, PACKAGE = "bamlssMVN")
}

# choose `mu_score_mvnchol_ref` or `mu_score_mvnchol_C` for computing mu-scores
phi_score_mvnmodchol <- phi_score_mvnmodchol_C


## BEGIN WITH HESSIAN

mu_hess_mvnmodchol_ref <- function(y, par, j) {
  n <- nrow(y) # number of observations
  k <- ncol(y) # dimension of gaussian distribution
  mu <- paste0("mu", seq_len(k))
  innov <- paste0("innov", seq_len(k))
  phi <- utils::combn(seq_len(k), 2, function(x) paste0("phi", x[1], x[2]))
  par_full <- do.call("cbind", par)
  dl_dmu <- vector(length = n)
  for (ni in 1:n) {
    dl_dmu[ni] <- 1 / par[[paste0("innov", j)]][ni]
    if (j < k) {
      for (jj in (j + 1):k) {
        dl_dmu[ni] <- dl_dmu[ni] + par[[paste0("phi", j, jj)]][ni] ^ 2 / 
		par[[paste0("innov", jj)]][ni]
      }
    } 
  }
  return(dl_dmu)
}      	

# choose `mu_score_mvnchol_ref` or `mu_score_mvnchol_C` for computing mu-scores
mu_hess_mvnmodchol <- mu_hess_mvnmodchol_ref



innov_hess_mvnmodchol_ref <- function(y, par, j) {
  n <- nrow(y) # number of observations
  k <- ncol(y) # dimension of gaussian distribution
  mu <- paste0("mu", seq_len(k))
  innov <- paste0("innov", seq_len(k))
  phi <- utils::combn(seq_len(k), 2, function(x) paste0("phi", x[1], x[2]))
  par_full <- do.call("cbind", par)
  y_til <- as.matrix(y - subset(par_full, select = mu))
  dl_dinnov <- vector(length = n)
  for (ni in 1:n) {
    y_tild <- y_til[ni, ]
    Phimat <- matrix(0, nrow = k, ncol = k) # initialise -T matrix
    diag(Phimat) <- -1
    temp <- utils::combn(k, 2)
    for (l in 1:length(phi)) { # assign off diagonal values
      ii <- temp[1, l]
      jj <- temp[2, l]
      Phimat[jj, ii] <- par[[paste0("phi", ii, jj)]][ni]
    }

    dl_dinnov[ni] <- 0.5 / par[[paste0("innov", j)]][ni] * 
	               sum(y_tild[1:j] * Phimat[j, 1:j]) ^ 2
  }
  return(dl_dinnov)
}

# choose `mu_score_mvnchol_ref` or `mu_score_mvnchol_C` for computing mu-scores
innov_hess_mvnmodchol <- innov_hess_mvnmodchol_ref





# #' @param i dimension of parameter
phi_hess_mvnmodchol_ref <- function(y, par, i, j) {
  n <- nrow(y) # number of observations
  k <- ncol(y) # dimension of gaussian distribution
  mu <- paste0("mu", seq_len(k))
  innov <- paste0("innov", seq_len(k))
  phi <- utils::combn(seq_len(k), 2, function(x) paste0("phi", x[1], x[2]))
  par_full <- do.call("cbind", par)
  y_til <- as.matrix(y - subset(par_full, select = mu))
  dl_dphi <- vector(length = n)
  for (ni in 1:n) {
    dl_dphi[ni] <- y_til[ni, i] ^ 2 / par[[paste0("innov", j)]][ni]
  }
  return(dl_dphi) 
}

# choose `mu_score_mvnchol_ref` or `mu_score_mvnchol_C` for computing mu-scores
phi_hess_mvnmodchol <- phi_hess_mvnmodchol_ref



