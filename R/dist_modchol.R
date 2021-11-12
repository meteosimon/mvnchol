## distree family for modified Choesky
## compare with dist_mvnorn l. 2208 in families.R of disttree pkg
#' @param k Integer, the dimension of the observation.
#' @param r Integer, the number of off-diagonals to model (AD-r covariance).
dist_mvn_modchol <- function(k, r = k - 1L, ...) {

    # --- set names of distributional parameters ---
    nms_mu <- paste0("mu_", seq_len(k))
    nms_innov <- paste0("innov_", seq_len(k))
    combns <- utils::combn(k, 2)
    combns_cut <- combns[, which((combns[2, ] - combns[1, ]) <= r)]
    nms_phi <- paste("phi", combns_cut[1, ], combns_cut[2, ], sep = "_")
    k_phi <- ncol(combns_cut) # number of phi parameters
    k_all <- k + k + k_phi    ## number of all parameters
    nms_par <- c(nms_mu, nms_innov, nms_phi)
    nms_eta <- c(nms_mu, paste0("log(", nms_innov, ")"), nms_phi)

    ## density
    ddist <-  function(y, eta, log = TRUE, weights = NULL, sum = FALSE) {
        n <- nrow(y) # number of observations
        term_1 <- -k / 2 * log(2 * pi)
        term_2 <- -1 / 2 * sum(eta[paste0("log(", nms_innov, ")")])
        term_3 <- numeric(n) # initialise term_3 vector
        y_til <- y - matrix(rep(eta[nms_mu], n), ncol = k, byrow = TRUE)
        for (ni in seq_len(n)) {
            y_tild <- y_til[ni, ]
            L_inv <- matrix(0, nrow = k, ncol = k) # initialise L^-1 matrix
            for (l in seq_len(k_phi)) { # assign off diagonal values    
                i <- combns_cut[1, l]
                j <- combns_cut[2, l]
                L_inv[j, i] <- -eta[[paste0("phi_", i, "_", j)]] /
                    sqrt(exp(eta[[paste0("log(innov_", j, ")")]]))
            }
            for (i in seq_len(k)) {
                L_inv[i, i] <- 1/sqrt(exp(eta[[paste0("log(innov_", i, ")")]]))
            }
            term_3[ni] <- -1 / 2 * norm(L_inv %*% y_tild, type = "2") ^ 2
        }
        ll <- term_1 + term_2 + term_3
        if (isTRUE(log)) {
            if (isFALSE(sum)) {
                return(ll)
            } else {
                return(sum(ll))
            }
        } else {
            if (isFALSE(sum)) {
                return(exp(ll))
            } else {
                message("You are summing the (non-log) likelihoods!!!")
                return(sum(exp(ll)))
            }
        }
    }

    ## scores
    sdist <- function(y, eta, weights = NULL, sum = FALSE) {
        n <- nrow(y)
        scores_mat <- matrix(ncol = k_all, nrow = n)  
        colnames(scores_mat) <- nms_eta
        y_til <- y - matrix(rep(eta[nms_mu], n), ncol = k, byrow = TRUE)

        for (ni in seq_len(n)) {               
            L_inv <- matrix(0, nrow = k, ncol = k)
            Phi_mat <- matrix(0, nrow = k, ncol = k)
            diag(Phi_mat) <- -1  
            for (l in seq_len(k_phi)) {    
                i <- combns_cut[1, l]
                j <- combns_cut[2, l]
                L_inv[j, i] <- -eta[[paste0("phi_", i, "_", j)]] /
                    sqrt(exp(eta[[paste0("log(innov_", j, ")")]]))
                Phi_mat[j, i] <- eta[[paste0("phi_", i, "_", j)]] 
            }
            for (i in seq_len(k)) { 
                L_inv[i, i] <- 1/sqrt(exp(eta[[paste0("log(innov_", i, ")")]]))
            }
            Sigma_inv <- t(L_inv) %*% L_inv

            # mu scores in first k columns
            for (j in seq_len(k)) {  
                scores_mat[ni, j] <- sum(Sigma_inv[j, ] * y_til[ni,])
            }

            # log(innov) scores in next k columns    
            for (j in seq_len(k)) {  
                scores_mat[ni, j + k] <-
                    -0.5 + 1 / (2 * exp(eta[[paste0("log(innov_", j, ")")]])) *
                    sum(y_til[ni, 1:j] * Phi_mat[j, 1:j])^2
            }

            # phi scores in columns 2k, ..., k_all
            for (l in seq_len(k_phi)) {
                i <- combns_cut[1, l]
                j <- combns_cut[2, l]
                scores_mat[ni, 2*k + l] <- - y_til[ni, i] / 
                    exp(eta[[paste0("log(innov_", j, ")")]]) *
                    sum(y_til[ni, 1:j] * Phi_mat[j, 1:j])
            }
        }
        if (isFALSE(sum)) {
            return(scores_mat)
        } else {
            return(colSums(scores_mat))
        }
    }

    ## links
    link <- c(
        rep("identity", k),
        rep("log", k),
        rep("identity", k_phi)
    )

    linkfun <- function(par) {
        eta <- par
        eta[seq_len(k) + k] <- log(par[seq_len(k) + k])
        names(eta) <- nms_eta
        return(eta)
    }

    linkinv <- function(eta) {
        par <- eta
        par[seq_len(k) + k] <- exp(eta[seq_len(k) + k])
        names(par) <- nms_par
        return(par)
    }
    
    linkinvdr <- function(eta) {
        dpardeta <- rep(1, k_all)
        dpardeta[seq_len(k) + k] <- exp(eta[seq_len(k) + k])
        names(dpardeta) <- nms_par
        return(dpardeta)
    }

    ## MLE
    startfun <- function(y, weights = NULL) {
        n <- nrow(y)
        if(is.null(weights) || (length(weights) == 0L)) {
            starteta <- rep(-999, k_all)
            names(starteta) <- nms_eta 
            starteta[nms_mu] <- colMeans(y) # Estimates for mus
          
            y_til <- y - matrix(rep(starteta[nms_mu], n), ncol = k, byrow = TRUE) 
            starteta[["log(innov_1)"]] <- var(y_til[, 1])
            for (i in 2:k) {
                X <- y_til[, max(1, r-k):(i-1)] 
                beta_vec <- solve(t(X) %*% X) %*% t(X) %*% y_til[, i]            
                starteta[paste0("phi_", max(1, r-k):(i-1) ,"_", i)] <- beta_vec
    
                starteta[paste0("log(innov_", i, ")")] <-
                    log(var(y_til[, i] - as.vector(X %*% beta_vec)))
            } 
        } else {
            stop("Weights not implemented for startfun")
        }
        return(starteta)
    }
  
    mle <- TRUE

    dist_list <- list(family.name = "MVN Cholesky",
        ddist = ddist,
        sdist = sdist,
        link = link,
        linkfun = linkfun,
        linkinv = linkinv,
        linkinvdr = linkinvdr,
        startfun = startfun,
        mle = mle,
        gamlssobj = FALSE,
        censored = FALSE
    )
  
    # Return family object
    class(dist_list) <- "disttree.family"
    return(dist_list)
}




