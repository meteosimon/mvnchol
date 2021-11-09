## distree family for modified Choesky
## compare with dist_mvnorn l. 2208 in families.R of disttree pkg
dist_mvn_modchol <- function(k, ...) {
    stop("distree family for modified Choesky is not implemented yet!")

    ## klt <- k * (k - 1L) / 2L  ## number of elements in the lower triangle of Sigma (and similar matrices)
    ## kn <- k + k + klt ## number of all distributional parameters

    # --- set names of distributional parameters ---
    mu <- paste0("mu_", seq_len(k))
    innov <- paste0("innov_", seq_len(k))
    phi <- utils::combn(seq_len(k), 2, function(x) paste("phi", x[1], x[2], sep = "_"))
    k_phi <- k * (k-1) / 2    ## number of phi parameters
    k_all <- k + k + k_phi    ## number of all parameters
    nms_par <- c(mu, innov, phi)
    nms_eta <- c(mu, paste0("log(", innov, ")"), phi)

    ## TODO: density
    ## Input:
    ## y ...... observed samples (matrix)
    ## eta .... named numeric vector (one entry for each distributional parameter)
    ## Output:
    ## numeric vector with contrbutions to density 
    ddist <-  function(y, eta, log = TRUE, weights = NULL, sum = FALSE) {

    }

    ## TODO: scores
    ## Output:
    ## matrix (n x k) with scores contrbutions
    sdist <- function(y, eta, weights = NULL, sum = FALSE) {

    }

    ## links TODO: must be conditioned on k
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

    ## TODO: startfun should give MLE for all distributional parameters
    ##       (code comes from univariate gaussian dittree family)
    ## Output:
    ## vector with ML estimates of distributional parameters
    startfun <- function(y, weights = NULL) {
        if(is.null(weights) || (length(weights) == 0L)) {
            mu <- mean(y)
            sigma <- sqrt(1/length(y) * sum((y - mu)^2))
        } else {
            mu <- weighted.mean(y, weights)
            sigma <- sqrt(1/sum(weights) * sum(weights * (y - mu)^2))
        }
        starteta <- c(mu, log(sigma))
        names(starteta) <- c("mu", "log(sigma)")
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




