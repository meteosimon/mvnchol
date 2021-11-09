## distree family for modified Choesky
dist_mvn_modchol <- function(k, ...) {
    stop("distree family for modified Choesky is not implemented yet!")

    ## TODO: density
    ddist <-  function(y, eta, log = TRUE, weights = NULL, sum = FALSE) {

    }

    ## TODO: scores
    sdist <- function(y, eta, weights = NULL, sum = FALSE) {

    }

    ## links TODO: must be conditioned on k
    link <- c("identity", "log")

    linkfun <- function(par) {
        eta <- c(par[1], log(par[2]))
        names(eta) <- c("mu", "log(sigma)")
        return(eta)
    }

    linkinv <- function(eta) {
        par <- c(eta[1], exp(eta[2]))
        names(par) <- c("mu", "sigma")
        return(par)
    }
    
    linkinvdr <- function(eta) {
        dpardeta <- c(1, exp(eta[2]))
        names(dpardeta) <- c("mu", "sigma")
        return(dpardeta)
    }

    ## TODO: startfun should give MLE for all distributional parameters
    ##       (code comes from univariate gaussian dittree family)
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




