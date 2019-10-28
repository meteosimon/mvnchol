## rhogit linkfun
linkinv.rhogit <- function(eta) { eta / sqrt(1 + eta^2) }
linkfun.rhogit <- function(mu) { mu / sqrt(1 - mu^2) }

## PITs
hist.pit <- function(p, ...) {
  h <- hist(p, freq=FALSE, col="lightgray", las=1, ...)
  abline(h=1, col=2, lty=2, lwd=2)

  invisible(h)
}

pit.bamlss.mvnorm <- function(object, ltime="012", newdata=NULL, 
                              mstop=NULL, ...) {
  response <- paste0("T.",ltime)
  ind1 <- grep(response, names(object$model.frame))
  mu <- paste0("mu",ind1)
  sig <- paste0("sigma",ind1)

  if ( is.null(mstop) )
    mstop <- object$model.stats$optimizer$boost.summary$mstop

  fit <- predict(object, type="parameter", newdata=newdata,
                 mstop=mstop)

  if ( is.null(newdata) ) newdata <- object$model.frame
  ind2 <- which(names(newdata)==response)
  y <- newdata[,ind2]

  p <- pnorm(y, mean=fit[[mu]], sd=fit[[sig]])
  h <- hist.pit(p, main=paste0("MVN, lead time ",ltime,"h"), ...)

  invisible(h)
}

pit.bamlss.gaussian <- function(object, newdata=NULL, ...) {
  response <- attr(formula(object), "response.name")
  ltime <- substr(response, 3, 5)
  fit <- predict(object, type="parameter", newdata=newdata)

  if ( is.null(newdata) ) newdata <- object$model.frame
  ind2 <- which(names(newdata)==response)
  y <- newdata[,ind2]
  p <- pnorm(y, mean=fit$mu, sd=fit$sigma)
  h <- hist.pit(p, main=paste0("NGR, lead time ",ltime,"h"), ...)

  invisible(h)
}

pit.crch <- function(object, ...) {
  ## TODO: add option newdata for out-of-sample verification
  response <- all.vars(formula(object))[1]
  ltime <- substr(response, 3, 5)
  m <- predict(object, type="location")
  s <- predict(object, type="scale")

  p <- pnorm(object$model[[response]], mean=m, sd=s)
  h <- hist.pit(p, main=paste0("MVN, lead time ",ltime), ...)

  invisible(h)
}

## Sharpness Diagram
SharpnessDiagram <- function(bamlssobj, modellist, coverage=.9, ...) {

  k <- length(modellist)
  boxes <- NULL
  upper <- ( 1 + coverage ) / 2  
  lower <- ( 1 - coverage ) / 2  
  steps <- substr(names(modellist), 5, 7)

  fit <- predict(bamlssobj, type="parameter")
  snames <- paste0("sigma",seq(k))
  for ( j in seq(k) ) {
    sig <- fit[[snames[j]]]
    int <- qnorm(upper, mean=0, sd=sig) - 
           qnorm(lower, mean=0, sd=sig)
    boxes <- cbind(boxes, int)
  }

  for ( j in seq(k) ) {
    if ( class(modellist[[1]])=="crch" )
      sig <- predict(modellist[[j]], type="scale")

    if ( class(modellist[[1]])=="bamlss" )
      sig <- predict(modellist[[j]], type="parameter")$sigma

    int <- qnorm(upper, mean=0, sd=sig) - 
           qnorm(lower, mean=0, sd=sig)
    boxes <- cbind(boxes, int)
  }
  colnames(boxes) <- c(paste0("bamlss.",steps), 
                       paste0("crch.",steps))

  ## re-order matrix
    boxes <- boxes[,rep(1:k, each=2) + rep(c(0,k), k)]

  ## plot
    names <- rep(c("MVN", "NGR"), k)
    boxplot(boxes, names=names, col="lightgray", ...)
    abline(v=seq(2.5, by=2, length=k-1))
    axis(side = 3, at = seq(1.5, by=2, length=k),
         labels = steps, tick=FALSE, padj=1)

  invisible(boxes)
}

## Scoring Rules
## CRPS
crps.norm <- function(y, mu, sig) {
  z <- (y - mu)/sig
  rval <- z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 1/sqrt(pi)
  rval <- sig * rval
  return(rval)
}

crps.bamlss <- function(object, ltime="012", newdata=NULL, ...) {
  response <- paste0("T.",ltime)
  ind1 <- grep(response, names(object$model.frame))
  mu <- paste0("mu",ind1)
  sig <- paste0("sigma",ind1)

  fit <- predict(object, type="parameter", newdata=newdata, ...)

  if ( is.null(newdata) ) newdata <- object$model.frame
  ind2 <- which(names(newdata)==response)
  y <- newdata[,ind2]

  rval <- crps.norm(y, mu=fit[[mu]], sig=fit[[sig]])
  return(rval)
}

crps.bamlss.gaussian <- function(object, newdata=NULL) {
  response <- attr(formula(object), "response.name")
  ltime <- substr(response, 3, 5)
  fit <- predict(object, type="parameter", newdata=newdata)

  if ( is.null(newdata) ) newdata <- object$model.frame
  ind2 <- which(names(newdata)==response)
  y <- newdata[,ind2]

  rval <- crps.norm(y, mu=fit$mu, sig=fit$sigma)
  return(rval)
}

crps.crch <- function(object) {
  response <- all.vars(formula(object))[1]
  ltime <- substr(response, 3, 5)
  m <- predict(object, type="location")
  s <- predict(object, type="scale")

  rval <- crps.norm(object$model[[response]], mu=m, sig=s)
  return(rval)
}

## Energy Score
es.sample <- function(y, X, X_prime=NULL) {
  n <- nrow(y)
  if ( is.null(X_prime) )
    X_prime <- X[,sample.int(n)]

  diff <- X - y
  term1 <- mean(sqrt(rowSums(diff^2)))

  diff <- X - X_prime
  term2 <- .5*mean(sqrt(rowSums(diff^2)))

  rval <- term1 - term2
  return(rval)
}

es.mvnorm <- function(y, mu, Sigma, n=10000) {
  require("mvtnorm")
  X <- rmvnorm(n, mean=mu, sigma=Sigma)
  X_prime <- rmvnorm(n, mean=mu, sigma=Sigma)

  diff <- X - y
  term1 <- mean(sqrt(rowSums(diff^2)))

  diff <- X - X_prime
  term2 <- .5*mean(sqrt(rowSums(diff^2)))

  rval <- term1 - term2
  return(rval)
}

es.bamlss <- function(object, newdata=NULL, nmcmc=10000, ...) {
  respnames <- names(object$y)
  if ( is.null(newdata) ) {
    y <- as.matrix(object$y)
  } else {
    y <- as.data.frame(subset(newdata, select=respnames))
  }
  n <- nrow(y)
  k <- ncol(y)

  fit <- predict(object, type="parameter", newdata=newdata, ...)
  par <- do.call("cbind", fit)
  cn <- colnames(par)
  sj <- grep("sigma", cn)
  mj <- grep("mu", cn)
  rj <- grep("rho", cn)

  rval <- NULL
  if ( length(rj)==1 & k>2 ) {
    for ( j in seq(n) ) {
      mu <- par[j,mj]
      D <- diag(par[j,sj])
      R <- diag(1,k)
      for ( ii in 1:(k-1) ) {
        for ( jj in 2:k ) {
          if ( jj > ii ) {
            R[jj,ii] <- par[j,rj]^(jj-ii)
            R[ii,jj] <- R[jj,ii]
          }
        }
      }
      Sigma <- D %*% R %*% D
      rval[j] <- es.mvnorm(y[j,], mu, Sigma, n=nmcmc)      
    }
  } else {
    for ( j in seq(n) ) {
      mu <- par[j,mj]
      D <- diag(par[j,sj])
      R <- diag(1,k)
      R[lower.tri(R)] <- par[j,rj]
      R <- t(R)
      R[lower.tri(R)] <- par[j,rj]
      Sigma <- D %*% R %*% D
 
      rval[j] <- es.mvnorm(y[j,], mu, Sigma, n=nmcmc)
    }
  }
  return(rval)
}

es.bAR <- function(object, newdata=NULL, nmcmc=10000, ...) {
  respnames <- names(object$y)
  if ( is.null(newdata) ) {
    y <- as.matrix(object$y)
  } else {
    y <- as.data.frame(subset(newdata, select=respnames))
  }
  n <- nrow(y)
  k <- ncol(y)

  fit <- predict(object, type="parameter", newdata=newdata, ...)
  par <- do.call("cbind", fit)
  cn <- colnames(par)
  sj <- grep("sigma", cn)
  mj <- grep("mu", cn)
  rj <- grep("rho", cn)

  rval <- NULL
  if ( length(rj)==1 & k>2 ) {
    for ( j in seq(n) ) {
      mu <- par[j,mj]
      D <- diag(par[j,sj])
      R <- diag(1,k)
      for ( ii in 1:(k-1) ) {
        for ( jj in 2:k ) {
          if ( jj > ii ) {
            R[jj,ii] <- par[j,rj]^(jj-ii)
            R[ii,jj] <- R[jj,ii]
          }
        }
      }
      Sigma <- D %*% R %*% D
      rval[j] <- es.mvnorm(y[j,], mu, Sigma, n=nmcmc)      
    }
  } else {
    for ( j in seq(n) ) {
      mu <- par[j,mj]
      D <- diag(par[j,sj])
      R <- diag(1,k)
      R[lower.tri(R)] <- par[j,rj]
      R <- t(R)
      R[lower.tri(R)] <- par[j,rj]
      Sigma <- D %*% R %*% D
 
      rval[j] <- es.mvnorm(y[j,], mu, Sigma, n=nmcmc)
    }
  }
  return(rval)
}


es.bamlsslist <- function(modellist, newdata=NULL, ...) {
  k <- length(modellist)
  #n <- nrow(modellist[[1]]$y)
  y <- mu <- sig <- NULL
  for ( j in seq(k) ) {
    response <- attr(formula(modellist[[j]]), "response.name")
    if ( is.null(newdata) ) {
      y <- cbind(y, modellist[[j]]$y[[response]])
    } else {
      y <- cbind(y, as.numeric(newdata[,response]))
    }
    fit <- predict(modellist[[j]], type="parameter",
                   newdata=newdata, ...)
    mu <- cbind(mu, fit$mu)
    sig <- cbind(sig, fit$sigma) 
  }

  n <- nrow(y)
  rval <- NULL
  for ( j in seq(n) ) {
    rval[j] <- es.mvnorm(y[j,], mu[j,], Sigma=diag(sig[j,]^2))
  }
  return(rval)
}

es.ECC <- function(modellist, ECCdata, newdata=NULL, 
                   nens=51, ...) {
  k <- length(modellist)
  #n <- nrow(modellist[[1]]$y)
  y <- mu <- sig <- NULL
  for ( j in seq(k) ) {
    response <- attr(formula(modellist[[j]]), "response.name")
    if ( is.null(newdata) ) {
      y <- cbind(y, modellist[[j]]$y[[response]])
    } else {
      y <- cbind(y, as.numeric(newdata[,response]))
    }
    fit <- predict(modellist[[j]], type="parameter",
                   newdata=newdata, ...)
    mu <- cbind(mu, fit$mu)
    sig <- cbind(sig, fit$sigma)
  }

  if ( !is.null(newdata) )
    ECCdata <- subset(ECCdata, 
       subset= index(ECCdata) %in% index(newdata))

  n <- nrow(y)
  probs <- (1:nens) / (nens+1)
  nsample <- 1000
  rval <- NULL
  for ( j in seq(n) ) {
    ## margins
    E <- matrix(qnorm(rep(probs,k), mean=rep(mu[j,], each=nens),
                sd=rep(sig[j,], each=nens)), ncol=k, nrow=nens)
    
    ## re-order
    ed <- matrix(ECCdata[j,], ncol=k, nrow=nens)
    ed <- apply(ed, 2, order)
    for ( i in seq(k) ) E[,i] <- E[ed[,i],i]

    ## sample
    X <- E[sample.int(nens, nsample, replace=TRUE),]
    for ( i in seq(k) )
      X[,i] <- X[,i] + rnorm(nsample, sd=(sig[j,i]^2)/(nens-1))
    
    X_prime <- E[sample.int(nens, nsample, replace=TRUE),]
    for ( i in seq(k) )
      X_prime[,i] <- X_prime[,i] + rnorm(nsample,
                                sd=(sig[j,i]^2)/(nens-1))

    ## evaluate   
    rval[j] <- es.sample(y[j,], X=X, X_prime=X_prime)
  }
  return(rval)
}

es.crch <- function(crchlist) {
  k <- length(crchlist)
  n <- crchlist[[1]]$n
  y <- mu <- sig <- NULL
  for ( j in seq(k) ) {
    response <- all.vars(formula(crchlist[[j]]))[1]
    y <- cbind(y, crchlist[[j]]$model[[response]])
    mu <- cbind(mu, predict(crchlist[[j]], type="location"))
    sig <- cbind(sig, predict(crchlist[[j]], type="scale")) 
  }

  rval <- NULL
  for ( j in seq(n) ) {
    rval[j] <- es.mvnorm(y[j,], mu[j,], Sigma=diag(sig[j,]^2))
  }
  return(rval)
}

## Dawid and Sebastiani or LogLik
llscore.bamlss <- function(object, newdata=NULL, ...) {
  require("mvtnorm")
  respnames <- names(object$y)
  y <- subset(newdata, select=respnames)
  if ( is.null(newdata) ) y <- as.matrix(object$y)
  n <- nrow(y)
  k <- ncol(y)

  fit <- predict(object, type="parameter", newdata=newdata, ...)
  par <- do.call("cbind", fit)
  cn <- colnames(par)
  sj <- grep("sigma", cn)
  mj <- grep("mu", cn)
  rj <- grep("rho", cn)

  rval <- NULL
  if ( length(rj)==1 & k>2 ) {
    for ( j in seq(n) ) {
      mu <- par[j,mj]
      D <- diag(par[j,sj])
      R <- diag(1,k)
      for ( ii in 1:(k-1) ) {
        for ( jj in 2:k ) {
          if ( jj > ii ) {
            R[jj,ii] <- par[j,rj]^(jj-ii)
            R[ii,jj] <- R[jj,ii]
          }
        }
      }
      Sigma <- D %*% R %*% D
      rval[j] <- dmvnorm(y[j,], mu, Sigma, log=TRUE)      
    }
  } else {
    for ( j in seq(n) ) {
      mu <- par[j,mj]
      D <- diag(par[j,sj])
      R <- diag(1,k)
      R[lower.tri(R)] <- par[j,rj]
      R <- t(R)
      R[lower.tri(R)] <- par[j,rj]
      Sigma <- D %*% R %*% D
 
      rval[j] <- dmvnorm(y[j,], mu, Sigma, log=TRUE)
    }
  }
  return(rval)
}

llscore.bamlsslist <- function(modellist, newdata=NULL, ...) {
  require("mvtnorm")
  k <- length(modellist)
  #n <- nrow(modellist[[1]]$y)
  y <- mu <- sig <- NULL
  for ( j in seq(k) ) {
    response <- attr(formula(modellist[[j]]), "response.name")
    if ( is.null(newdata) ) {
      y <- cbind(y, modellist[[j]]$y[[response]])
    } else {
      y <- cbind(y, as.numeric(newdata[,response]))
    }
    fit <- predict(modellist[[j]], type="parameter",
                   newdata=newdata, ...)
    mu <- cbind(mu, fit$mu)
    sig <- cbind(sig, fit$sigma) 
  }

  n <- nrow(y)
  rval <- NULL
  for ( j in seq(n) ) {
    rval[j] <- dmvnorm(y[j,], mu[j,], sigma=diag(sig[j,]^2), log=TRUE)
  }
  return(rval)
}

llscore.crch <- function(crchlist) {
  require("mvtnorm")
  k <- length(crchlist)
  n <- crchlist[[1]]$n
  y <- mu <- sig <- NULL
  for ( j in seq(k) ) {
    response <- all.vars(formula(crchlist[[j]]))[1]
    y <- cbind(y, crchlist[[j]]$model[[response]])
    mu <- cbind(mu, predict(crchlist[[j]], type="location"))
    sig <- cbind(sig, predict(crchlist[[j]], type="scale")) 
  }

  rval <- NULL
  for ( j in seq(n) ) {
    rval[j] <- dmvnorm(y[j,], mu[j,], sigma=diag(sig[j,]^2), log=TRUE)
  }
  return(rval)
}

## Plot Cov or Cor Matrix
table.covmat <- function(object, inv=FALSE, prob=.5, ...) {
  k <- ncol(as.matrix(object$y))
  
  fit <- predict(object, type="parameter")
  par <- do.call("cbind", fit)
  cn <- colnames(par)
  rj <- grep("rho", cn)
   
  q <- apply(par[,rj], 2, quantile, probs=prob)
  mat <- diag(1,k)
  mat[lower.tri(mat)] <- q
  mat <- t(mat)
  mat[lower.tri(mat)] <- q

  if (inv) {
    mat <- chol2inv(chol(mat))
  }

  colnames(mat) <- names(object$y)
  rownames(mat) <- names(object$y)

  rval <- xtable(round(mat,2))
  return(rval)
}

## Case study

