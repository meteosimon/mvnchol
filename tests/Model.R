########################################################################
##
## Fit a MVNorm via Boosting with predictors for all  mu, sig
##  and rho.
##
##  T. Simon, JAN 2017
##
########################################################################

########################################################################
## Libraries
library("zoo")
library("bamlss")
########################################################################

########################################################################
## Funcitons
#source("functions.R")
########################################################################

########################################################################
## Parameters
cat(" *  Set Parameter...\n")
args <- commandArgs(TRUE)
Sys.setenv("TZ"="UTC")
k <- 6     ## Duration of prediction in days
steps <- seq(12, by=24, length=k)
loyear <- 2016 # as.numeric(args[1])
maxit <- 50000
nu <- 0.05
########################################################################

########################################################################
## Load Data
cat(" *  Load Data...\n")
id <- 11120
load(paste0("../data/ID",id,".TEMP.EPS.lead10.rda"))
data$doy <- as.numeric(format(index(data), "%j"))
data$year <- as.numeric(format(index(data), "%Y"))
data <- subset(data, subset= year!=2017)
#data <- subset(data, subset= year!=2016)
#testdata <- NULL
testdata <- subset(data, subset= year==loyear)
data <- subset(data, subset= year!=loyear)
########################################################################

########################################################################
## Formula for MVN
cat(" *  Define Formula...\n")
ltime <- sprintf(".%03d", steps)

mu <- NULL
sig <- NULL
for ( j in seq(k) ) {
 mu[j] <- paste0("T",ltime[j],
                 " ~ s(doy, bs='cc') + s(doy, bs='cc', by=m",ltime[j], ")")
 sig[j] <- paste0("sigma",j,
                 " ~ s(doy, bs='cc') + s(doy, bs='cc', by=logs",ltime[j], ")")
# mu[j] <- paste0("T",ltime[j],
#                 " ~ s(doy, bs='cc') + s(doy, bs='cc', by=m",ltime[j], ")")
# sig[j] <- paste0("sigma",j,
#                 " ~ s(doy, bs='cc') + s(doy, bs='cc', by=logs",ltime[j], ")")
}
mu[1] <- paste0(mu[1], " + T.lag12")
m <- paste(mu, collapse=", ")
s <- paste(sig, collapse=", ")

rho <- NULL
l <- 0
for ( j in seq(k) ) {
 for ( i in seq(k) ) {
  if ( i > j ) {   
   l <- l + 1
   rho[l] <- paste0("rho",j,i," ~ rr",ltime[j],ltime[i])
#   rho[l] <- paste0("rho",j,i," ~ s(r",ltime[j],ltime[i],")")
  }
 }
}
r <- paste(rho, collapse=", ")
formtext <- paste0("list(", m, ", ", s, ", ", r, ")")

f <- eval(parse(text=formtext))
fAR <- f[1:(2*k)]
fAR[[(2*k+1)]] <- rho ~ s(doy, bs="cc")
########################################################################

########################################################################
## Models
## With seasonal cycles for the mean parameters
cat(" *  Estimate MVN...\n")
b <- bamlss(f, family = gF("mvnorm", k = k), data = data,
             optimizer = boost, sampler = FALSE,
             maxit = maxit, nu = nu, scale.d = TRUE)

b6000 <- bamlss(f, family = gF("mvnorm", k = k), data = data,
             optimizer = boost, sampler = FALSE,
             maxit = 6000, nu = .05, scale.d = TRUE)

## AR
cat(" *  Estimate AR1...\n")
bAR <- bamlss(fAR, family = gF("mvnormAR1", k = k), data = data,
             optimizer = boost, sampler = FALSE,
             maxit = 25000, nu = 0.1, scale.d = TRUE)

## Reference NGRs
cat(" *  Estimate NGR...\n")
bnames <- sprintf("step%03d", steps)
b.ngr <- list()
for ( j in seq(k) ) {
#  fotext <- paste0("list(T",ltime[j]," ~ m",ltime[j],
#                   " + s(doy, bs='cc') + s(doy, bs='cc', by=m",ltime[j], ")")
#  if ( j==1 ) fotext <- paste0(fotext, " + T.lag12")
#  fotext2 <- paste0(", sigma ~ logs",ltime[j],
#                   " + s(doy, bs='cc') + s(doy, bs='cc', by=logs",ltime[j], "))")
#  fo <- eval(parse(text=paste0(fotext,fotext2)))

  fo <- list(f[[j]], update(f[[k+j]], sigma ~ .))

  b.ngr[[bnames[j]]] <- bamlss(fo, data=data,        
                         optimizer = boost, sampler = FALSE,
                         maxit = 5000, nu = 0.1, scale.d = TRUE)
}
########################################################################

########################################################################
## Save Models
cat(" *  Save Models...\n")
save(b, b6000, bAR, b.ngr, testdata,
     file=paste0("../output/Models_ID",id,"_k6_loy_",loyear,".Rda"))
cat("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
########################################################################


