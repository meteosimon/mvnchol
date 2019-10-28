########################################################################
##
## Load meteo data --- IBK obs and ECMWF-EPS --- and builds a
##  data.frame with ensemble mean, sd, and correlations of the
##  ensemble between consecutive lead times.
##  Here, only predictions for 12~UTC are considered.
##
##  T. Simon, JAN 2017
##
########################################################################

########################################################################
## Libraries
library("STAGE")
library("zoo")
rhogit <- function(r) { r / sqrt(1 - r^2) }
########################################################################

########################################################################
## Parameters
Sys.setenv("TZ"="UTC")
# id <- 11120 ## IBK
id <- 10147 ## Hamburg
n.pred <- 10     ## Duration of prediction in days
steps <- seq(12, by=24, length=n.pred)
begin <- as.POSIXct("2010-01-01")
end <- as.POSIXct("2017-01-01") + n.pred*24*3600
########################################################################

########################################################################
## Load obs
obs <- sybase(type = "synop", station = id, parameter = c("T"),
              begin = begin, end = end)
## pick values at high noon
obs <- subset(obs, subset= format(index(obs), "%H")=="12")
########################################################################

########################################################################
## Load dmo
d <- list()
dnames <- sprintf("step.%03d", steps)
for( j in seq_along(steps) ) {
  ## load data
  d[[dnames[j]]] <- dmo(model = "ECEPS", runhour = 0, step = steps[j],
           station = id, parameter = c("t2m"),
           begin = begin, end = end)
  d[[dnames[j]]]$t2m <- d[[dnames[j]]]$t2m - 273.15
}
########################################################################

########################################################################
## build data.frame

## obs
df <- obs
for ( j in seq(n.pred-1) ) {
  df <- cbind(df, lag(obs, k=j))
}
df <- cbind(df, lag(obs, k=-1))
cnames <- c(sprintf("T.%03d", steps), "T.lag12") 
colnames(df) <- cnames

## add doy, year
df$doy <- as.POSIXlt(index(obs))$yday + 1
df$year <- as.POSIXlt(index(obs))$year + 1900

index(df) <- index(df) - 12*3600

## means
mnames <- sprintf("m.%03d", steps) 
for ( j in seq(n.pred) ) {
  x <- tapply(d[[dnames[j]]]$t2m, d[[dnames[j]]]$rundate, mean)
  x <- zoo(x, order.by=as.POSIXct(names(x), "%Y%m%d", tz="UTC"))
  df <- merge(df, x)
  names(df)[ncol(df)] <- mnames[j]
}

## sds
snames <- sprintf("s.%03d", steps) 
for ( j in seq(n.pred) ) {
  x <- tapply(d[[dnames[j]]]$t2m, d[[dnames[j]]]$rundate, sd)
  x <- zoo(x, order.by=as.POSIXct(names(x), "%Y%m%d", tz="UTC"))
  df <- merge(df, x)
  names(df)[ncol(df)] <- snames[j]
}

## sds
logsnames <- sprintf("logs.%03d", steps) 
for ( j in seq(n.pred) ) {
  x <- log( df[,snames[j]] )
  df <- merge(df, x)
  names(df)[ncol(df)] <- logsnames[j]
}

## rhos
norder <- paste0("t2m.",0:50)
ltime <- sprintf(".%03d", steps)
for ( j in seq(n.pred) ) {
 for ( i in seq(n.pred) ) {
  if ( i > j ) {
   ## 
   x <- reshape(subset(d[[dnames[j]]],
                       select=c("rundate", "member", "t2m")),
                direction="wide", timevar="member", idvar="rundate")
   ind <- grep("t2m", names(x))
   x$rundate <- as.character(x$rundate)
   x <- zoo(x[,ind], order.by=as.POSIXct(x$rundate, "%Y%m%d", tz="UTC"))
   ind <- match(norder, colnames(x))
   x <- x[,ind]
   ##
   y <- reshape(subset(d[[dnames[i]]],
                       select=c("rundate", "member", "t2m")),
                direction="wide", timevar="member", idvar="rundate")
   ind <- grep("t2m", names(y))
   y$rundate <- as.character(y$rundate)
   y <- zoo(y[,ind], order.by=as.POSIXct(y$rundate, "%Y%m%d", tz="UTC"))
   ind <- match(norder, colnames(y))
   y <- y[,ind]
   ##
   xcor <- apply(cbind(x,y), 1, function(x) cor(x[1:51],x[52:102]))
   xcor <- zoo(xcor, order.by=as.POSIXct(names(xcor), tz="UTC"))
   df <- merge(df, xcor)
   names(df)[ncol(df)] <- paste0("r",ltime[j],ltime[i])
  }
 }
}

## rhogit transform
colid <- grep("r.", names(df))
rdata <- df[,colid]
names(rdata) <- paste0("r", names(rdata))
rdata <- rhogit(rdata)
df <- cbind(df, rdata)
rm(colid, rdata)
########################################################################

########################################################################
## save data
data <- na.omit(df)
save(data, file=paste0("../data/ID",id,".TEMP.EPS.lead", n.pred, ".rda"))
########################################################################

