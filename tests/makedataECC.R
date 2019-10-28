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
########################################################################

########################################################################
## Parameters
Sys.setenv("TZ"="UTC")
n.pred <- 6     ## Duration of prediction in days
steps <- seq(12, by=24, length=n.pred)
begin <- as.POSIXct("2010-01-01")
end <- as.POSIXct("2017-01-01") + n.pred*24*3600
########################################################################

########################################################################
## Load dmo
d <- NULL
df <- NULL
norder <- paste0("t2m.",0:50)
dnames <- sprintf("step.%03d", steps)
for( j in seq_along(steps) ) {
  ## load data
  d <- dmo(model="ECEPS", runhour=0, step=steps[j],
           station=11120, parameter=c("t2m"),
           begin=begin, end=end)
  d$t2m <- d$t2m - 273.15

  x <- reshape(subset(d, select=c("rundate", "member", "t2m")),
                direction="wide", timevar="member", idvar="rundate")
  ind <- grep("t2m", names(x))
  x$rundate <- as.character(x$rundate)
  x <- zoo(x[,ind], order.by=as.POSIXct(x$rundate, "%Y%m%d", tz="UTC"))
  ind <- match(norder, colnames(x))
  x <- x[,ind]
  names(x) <- sprintf("%03d.%d", steps[j], 0:50)

  if ( j==1 ) {
   df <- x 
  } else {
   df <- cbind(df, x)
  }
}
########################################################################

########################################################################
## compute order of data
########################################################################

########################################################################
## save data
data <- na.omit(df)
save(data, file=paste0("../data/IBK.TEMP.EPS.lead", n.pred, ".ORDER.rda"))
########################################################################

