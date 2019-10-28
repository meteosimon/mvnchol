## Example how to plot rescaled varying coefficient
library("bamlss")
#load("../output/Models_ID10147_k6_loy_2016.Rda")
load("../output/Models_ID11120_k6_loy_2016.Rda")

################################################################################
## Function to plot seasonal effect of varying coefficient
plotVarCoef <- function(object, xvar = "doy", yvar = "m.012", add = FALSE, ...) {
  sc <- attr(object$model.frame, "scale")
  rng <- range(object$model.frame[xvar])
  rng <- sc$center[xvar] + rng*sc$scale[xvar]

  x <- object$results$mu$s.effects[[paste0("s(",xvar,",by=",yvar,")")]]
  if ( add ) {
    lines(sc$center[xvar] + sort(x[[xvar]])*sc$scale[xvar], 
      x$Mean[order(x[[xvar]])] / sc$scale[yvar], ...)   
  } else {
    plot(sc$center[xvar] + sort(x[[xvar]])*sc$scale[xvar], 
      x$Mean[order(x[[xvar]])] / sc$scale[yvar], xlim=rng, type = "l", ...)
  }
}

## Function to plot seasonal effect of 'intercept'
plotSeason <- function(object, xvar = "doy", add = FALSE, ...) {
  sc <- attr(object$model.frame, "scale")
  rng <- range(object$model.frame[xvar])
  rng <- sc$center[xvar] + rng*sc$scale[xvar]

  x <- object$results$mu$s.effects[[paste0("s(",xvar,")")]]
  if ( add ) {
    lines(sc$center[xvar] + sort(x[[xvar]])*sc$scale[xvar], 
      x$Mean, ...)   
  } else {
    plot(sc$center[xvar] + sort(x[[xvar]])*sc$scale[xvar], 
      x$Mean, xlim=rng, type = "l", ...)
  }
}

################################################################################

##
cols <- rainbow_hcl(6)
yvars <- paste0("m.", sprintf("%03d", seq(12, 132, by = 24)))
plotVarCoef(b.ngr[[1]], yvar = "m.012", ylim = c(.5, 1.1),
            col = cols[1], lwd = 2, xaxs = "i", axes = FALSE,
            xlab = "Time of year [month]", ylab = "")
for ( i in 2:6 )
  plotVarCoef(b.ngr[[i]], yvar = yvars[i], add = TRUE, col = cols[i], lwd = 2)

title(main=expression(f[2](season)*" - Innsbruck"), line=2, cex.main=1)
legend("topright", col = cols, lwd = 2, bty = "n",
       legend = paste("day", 1:6))
axis(2, las = 1)

ats <- c(as.POSIXlt(paste0("2011-",sprintf("%02d",1:12),"-01"),
                     tz="UTC")$yday+1, 366)
axis(1, at=ats, labels=FALSE)
ats <- as.POSIXlt(paste0("2011-",sprintf("%02d",1:12),"-15"), tz="UTC")$yday+1
labels <- c("J","F","M","A","M","J","J","A","S","O","N","D")
axis(1, at=ats[seq(1,11,2)], tick=FALSE, cex.axis=1, padj=-1, 
     labels=labels[seq(1,11,2)])
axis(1, at=ats[seq(2,12,2)], tick=FALSE, cex.axis=1, padj=-1, 
     labels=labels[seq(2,12,2)])
box()

