dT <- .1
T <- seq(-20, 40, by=dT)
T0 <- 20
demandCurve <- (T - T0)^2 / exp((2/3) - (1/45)*T) 
demandCurve <- demandCurve / max(demandCurve)
plot(T, demandCurve, type="l")
predictive <- dnorm(T, mean=-5, sd=4)
lines(T, predictive, col=2)
energyDemand <- sum(predictive * demandCurve) * dT
