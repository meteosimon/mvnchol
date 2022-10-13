[**BAMLSS**](http://www.bamlss.org/) families for [multivariate Gaussion
regression](https://doi.org/10.1016/j.ecosta.2022.03.001) models. The
multivariate Gaussian distribution is reparameterized using the
(modified) Cholesky decomposition of the covariance matrix.

Installation
------------

You can either install it directly on your system by using
[**remotes**](https://remotes.r-lib.org/) or
[**devtools**](https://devtools.r-lib.org/) in an R session:

    remotes::install_github("meteosimon/mvnchol")

Or clone it using [git](https://git-scm.com/) and a git client for your
system. And build and INSTALL it after cloning, e.g. in a bash:

    git clone git@github.com:meteosimon/mvnchol.git
    R CMD build mvnchol
    R CMD INSTALL mvnchol_0.3.0.tar.gz

First Model
-----------

    library(bamlss)

    ## Loading required package: coda

    ## Loading required package: colorspace

    ## Loading required package: mgcv

    ## Loading required package: nlme

    ## This is mgcv 1.8-40. For overview type 'help("mgcv-package")'.

    ## 
    ## Attaching package: 'bamlss'

    ## The following object is masked from 'package:mgcv':
    ## 
    ##     smooth.construct

    library(mvnchol)
    data(simdata)
    d <- simdata$d
    f <- make_formula(y1 | y2 | y3 ~ s(x) | s(x) | s(x))
    b <- bamlss(f, family = mvnchol_bamlss(k = 3), data = d)

    ## AICc 9796.373 logPost -5012.88 logLik -4850.79 edf 46.275 eps 1.0000 iteration   1
    ## AICc 9471.761 logPost -4863.12 logLik -4674.34 edf 59.666 eps 1.0564 iteration   2
    ## AICc 9432.365 logPost -4843.92 logLik -4654.78 edf 59.539 eps 0.2204 iteration   3
    ## AICc 9423.675 logPost -4840.18 logLik -4650.00 edf 59.947 eps 0.2975 iteration   4
    ## AICc 9421.974 logPost -4839.73 logLik -4648.88 edf 60.206 eps 0.0972 iteration   5
    ## AICc 9421.608 logPost -4839.75 logLik -4648.59 edf 60.300 eps 0.1687 iteration   6
    ## AICc 9421.519 logPost -4839.81 logLik -4648.51 edf 60.337 eps 0.0098 iteration   7
    ## AICc 9421.465 logPost -4839.82 logLik -4648.48 edf 60.336 eps 0.0027 iteration   8
    ## AICc 9421.456 logPost -4839.82 logLik -4648.47 edf 60.340 eps 0.0008 iteration   9
    ## AICc 9421.452 logPost -4839.83 logLik -4648.47 edf 60.341 eps 0.0002 iteration  10
    ## AICc 9421.450 logPost -4839.83 logLik -4648.47 edf 60.342 eps 0.0000 iteration  11
    ## AICc 9421.450 logPost -4839.83 logLik -4648.47 edf 60.342 eps 0.0000 iteration  11
    ## elapsed time:  3.75sec
    ## Starting the sampler...
    ## 
    ## |                    |   0%  1.14min
    ## |*                   |   5%  1.04min  3.27sec
    ## |**                  |  10% 58.11sec  6.46sec
    ## |***                 |  15% 54.07sec  9.54sec
    ## |****                |  20% 52.46sec 13.11sec
    ## |*****               |  25% 50.73sec 16.91sec
    ## |******              |  30% 48.31sec 20.70sec
    ## |*******             |  35% 45.56sec 24.53sec
    ## |********            |  40% 42.54sec 28.36sec
    ## |*********           |  45% 39.33sec 32.18sec
    ## |**********          |  50% 36.00sec 36.00sec
    ## |***********         |  55% 32.49sec 39.71sec
    ## |************        |  60% 28.98sec 43.48sec
    ## |*************       |  65% 25.43sec 47.23sec
    ## |**************      |  70% 21.93sec 51.17sec
    ## |***************     |  75% 18.46sec 55.38sec
    ## |****************    |  80% 14.88sec 59.52sec
    ## |*****************   |  85% 11.19sec  1.06min
    ## |******************  |  90%  7.48sec  1.12min
    ## |******************* |  95%  3.74sec  1.19min
    ## |********************| 100%  0.00sec  1.25min

    plot(b, model = "lambda12")

![](README_files/figure-markdown_strict/unnamed-chunk-1-1.png)
