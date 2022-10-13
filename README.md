[**BAMLSS**](http://www.bamlss.org/) families for [multivariate Gaussion
regression](https://arxiv.org/abs/2102.13518) models. The multivariate
Gaussian distribution is reparameterized using the (modified) Cholesky
decomposition of the covariance matrix.

Installation
------------

This package is hosted on the gitlab server of the Universitaet
Innsbruck. The repository has public access.

-   You can either install it directly on your system by using
    [**remotes**](https://remotes.r-lib.org/) or
    [**devtools**](https://devtools.r-lib.org/) in an R session:

<!-- -->

    remotes::install_git("https://git.uibk.ac.at/c4031039/mvnchol")

-   Or clone it using [git](https://git-scm.com/) and a git client for
    your system. And build and INSTALL it after cloning, e.g. in a bash:

<!-- -->

    git clone https://git.uibk.ac.at/c4031039/mvnchol.git
    R CMD build mvnchol
    R CMD INSTALL mvnchol_0.3.0.tar.gz

-   Or download the source code as tarball
    [mvnchol-master.tar.gz](https://git.uibk.ac.at/c4031039/mvnchol/-/archive/master/mvnchol-master.tar.gz).
    Here, you also have to build and INSTALL the package after the
    download.

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
    ## elapsed time:  3.65sec
    ## Starting the sampler...
    ## 
    ## |                    |   0%  1.05min
    ## |*                   |   5%  1.05min  3.31sec
    ## |**                  |  10% 57.50sec  6.39sec
    ## |***                 |  15% 54.04sec  9.54sec
    ## |****                |  20% 52.95sec 13.24sec
    ## |*****               |  25% 51.24sec 17.08sec
    ## |******              |  30% 48.83sec 20.93sec
    ## |*******             |  35% 46.00sec 24.77sec
    ## |********            |  40% 42.76sec 28.51sec
    ## |*********           |  45% 39.50sec 32.31sec
    ## |**********          |  50% 36.19sec 36.19sec
    ## |***********         |  55% 32.73sec 40.00sec
    ## |************        |  60% 29.29sec 43.93sec
    ## |*************       |  65% 25.71sec 47.74sec
    ## |**************      |  70% 22.16sec 51.71sec
    ## |***************     |  75% 18.60sec 55.80sec
    ## |****************    |  80% 14.99sec 59.96sec
    ## |*****************   |  85% 11.29sec  1.07min
    ## |******************  |  90%  7.53sec  1.13min
    ## |******************* |  95%  3.77sec  1.19min
    ## |********************| 100%  0.00sec  1.26min

    plot(b, model = "lambda12")

![](README_files/figure-markdown_strict/unnamed-chunk-1-1.png)
