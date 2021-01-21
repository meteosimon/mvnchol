bamlssMVN
=========

[**BAMLSS**](http://www.bamlss.org/) families for modelling multivariate
normal distributions reparameterized using the (modified) Cholesky
decomposition of the covariance matrix.

Installation
------------

This package is hosted on the gitlab server of the Universitaet
Innsbruck. The repository has public access.

-   You can either install it directly on your system by using
    [**devtools**](https://devtools.r-lib.org/) in an R session:

<!-- -->

    remotes::install_git("https://git.uibk.ac.at/c4031039/bamlss_mvn")

-   Or clone it using [git](https://git-scm.com/) and a git client for
    your system. And build and INSTALL it after cloning, e.g. in a bash:

<!-- -->

    git clone https://git.uibk.ac.at/c4031039/bamlss_mvn.git
    R CMD build bamlss_mvn
    R CMD INSTALL bamlssMVN_0.1.3.tar.gz

-   Or download the source code as tarball
    [bamlss\_mvn-master.tar.gz](https://git.uibk.ac.at/c4031039/bamlss_mvn/-/archive/master/bamlss_mvn-master.tar.gz).
    Here, you also have to build and INSTALL the package after the
    download.

First Model
-----------

    library(bamlss)

    ## Loading required package: coda

    ## Loading required package: colorspace

    ## Loading required package: mgcv

    ## Loading required package: nlme

    ## This is mgcv 1.8-33. For overview type 'help("mgcv-package")'.

    ## 
    ## Attaching package: 'bamlss'

    ## The following object is masked from 'package:mgcv':
    ## 
    ##     smooth.construct

    library(bamlssMVN)
    data(simdata)
    d <- simdata$d
    f <- make_formula(y1 | y2 | y3 ~ s(x) | s(x) | s(x))
    b <- bamlss(f, family = mvnchol_bamlss(k = 3), data = d)

    ## AICc 9796.373 logPost -5012.88 logLik -4850.79 edf 46.275 eps 1.0000 iteration   1
    ## AICc 9472.640 logPost -4863.28 logLik -4674.78 edf 59.666 eps 1.0557 iteration   2
    ## AICc 9433.197 logPost -4844.09 logLik -4655.19 edf 59.545 eps 0.2736 iteration   3
    ## AICc 9424.514 logPost -4840.35 logLik -4650.41 edf 59.956 eps 0.3922 iteration   4
    ## AICc 9422.810 logPost -4839.90 logLik -4649.28 edf 60.214 eps 0.1256 iteration   5
    ## AICc 9422.444 logPost -4839.92 logLik -4649.00 edf 60.308 eps 0.0353 iteration   6
    ## AICc 9422.355 logPost -4839.98 logLik -4648.92 edf 60.344 eps 0.0091 iteration   7
    ## AICc 9422.301 logPost -4839.99 logLik -4648.89 edf 60.344 eps 0.0030 iteration   8
    ## AICc 9422.290 logPost -4839.99 logLik -4648.88 edf 60.347 eps 0.0008 iteration   9
    ## AICc 9422.286 logPost -4839.99 logLik -4648.88 edf 60.348 eps 0.0002 iteration  10
    ## AICc 9422.284 logPost -4839.99 logLik -4648.88 edf 60.349 eps 0.0000 iteration  11
    ## AICc 9422.284 logPost -4839.99 logLik -4648.88 edf 60.349 eps 0.0000 iteration  11
    ## elapsed time:  4.30sec
    ## Starting the sampler...
    ## 
    ## |                    |   0%  1.21min
    ## |*                   |   5%  1.08min  3.42sec
    ## |**                  |  10%  1.01min  6.71sec
    ## |***                 |  15% 56.00sec  9.88sec
    ## |****                |  20% 54.37sec 13.59sec
    ## |*****               |  25% 52.59sec 17.53sec
    ## |******              |  30% 50.14sec 21.49sec
    ## |*******             |  35% 47.33sec 25.49sec
    ## |********            |  40% 44.25sec 29.50sec
    ## |*********           |  45% 40.77sec 33.36sec
    ## |**********          |  50% 37.29sec 37.29sec
    ## |***********         |  55% 33.72sec 41.21sec
    ## |************        |  60% 30.10sec 45.16sec
    ## |*************       |  65% 26.46sec 49.13sec
    ## |**************      |  70% 22.76sec 53.10sec
    ## |***************     |  75% 19.03sec 57.10sec
    ## |****************    |  80% 15.24sec  1.02min
    ## |*****************   |  85% 11.44sec  1.08min
    ## |******************  |  90%  7.64sec  1.15min
    ## |******************* |  95%  3.83sec  1.21min
    ## |********************| 100%  0.00sec  1.28min

    plot(b, model = "lambda12")

![](README_files/figure-markdown_strict/unnamed-chunk-1-1.png)
