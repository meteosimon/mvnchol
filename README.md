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
    [**remotes**](https://remotes.r-lib.org/) or
    [**devtools**](https://devtools.r-lib.org/) in an R session:

<!-- -->

    remotes::install_git("https://git.uibk.ac.at/c4031039/mvnchol")

-   Or clone it using [git](https://git-scm.com/) and a git client for
    your system. And build and INSTALL it after cloning, e.g. in a bash:

<!-- -->

    git clone https://git.uibk.ac.at/c4031039/mvnchol.git
    R CMD build mvnchol
    R CMD INSTALL mvnchol_0.2.0.tar.gz

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

    ## This is mgcv 1.8-33. For overview type 'help("mgcv-package")'.

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

    ## AICc 9422.284 logPost -4839.99 logLik -4648.88 edf 60.349 eps 0.0000 iteration  11
    ## elapsed time:  4.15sec
    ## Starting the sampler...
    ## 
    ## |********************| 100%  0.00sec  1.26min

    plot(b, model = "lambda12")

![](README_files/figure-markdown_strict/unnamed-chunk-1-1.png)
