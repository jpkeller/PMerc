
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PMerc

Pariculate Matter (PM) Exposure-Response Curves (ERCs)

This R provides a lean set of code for computing specific values from an
exposure-response curve developed using the
[`bercs`](http://www.github.com/jpkeller/bercs) R package. This package
contains two elements:

  - **Data** objects containing the minimal set of information needed to
    calculate the full ERC or any desired values from it. Data files are
    provided for published ERC curves (see examples below).

  - **Functions**. Currently, there are four R functions included
    (including `compute_OR2()`), which are derived from functions of
    similar name from the `bercs` package (e.g. `compute_OR`). They are
    included here in `PMerc` to facilitate calculating ERC values
    without requiring the full `bercs` package, which requires compiled
    C code for the STAN model objects. This package is recommended only
    for calculating specific risk measures from the accompanying data
    files. If developing your own models, it is recommended to use the
    functions in the `bercs` package and not these. It is possible that
    these functions may be removed at a future date.

To install `PMerc`, use the following commands:

    devtools::install_github("jpkeller/PMerc")

## Example: PM-ALRI Curve Values

``` r
library(PMerc)
library(splines2)

data(nepal_pm_alri)

compute_OR2(
    expsequence = c(35, 37.5, 50, 75, 100, 150, 200, 400),
    ref_exposure=50,
    ciband=0.95,
    beta_post=nepal_pm_alri$posterior_params$beta,
    bs_post=nepal_pm_alri$posterior_params$bS,
    xdf=nepal_pm_alri$model_data$xdf,
    nS=nepal_pm_alri$model_data$S,
    Mx=nepal_pm_alri$model_data$Mx,
    Mx_attributes = nepal_pm_alri$model_data$Mx_attributes)
```

    ##   exposure logOR_mean   logOR_low logOR_high   OR_mean    OR_low  OR_high study
    ## 1     35.0 -0.8275990 -2.19826068  0.1925011 0.4370975 0.1109960 1.212278     1
    ## 2     37.5 -0.7724221 -2.03520848  0.1714126 0.4618929 0.1306532 1.186980     1
    ## 3     50.0  0.0000000  0.00000000  0.0000000 1.0000000 1.0000000 1.000000     1
    ## 4     75.0  0.6184297  0.05759951  1.1865671 1.8560112 1.0592907 3.275816     1
    ## 5    100.0  0.9871702  0.42800539  1.5537510 2.6836296 1.5341944 4.729176     1
    ## 6    150.0  1.2221348  0.63530381  1.8084508 3.3944265 1.8875955 6.100989     1
    ## 7    200.0  1.1964850  0.62168084  1.7741433 3.3084673 1.8620552 5.895229     1
    ## 8    400.0  0.9678931  0.39970055  1.5364891 2.6323924 1.4913780 4.648242     1

## Community guidelines

If you have a bug to report, are having technical issues, or want to
recommend features, please open a [Github
Issue](https://github.com/jpkeller/bercs/issues).
