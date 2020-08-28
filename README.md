
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PMerc

Pariculate Matter (PM) Exposure-Response Curves (ERCs)

This R provides a lean set of code for computing specific values from an
exposure-response curve developed using the
[`bercs`](www.github.com/jpkeller/bercs) R package. This package
contains two elements:

  - **Data** objects containing the minimal set of information needed to
    calculate the full ERC or any desired values from it. Data files are
    provided for published ERC curves (see examples below).

  - **Functions**. Currently, there are three R functions included,
    which are exact copies of functions of the same name from the
    `bercs` package. They are included here in `PMerc` to facilitate
    calculating ERC values without requiring the full `bercs` package
    (which requires compiled C code for the STAN model objects). This
    package is recommended only for calculating specific risk measures
    from the accompanying data files. If developing your own models, it
    is recommended to use the functions in the `bercs` package and not
    these. It is possible that these functions may be removed at a
    future date.

To install `PMerc`, use the following commands:

    devtools::install_github("jpkeller/PMerc")

## Example: PM-ALRI Curve Values

``` r
library(PMerc)
library(splines2)

data(nepal_pm_alri)

erc <- compute_ERC2(
                      expsequence = c(35, 37, 50, 75, 100, 150, 200),
                      ref_exposure=50,
                      inclInterceptUncertainty=T,
                      ciband=0.95,
                      beta_post=nepal_pm_alri$posterior_params$beta,
                      bs_post=nepal_pm_alri$posterior_params$bS,
                      xdf=nepal_pm_alri$model_data$xdf,
                      nS=nepal_pm_alri$model_data$S,
                      Mx=nepal_pm_alri$model_data$Mx,
                      intercept_prop="equal",
                      Mx_attributes = nepal_pm_alri$model_data$Mx_attributes)
erc
```

    ##         mean        low      high exposure
    ## 1 -0.8275990 -2.3284909 0.4183978       35
    ## 2 -0.7911737 -2.2202256 0.4104241       37
    ## 3  0.0000000 -0.5675226 0.5559417       50
    ## 4  0.6184297  0.1519882 1.0675907       75
    ## 5  0.9871702  0.5811591 1.3872459      100
    ## 6  1.2221348  0.8169162 1.6326992      150
    ## 7  1.1964850  0.8451832 1.5635276      200

These values are on the log-odds scale. To compute odds ratios,
exponentiate:

``` r
library(dplyr)
erc %>%
    mutate(OR=exp(mean),
           ORlow=exp(low),
           ORhigh=exp(high))
```

    ##         mean        low      high exposure        OR      ORlow   ORhigh
    ## 1 -0.8275990 -2.3284909 0.4183978       35 0.4370975 0.09744269 1.519525
    ## 2 -0.7911737 -2.2202256 0.4104241       37 0.4533124 0.10858461 1.507457
    ## 3  0.0000000 -0.5675226 0.5559417       50 1.0000000 0.56692821 1.743582
    ## 4  0.6184297  0.1519882 1.0675907       75 1.8560112 1.16414656 2.908364
    ## 5  0.9871702  0.5811591 1.3872459      100 2.6836296 1.78810985 4.003808
    ## 6  1.2221348  0.8169162 1.6326992      150 3.3944265 2.26350881 5.117670
    ## 7  1.1964850  0.8451832 1.5635276      200 3.3084673 2.32840440 4.775638

## Community guidelines

If you have a bug to report, are having technical issues, or want to
recommend features, please open a [Github
Issue](https://github.com/jpkeller/bercs/issues).
