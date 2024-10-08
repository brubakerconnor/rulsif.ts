
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rulsif.ts

<!-- badges: start -->

<!-- badges: end -->

`rulsif.ts` implements relative unconstrained least squares importance
fitting for the task of detecting change points in time series data.
Below is a simple example of the use of this package.

``` r
series <- c(
    rnorm(50, mean = 0, sd = 0.3),
    rnorm(25, mean = 8, sd = 1),
    rnorm(75, mean = 3, sd = 0.6),
    rnorm(25, mean = 1, sd = 0.8),
    rnorm(100, mean = -5, sd = 1.5),
    rnorm(100, mean = -5, sd = 0.2),
    rnorm(50, mean = -2.5, sd = 0.4),
    rnorm(50, mean = 2, sd = 1.2)
)
d <- ts_detect(series, window_size = 3, step = 10, make_plot = TRUE)
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

## Installation

You can install the the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("brubakerconnor/rulsif.ts")
```

## References

  - Song Liu, et al. “Change-point detection in time-series data by
    relative density-ratio estimation”. In: *Neural Networks* 43 (2013),
    pp. 72-83. ISSN: 0893-6080. DOI: `10.1016/j.neunet.2013.01.012.`

  - Masashi, et al. “Direct importance estimation for covariance shift
    adaptation.” In: *Annals of the Institute of Statistical
    Mathematics* 60.4 (2008), pp. 699-746. ISSN: 0020-3157. DOI:
    `10.1007/s10463-008-0197-x.`
