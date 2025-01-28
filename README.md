# hurdleNB
[![R-CMD-check](https://github.com/muskuloes/hurdleNB/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/muskuloes/hurdleNB/actions/workflows/R-CMD-check.yaml)

Hurdle negative binomial family for [`mgcv`](https://cran.r-project.org/web/packages/mgcv/index.html).

### Install
```r
# using remotes
# install.packages("remotes")
remotes::install_github("muskuloes/hurdleNB")
```

### Example

```r
library(mgcv)
library(hurdleNB)

set.seed(1)
n <- 400
dat <- gamSim(1, n = n)
dat$y <- rhurdleNB(dat$f / 4 - 1)

m <- gam(y ~ s(x0) + s(x1) + s(x2) + s(x3),
  data = dat, family = hurdleNB()
)
```
