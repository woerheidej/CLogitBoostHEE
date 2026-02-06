# CLogitBoostHEE

CLogitBoostHEE is a package that provides machine learning methods for analyzing heterogeneous exposure effects in matched case-control studies.

For further information on CLogitBoostHEE see: "Missing"

# Installation

The package can be installed in R via

``` r
install.packages("devtools")
devtools::install_github("woerheidej/CLogitBoostingHEE")
```

# Example CLogitBoostHEE

``` r
library(CLogitBoostingHEE)
data(sim_data)

set.seed(1899)
sim_results <- CLogitBoostingHEE(
  stroke_data,
  exposure = "X",
  strata = "strata",
  outcome = "y",
  matching = c("s", "d"),
  q = 5,
  B = 50,
  PFER = 0.5,
  nu = 1,
  mstop = 2000,
  sampling_type = "SS",
  assumption = "r-concave",
  n_cores = 1
))

plot(sim_results)
```
vi <- varimp(illu.rf)
plot(vi)
```
