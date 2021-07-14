
<!-- README.md is generated from README.Rmd. Please edit that file -->

# VIVID: Visualisation of Variable Importance Differences

The VIVID (Variability of Variable Importance Differences) package
implements both a feature selection method and visualization for complex
data.

## Installation

``` r
library(devtools)
devtools::install_github("SmithConnor/VividR")
install.packages('ropls')
install.packages('furrr')
```

## A quick example

``` r
library('VIVID')
library('ropls')

data("sacurine") #Load sacurine dataset from the 'ropls' package

dat <- sacurine$dataMatrix
outcomes <- sacurine$sampleMetadata$gender

vivid.sacurine <- vivid(x = dat,
                        y = outcomes,
                        bootstraps = 50,
                        cores = parallel::detectCores() - 1,
                        seed = 1234567,
                        lambda = 'lambda.1se',
                        compareMethod = 'BIC')
```

From the feature selection algorithm using the BIC we have selected:

``` r
vivid.sacurine$optFeatures
#> [1] "X.gamma.Glu.Leu.Ile"          "Glu.Val"                     
#> [3] "Gluconic.acid.and.or.isomers" "Malic.acid"                  
#> [5] "N2.Acetylaminoadipic.acid"    "Oxoglutaric.acid"            
#> [7] "p.Hydroxyhippuric.acid"       "Pantothenic.acid"            
#> [9] "Testosterone.glucuronide"
```
