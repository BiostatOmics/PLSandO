# PLSandO
PLSandO (PLS and Others) is a comprehensive R package developed for advanced multivariate modeling. PLSandO implements a wide range of methods including Principal Component Analysis (PCA), Partial Least Squares (PLS), and PLS Discriminant Analysis (PLS-DA), with extensive support for Multi-block (MB-PCA, MB-PLS, MB-PLSDA) and multilevel data structures. It is specifically designed to facilitate the exploration and interpretation of high-dimensional data, enabling researchers to move from raw data to interpretable models within a single framework.

## Installation

### Installing PLSandO

Currently, the package can be installed directly from GitHub using the `devtools` R package:

    install.packages("devtools")
    devtools::install_github("BiostatOmics/PLSandO")


Although all dependencies should install automatically, you can manually install them if the process fails:

* ggplot2
* naniar
* reshape2
* patchwork
* future
* furrr
* purrr
* kneedle
* biostatcolors

Note that the last two dependencies are available only on GitHub and must be installed using:

```r
devtools::install_github('BiostatOmics/biostatcolors')
devtools::install_github('etam4260/kneedle')
```

## Usage

You can find the User Guide for the package in the vignettes folder or access it directly [here](https://github.com/BiostatOmics/PLSandO/blob/master/vignettes/tutorial.html). 

