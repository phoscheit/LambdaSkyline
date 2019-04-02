# LambdaSkyline

The `LambdaSkyline` package provides functions to manipulate multifurcating trees (as `ape` `phylo` objects) and
compute their skyline plot (originally introduced by Pybus et al.), using Lambda-coalescents (specifically, Beta(2-alpha,alpha)-coalescents) as model. The package also provides functions for simulating trees under the Beta-coalescent, using the results in (Schweinsberg 2003).

## Installation

To install `LambdaSkyline`, make sure the `devtools` library is installed, then type the following commands:

```
library(devtools)
devtools::install_github("phoscheit/LambdaSkyline")
```

## Usage

The main function provided by the package is `skyline.multi()`, which extends the `skyline` functions present in `ape` to enable the computation of skyline plots based on Lambda-coalescents, specifically Beta(2-alpha,alpha)-coalescents. The function takes a `phylo` tree as input, as well as the alpha parameter between 0 and 2, and returns the skyline plot as a `skyline` object (which can then be plotted using `ape`'s `plot.skyline` function). 

Strictly bifurcating trees have zero likelihood under any Beta(2-alpha,alpha)-coalescent model with 0 < alpha < 2. Such trees, especially containing branches with zero length, can be converted into multifurcating trees by the `di2multi.cons()` function, which collapses branches shorter than a given length. 

The function `betacoal.maxlik()` performs maximum-likelihood estimation of the alpha parameter on a given multifurcating tree. 