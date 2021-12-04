# spruce <img src="inst/logo/spruce_transp.png" align="right" width="115" />
## A suite of Bayesian multivariate finite mixture models for clustering single cell spatial transcriptomics data. 

The `spruce` package is a robust and comprehensive tool for analyzing single cell spatial transcriptomics data using Bayesian multivariate finite mixture models. `spruce` accommodates multiple gene expression distributions, namely multivariate normal (MVN) and multivariate skew-normal (MSN). The MVN model, which is simpler and more computationally efficient, is suitable for modeling dimension reduction features as provided from approaches like principal components analysis, while the MSN model should be used when modeling normalized gene expression features directly, as conversion of over-dispersed count data to continuous features results in inherently right-skewed features. `spruce` also allows for a range of methods for accomodating spatial correlation across a tissue sample, including spatially correlated CAR/MCAR random intercepts or spatial prior smoothing in cluster indicators via a Potts similar to [BayesSpace](https://www.nature.com/articles/s41587-021-00935-2). 

`spruce` can also be used as a comprehensive data simulation tool for power analysis and the design of spatial transcriptomics experiments. 

While `spruce` is intended for modeling single-sample HST experiments, we have developed [maple](https://github.com/carter-allen/maple) for extension of this methodology to multi-sample HST data.

## Installation

```
library(devtools)
install_github('carter-allen/spruce')
```

## Vignettes

