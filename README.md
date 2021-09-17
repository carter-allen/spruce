# spruce
## A suite of Bayesian multivariate finite mixture models for clustering single cell spatial transcriptomics data. 

The `spruce` package is a robust and comprehensive tool for analyzing single cell spatial transcriptomics data using Bayesian multivariate finite mixture models. We accommodate multiple gene expression distributions, including multivariate normal (MVN) and multivariate skew-normal (MSN). We also allow for a range of spatial assumptions that rely on spatially correlated CAR models for random intercepts and sptial prior smoothing in cluster indicators. 

`spruce` can also be used as a comprehensive data simulation tool for power analysis and the design of spatial transcriptomics experiments. 

## Installation

```
library(devtools)
install_github('carter-allen/spruce')
```

