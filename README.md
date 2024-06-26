# LatentSNA tutorial
## Overview


The LatentSNA model is specifically designed for identifying imaging biomarkers using brain connectivity data. It offers robust tools for unbiased estimation of imaging biomarkers' influence on behavior variants, quantification of the uncertainty and evaluation of the likelihood of the estimated biomarker effects against chance, brain-behavior prediction in novel samples for both connectivity and behaviors.



## Installation

This toolkit is implemented in R. Follow these steps for setup:


1. Clone or download the repository to your local machine.
2. Open R and navigate to the directory containing the toolkit.
3. Run the latentSNA model.



## Required Data


The toolkit is designed to analyze data consisting of brain connectivity and individual outcomes. Specifically, it requires two key files:

* `X`: a list of $V \times V$ brain connectivity data.
* `Y`: a matrix of $N \times P$ individual outcome data.

Additionally, simulated example data is available in the directory `data/X.RData` and `data/Y.RData` for demonstration purposes.


## Key parameters

* `W` a matrix of $N \times Q$ covariates for the connectivity data.
* `H` a matrix of $N \times Q1$ covariates for the attribute data.
* `nscan` number of iterations of the Markov chain (beyond burn-in)
* `burn` burn-in for the Markov chain

Note that sufficient burn-in is need to reach optimal covariance parameter estimates. See details in the method paper. 

## Usage


The main functionality of the LatentSNA Toolkit is encapsulated in the `latentSNA.r` script, which performs MCMC estimation of each of the unknown quantities. Toy Example data are in the `data` folder. The `latentSNA.r` can do the following connectivity data analysis:

### Run LatentSNA model

```{r s1, eval=FALSE}
library(latentSNA)

setwd("/Users/selena/Desktop/github/LatentSNA/data")

load(file='X.rda')
load(file='Y.rda')


model1=latentSNA(X, Y,W=NULL, H=NULL,
                   seed = 1, nscan = 10, burn = 1, odens = 1,
                   prior=list())

```


### Obtain unbiased estimates of covariance parameters.

From the saved model results, we can obtain estimated $N \times V$ brain latent connectivity as `model1$UPM`, posterior samples of the covariance parameters can be obtained via `model1$COV`, etc. See details in the pacakge documentation. 

