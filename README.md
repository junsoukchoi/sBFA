sBFA: Sparse Bayesian Factor Models with the Dirichlet-Laplace Priors
================
Junsouk Choi

- <a href="#installation" id="toc-installation">Installation</a>
- <a href="#usage" id="toc-usage">Usage</a>

The R package `sBFA` implements an efficient data-augmented Gibbs
sampler for the sparse Bayesian factor models proposed by Pati et al.,
(2014) that employ the Dirichlet-Laplace priors. The implementation of
this sampler aligns with the approach detailed by Bhattacharya et
al. (2015).

## Installation

To install the latest version from Github, use

``` r
library(devtools)
devtools::install_github("junsoukchoi/sBFA")
```

## Usage

The following example describes how to use the package `sBFA` to do the
posterior inference for the sparse Bayesian factor models with the
Dirichlet-Laplace priors. The function `sBFA_DL` takes data with the
number of factor to be fitted and some arguments needed for the Gibbs
sampler, and returns MCMC samples from the posterior distribution of the
sparse Bayesian factor models with the Dirichlet-Laplace priors.

``` r
library(sBFA)
set.seed(7)

# set the sample size n, dimension p, and the number of factors q
n = 100
p = 100
q = floor(log(p))

# generate true Lambda
# s = log(p) nonzero elements are drawn uniformly between 1 and 2 per columns
Lambda = matrix(0, p, q)
for (k in 1 : q)
{
   id_nonzero = sample(1 : p, floor(log(p)))
   Lambda[id_nonzero, k] = runif(floor(log(p)), 1, 2)
}

# set true Sigma to be identity
Sigma = diag(1, p)

# generate data from a sparse factor model with the given Lambda and Sigma
Y = matrix(NA, n, p)
U = matrix(NA, n, q)
for (i in 1 : n)
{
   U[i, ] = rnorm(q)
   Y[i, ] = Lambda %*% U[i, ] + sqrt(Sigma) %*% rnorm(p)
}

# choose the values of hyperparameters 
priors = list()
priors$Sigma = c(0.1, 0.1)
priors$Phi   = 0.5

# obtain starting values for MCMC from the prior distribution
starting = list()
starting$Phi    = matrix(1, p, q)
starting$Phi    = starting$Phi / sum(starting$Phi)
starting$tau    = rgamma(1, shape = p * q * priors$Phi, rate = 0.5)
starting$Psi    = matrix(rexp(p * q, rate = 0.5), p, q)
starting$Lambda = matrix(rnorm(p * q, sd = c(sqrt(starting$Psi * starting$tau^2 * starting$Phi^2))), p, q)
starting$U      = matrix(rnorm(n * q), n, q)
starting$Sigma  = diag(1, p)

# run MCMC for the sparse Bayesian factor models with the Dirichlet-Laplace priors
out = sBFA_DL(Y, q, starting, priors)

# calculate the posterior mean of the loading matrix
Lambda_est = apply(out$Lambda, c(1, 2), mean)

# visualize the true loading matrix and the posterior mean of the loading matrix
library(gplots)
library(RColorBrewer)
Colors=rev(brewer.pal(7, "RdBu"))
Colors=colorRampPalette(Colors)(100)
breaks = seq(-2.5, 2.5, length.out = 101)
pdf(file = "Lambda_true.pdf", width = 10, height = 8)
heatmap.2(Lambda, col = Colors, density.info = "none", dendrogram = 'none', Colv = FALSE, Rowv =FALSE, trace = "none", breaks = breaks, 
          keysize = 1, key.par = list(mar=c(3.5,0,3,4)), lmat=rbind(c(5,4,2),c(6,1,3)), lhei=c(1,5), lwid=c(1,10,1))
dev.off()
pdf(file = "Lambda_est.pdf", width = 10, height = 8)
heatmap.2(Lambda_est, col = Colors, density.info = "none", dendrogram = 'none', Colv = FALSE, Rowv =FALSE, trace = "none", breaks = breaks, 
          keysize = 1, key.par = list(mar=c(3.5,0,3,4)), lmat=rbind(c(5,4,2),c(6,1,3)), lhei=c(1,5), lwid=c(1,10,1))
dev.off()
```
