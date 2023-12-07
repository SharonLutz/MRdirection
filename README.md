# MRdirection
Examines MR approaches to detect the directionality between phenotype 1 X and phenotype 2 Y through simulation studies.

## Installation
```
install.packages("devtools")
devtools::install_github("xue-hr/MRCD")
devtools::install_github("MRCIEU/TwoSampleMR")
devtools::install_github("xue-hr/MRcML")
devtools::install_github("xue-hr/BiDirectCausal")

devtools::install_github("SharonLutz/MRdirection")
```

## Input
First, the SNPs, GX and GY, are generated from a binomial distribution for max(nX, nY) subjects (input nX, nY) for a given vector of minor allele frequencies (input MAF_GX and MAF_GY).
The true phenotype 1 (Xtrue) is generated from a normal distribution with the variance (input varX) and the mean as follows:

E\[Xtrue \] = $\gamma_0$ + GX$\gamma_GX$ + GY$\gamma_GY$

All of these values are inputted by the user (i.e. the intercept gamma0, and the vectors of genetic effect sizes gammaGX and gammaGY). If there is no measurement error (input measurementError=F), then X=Xtrue. If there is measurement error (input measurementError=T), then the measured phenotype 1 X is generated from the true phenotype 1 Xtrue such that:

E\[Xtrue \] = $\delta_0$ + $\delta_{X}X_{true}$





## Output

## Example


## References
