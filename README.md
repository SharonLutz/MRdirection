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

E\[Xtrue \] = $\gamma_0$ + $GX\gamma_{GX}$ + $GY\gamma_{GY}$

where if $\gamma_{GY}$ is nonzero then the SNPs GY are pleiotropic.

All of these values are inputted by the user (i.e. the intercept gamma0, and the vectors of genetic effect sizes gammaGX and gammaGY). If there is no measurement error (input measurementError=F), then X=Xtrue. If there is measurement error (input measurementError=T), then the measured phenotype 1 X is generated from the true phenotype 1 Xtrue such that:

E\[Xtrue \] = $\delta_0$ + $\delta_{X}X_{true}$

where $\delta_0$ and $\delta_{X}$ are inputted by the user (input delta0, deltaX). Phenotype 2 Y is generated from a normal distribution with the variance (input varY) and the mean as follows:

E\[Y1 \] = $\beta_0$ + $\beta_{X}X_{true}$ + $GY\beta_{GY}$ + $GX\beta_{GX}$

where if $\beta_{GX}$ is nonzero then the SNPs GX are pleiotropic.

All of these values are inputted by the user (i.e. the intercept beta0, the effect of phenotype 1 Xtrue on phenotype 2 as betaX, the vector of the effect of the SNPs GX directly on phenotype 2 as betaGX, and the vector of the effect of the SNPs GY directly on phenotype 2 as betaGY).

If there is unmeasured confounding (unmeasuredConfounding=T) between the exposure X and the outcome Y, then the unmeasured confounder U is generated from a normal distribution with user specified mean and variance (i.e. meanU, varU). Then, the exposure X and outcome Y are generated such that

E\[Xtrue \] = $\gamma_0$ + $GX\gamma_{GX}$ + $GY\gamma_{GY}$ + $\gamma_{U}U$

E\[Y1 \] = $\beta_0$ + $\beta_{X}X_{true}$ + $GY\beta_{GY}$ + $GX\beta_{GX}$ + $\beta_{U}U$





## Output

## Example


## References
