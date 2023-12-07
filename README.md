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
First, the SNPs, $G_X$ and $G_Y$, are generated from a binomial distribution for max(nX, nY) subjects (input nX, nY) for a given vector of minor allele frequencies (input MAF_GX and MAF_GY).
The true phenotype 1 (Xtrue) is generated from a normal distribution with the variance (input varX) and the mean as follows:

E\[Xtrue \] = $\gamma_0$ + $G_X*\gamma_{G_X}$ + $G_Y*\gamma_{G_Y}$

where if $\gamma_{G_Y}$ is nonzero then the SNPs $G_Y$ are pleiotropic.

All of these values are inputted by the user (i.e. the intercept gamma0, and the vectors of genetic effect sizes gammaGX and gammaGY). If there is no measurement error (input measurementError=F), then $X=X_true$. If there is measurement error (input measurementError=T), then the measured phenotype 1 X is generated from the true phenotype 1 $X_true$ such that:

E\[Xtrue \] = $\delta_0$ + $\delta_{X}*X_{true}$

where $\delta_0$ and $\delta_{X}$ are inputted by the user (input delta0, deltaX). Phenotype 2 Y is generated from a normal distribution with the variance (input varY) and the mean as follows:

E\[Y1 \] = $\beta_0$ + $\beta_{X}X_{true}$ + $G_Y\beta_{G_Y}$ + $GX\beta_{G_X}$

where if $\beta_{G_X}$ is nonzero then the SNPs $\G_X$ are pleiotropic.

All of these values are inputted by the user (i.e. the intercept beta0, the effect of phenotype 1 Xtrue on phenotype 2 as betaX, the vector of the effect of the SNPs $G_X$ directly on phenotype 2 as betaGX, and the vector of the effect of the SNPs $G_Y$ directly on phenotype 2 as betaGY).

If there is unmeasured confounding (unmeasuredConfounding=T) between the exposure X and the outcome Y, then the unmeasured confounder U is generated from a normal distribution with user specified mean and variance (i.e. meanU, varU). Then, the exposure X and outcome Y are generated such that

E\[Xtrue \] = $\gamma_0$ + $G_X\gamma_{G_X}$ + $G_Y\gamma_{G_Y}$ + $\gamma_{U}U$

E\[Y1 \] = $\beta_0$ + $\beta_{X}X_{true}$ + $G_Y\beta_{G_Y}$ + $G_X\beta_{G_X}$ + $\beta_{U}U$





## Output

## Example


## References
