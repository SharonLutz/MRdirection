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
First, the SNPs ($G_X$ and $G_Y$) are generated from a binomial distribution for the max of nX and nY subjects (input nX, nY) for a given vector of minor allele frequencies (input MAF_GX and MAF_GY).
The true phenotype 1 (Xtrue) is generated from a normal distribution with the variance (input varX) and the mean as follows:

E\[Xtrue \] = $\gamma_0$ + $G_X*\gamma_{G_X}$ + $G_Y*\gamma_{G_Y}$

where if $\gamma_{G_Y}$ is nonzero then the SNPs $G_Y$ are pleiotropic.

All of these values are inputted by the user (i.e., the intercept gamma0, and the vectors of genetic effect sizes gammaGX and gammaGY). If there is no measurement error (input measurementError=F), then $X=Xtrue$. If there is measurement error (input measurementError=T), then the measured phenotype 1 X is generated from the true phenotype 1 $X_{true}$ such that:

E\[Xtrue \] = $\delta_0$ + $\delta_{X}X_{true}$

where $\delta_0$ and $\delta_{X}$ are inputted by the user (input delta0, deltaX). Phenotype 2 Y is generated from a normal distribution with the variance (input varY) and the mean as follows:

E\[Y1 \] = $\beta_0$ + $\beta_{X}X_{true}$ + $G_Y\beta_{G_Y}$ + $G_X\beta_{G_X}$

where if $\beta_{G_X}$ is nonzero then the SNPs $G_X$ are pleiotropic.

All of these values are inputted by the user (i.e. the intercept beta0, the effect of phenotype 1 Xtrue on phenotype 2 as betaX, the vector of the effect of the SNPs $G_X$ directly on phenotype 2 as betaGX, and the vector of the effect of the SNPs $G_Y$ directly on phenotype 2 as betaGY).

If there is unmeasured confounding (unmeasuredConfounding=T) between the exposure X and the outcome Y, then the unmeasured confounder U is generated from a normal distribution with user specified mean and variance (i.e. meanU, varU). Then, the exposure X and outcome Y are generated such that

E\[Xtrue \] = $\gamma_0$ + $G_X\gamma_{G_X}$ + $G_Y\gamma_{G_Y}$ + $\gamma_{U}U$

E\[Y1 \] = $\beta_0$ + $\beta_{X}X_{true}$ + $G_Y\beta_{G_Y}$ + $G_X\beta_{G_X}$ + $\beta_{U}U$

After the data are generated, then the MRdirection function runs the specified approaches to determine if the measured exposure X causes the outcome Y. 

The approaches can be specified as the following (input: runMethods):
"All": all approaches will be run.

"MRSteigerMethods": runs the MR Steiger approach.

"CDMethods": runs the CD approaches.

"BidirectionalMethods": runs the twoSampleMR approaches.

The following approaches can be individually specified in a list as input for runMethods: "MRS.ivw", "MRS.wMedian", "MRS.Egger", "CDRatio", "CDEgger", "CDgls", "tsmr_IVW", "tsmr_Egger", "tsmr_weighted_median", "tsmr_uwr", "tsmr_simple_median", "tsmr_pwm", "tsmr_meta_fixed_simple", "tsmr_IVW_mre", "tsmr_IVW_fe", "tsmr_Egger_Boot", "BDCDcML.S.DP", "BDMRcML.S.DP", "BDCD.Ratio.S", "BDCD.Egger.S"


## Output
This function outputs matrices of the proportion of simulations where Case 1 (X->Y), Case 2 (Y->X), and Case 3 (inconclusive) are returned for each of the specified approaches. The matrices are also saved to the working directory.

E\[Y1 \] = $\beta_0$ + $\beta_{X}X_{true}$ + $G_Y\beta_{G_Y}$ + $G_X\beta_{G_X}$

## Example
Consider an example with 1000 subjects for both X and Y (input nX=1000 and nY=1000) with a MAF_GX of 50 (input MAF_GX=0.5) and MAF_GY of 50 (input MAF_GY=0.5). Consider no pleiotropy, measurement error, or unmeasured counfouding (input measurementError = F, unmeasuredConfounding=F, betaGX=0, gammaGY=0). Then, let X be generated from a normal distribution with a variance of 1 (input varX = 1) and mean such that E\[Xtrue \] = 0 + 0.2\*GX (input gamma0=0, gammaGX=0.2). 
Y1 is generated from a normal distribution with a variance of 0.2 (input varY = 0.2) and mean such that E\[Y1 \] = 0 + $\beta_{X}X$ + 0.2\*GY (input beta0 = 0, betaGY=0.2) and beta_{X} varies from 0 to 2 (betaX = c(seq(from = 0, to = 0.5, by=0.1),seq(from = 0.75, to = 2, by=0.25))). The R code to run this example is given below.
```
library(MRdirection)
results<-MRdirection(nSim = 100, nX = 1000, nY=1000, MAF_GX = 0.5, MAF_GY = 0.5, gamma0 = 0, gammaGX = 1, gammaGY = 0, varX = 1, measurementError = F, beta0 = 0, betaX = c(seq(from = 0, to = 0.5, by=0.1),seq(from = 0.75, to = 2, by=0.25)), betaGX=rep(0,10), betaGY=rep(0.2,10), unmeasuredConfounding=F, varY = 1, sig.level = 0.05, SEED = 1, runMethods="All")
```

```
round(results$matrix,2)
```

## References
