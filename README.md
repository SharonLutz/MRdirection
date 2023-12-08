# MRdirection
Examines MR approaches to detect the directionality between phenotype 1 X and phenotype 2 Y through simulation studies.

## Installation
```
install.packages("devtools")
```
Note: if you already have devtools installed, you don't need to install it again.

```
devtools::install_github("xue-hr/MRCD")
```
Note: 
```
devtools::install_github("MRCIEU/TwoSampleMR")
```

```
devtools::install_github("xue-hr/MRcML")
```

```
devtools::install_github("xue-hr/BiDirectCausal")
```
Note: MRCD must be installed first

```
devtools::install_github("SharonLutz/MRdirection")
```

## Input
First, the SNPs ($G_X$ and $G_Y$) are generated from a binomial distribution for the max of nX and nY subjects (input nX, nY) for a given vector of minor allele frequencies (input MAF_GX and MAF_GY).
The true phenotype 1 (Xtrue) is generated from a normal distribution (input contX=T) or Bernoulli distribution (input contX=F) with the variance (input varX) and the mean as follows:

E\[Xtrue\] = $\gamma_0$ + $G_X*\gamma_{G_X}$ + $G_Y*\gamma_{G_Y}$

where if $\gamma_{G_Y}$ is nonzero then the SNPs $G_Y$ are pleiotropic.

All of these values are inputted by the user (i.e., the intercept gamma0, and the vectors of genetic effect sizes gammaGX and gammaGY). If there is no measurement error (input measurementError=F), then $X=Xtrue$. If there is measurement error (input measurementError=T), then the measured phenotype 1 X is generated from the true phenotype 1 $X_{true}$ such that:

E\[X\] = $\delta_0$ + $\delta_{X}X_{true}$

where $\delta_0$ and $\delta_{X}$ are input by the user (input delta0, deltaX). Phenotype 2 Y is generated from a normal distribution (contY=T) or Bernoulli distribution (contY=F) with the variance (input varY) and the mean as follows:

E\[Y1\] = $\beta_0$ + $\beta_{X}X_{true}$ + $G_Y\beta_{G_Y}$ + $G_X\beta_{G_X}$

where if $\beta_{G_X}$ is nonzero then the SNPs $G_X$ are pleiotropic.

All of these values are inputted by the user (i.e. the intercept beta0, the effect of phenotype 1 Xtrue on phenotype 2 as betaX, the vector of the effect of the SNPs $G_X$ directly on phenotype 2 as betaGX, and the vector of the effect of the SNPs $G_Y$ directly on phenotype 2 as betaGY).

If there is unmeasured confounding (unmeasuredConfounding=T) between the exposure X and the outcome Y, then the unmeasured confounder U is generated from a normal distribution with user specified mean and variance (i.e., meanU, varU). Then, the exposure X and outcome Y are generated such that

E\[Xtrue\] = $\gamma_0$ + $G_X\gamma_{G_X}$ + $G_Y\gamma_{G_Y}$ + $\gamma_{U}U$

E\[Y1\] = $\beta_0$ + $\beta_{X}X_{true}$ + $G_Y\beta_{G_Y}$ + $G_X\beta_{G_X}$ + $\beta_{U}U$

where $\gamma_{U}$ and $\beta_{U}$ are input by the user (i.e., gammaU and betaU).

After the data are generated, then the MRdirection function runs the specified approaches to determine which of the following three cases each approach concludes: case 1: the exposure X causes the outcome Y; case 2: the outcome Y causes the exposure X; case 3: inconclusive.

The following approaches can be specified by the user (input: runMethods):
"All": all approaches will be run. This is the default.

"MRSteigerMethods": runs the MR Steiger approach.

"CDMethods": runs the CD approaches.

"BidirectionalMethods": runs the twoSampleMR approaches.

The following approaches can be individually specified in a list as input for runMethods: "MRS.ivw", "MRS.wMedian", "MRS.Egger", "CDRatio", "CDEgger", "CDgls", "tsmr_IVW", "tsmr_Egger", "tsmr_weighted_median", "tsmr_uwr", "tsmr_simple_median", "tsmr_pwm", "tsmr_meta_fixed_simple", "tsmr_IVW_mre", "tsmr_IVW_fe", "tsmr_Egger_Boot", "BDCDcML.S.DP", "BDMRcML.S.DP", "BDCD.Ratio.S", "BDCD.Egger.S"


## Output
This function outputs matrices of the proportion of simulations where Case 1 (X->Y), Case 2 (Y->X), and Case 3 (inconclusive) are returned for each of the specified approaches. The matrices are also saved to the working directory.


## Example
Consider an example with 1000 subjects for both X and Y (input nX=1000 and nY=1000) and 10 SNPs each for X and Y with a MAF_GX of 50 (input MAF_GX=rep(0.5,10)) and MAF_GY of 50 (input MAF_GY=rep(0.5,10)). Consider no pleiotropy, measurement error, or unmeasured counfouding (input measurementError = F, unmeasuredConfounding=F, betaGX=0, gammaGY=0). Then, let X be generated from a normal distribution (input contX=T) with a variance of 1 (input varX = 1) and mean such that E\[Xtrue \] = 0 + 0.2\*GX (input gamma0=0, gammaGX=0.2). 
Y1 is generated from a normal distribution (input contY=T) with a variance of 0.2 (input varY = 0.2) and mean such that E\[Y1 \] = 0 + $\beta_{X}X$ + 0.2\*GY (input beta0 = 0, betaGY=0.2) and beta_{X} varies from 0 to 2 (input betaX = c(seq(from = 0, to = 0.5, by=0.25),seq(from = 0.75, to = 2, by=0.25))). All available MR approaches will be examined (input runMethods="All"). The R code to run this example is given below. 
```
library(MRdirection)
results<-MRdirection(nSim = 10, nX = 1000, nY=1000, MAF_GX = rep(0.5,10), MAF_GY = rep(0.5,10),
gamma0 = 0, gammaGX = rep(0.2,10), gammaGY = rep(0,10), varX = 1, measurementError = F, beta0 = 0,
betaX = c(seq(from = 0, to = 0.5, by=0.25),seq(from = 0.75, to = 2, by=0.25)),
betaGX=rep(0,10), betaGY=rep(0.2,10), unmeasuredConfounding=F, varY = 1,
sig.level = 0.05, SEED = 1, runMethods="All")
```

The function outputs the following matrices where each row corresponds to $\beta_{X}$ (input betaX). The proportion of simulations that concluded case 1 (X->Y) for the specified approaches are given in the first matrix (results$mat_total1). The proportion of simulations that concluded case 2 (Y->X) for the specified approaches are given in the second matrix (results$mat_total2). The proportion of simulations that concluded case 3 (inconclusive) for the specified approaches are given in the third matrix (results$mat_total3).

```
round(results$mat_total1, 2)

      betaX MRS.ivw MRS.wMedian MRS.Egger CDRatio CDEgger CDgls tsmr_IVW tsmr_Egger tsmr_weighted_median
 [1,]  0.00     0.2         0.1       0.0     0.1       0     0      0.1        0.0                  0.1
 [2,]  0.25     0.9         0.8       0.0     0.2       0     0      0.9        0.0                  0.8
 [3,]  0.50     1.0         0.9       0.0     0.7       0     0      0.8        0.0                  0.8
 [4,]  0.75     1.0         1.0       0.1     0.7       0     0      0.9        0.0                  0.8
 [5,]  1.00     1.0         1.0       0.4     0.5       0     0      0.9        0.3                  1.0
 [6,]  1.25     0.8         0.8       0.2     0.1       0     0      1.0        0.2                  1.0
 [7,]  1.50     0.2         0.2       0.2     0.0       0     0      0.8        0.4                  0.9
 [8,]  1.75     0.0         0.0       0.0     0.0       0     0      0.9        0.4                  1.0
 [9,]  2.00     0.0         0.0       0.0     0.0       0     0      0.6        0.3                  0.7
      tsmr_uwr tsmr_simple_median tsmr_pwm tsmr_meta_fixed_simple tsmr_IVW_mre tsmr_IVW_fe tsmr_Egger_Boot
 [1,]        0                0.1      0.1                    0.1          0.1         0.1             0.0
 [2,]        0                0.7      0.8                    1.0          0.9         1.0             0.0
 [3,]        0                0.8      0.8                    0.8          0.7         0.8             0.0
 [4,]        0                1.0      0.9                    0.8          0.8         0.8             0.0
 [5,]        0                1.0      1.0                    0.9          0.8         0.9             0.2
 [6,]        0                1.0      1.0                    1.0          0.9         1.0             0.2
 [7,]        0                1.0      0.9                    0.8          0.8         0.8             0.7
 [8,]        0                0.9      1.0                    0.9          0.7         0.9             0.6
 [9,]        0                1.0      0.8                    0.6          0.4         0.6             0.5
      BDCDcML.S.DP BDMRcML.S.DP BDCD.Ratio.S BDCD.Egger.S
 [1,]          0.2          0.0          0.3          0.1
 [2,]          0.1          0.0          0.1          0.0
 [3,]          0.6          0.3          0.4          0.0
 [4,]          0.7          0.5          0.2          0.0
 [5,]          0.7          0.7          0.7          0.0
 [6,]          0.7          0.6          0.5          0.0
 [7,]          0.9          0.9          0.9          0.0
 [8,]          1.0          1.0          0.9          0.0
 [9,]          1.0          1.0          1.0          0.0
```

```
round(results$mat_total2, 2)

      betaX MRS.ivw MRS.wMedian MRS.Egger CDRatio CDEgger CDgls tsmr_IVW tsmr_Egger tsmr_weighted_median
 [1,]  0.00       0           0         0       0       0     0      0.1        0.0                  0.1
 [2,]  0.25       0           0         0       0       0     0      0.0        0.0                  0.0
 [3,]  0.50       0           0         0       0       0     0      0.0        0.0                  0.0
 [4,]  0.75       0           0         0       0       0     0      0.0        0.2                  0.0
 [5,]  1.00       0           0         0       0       0     0      0.0        0.2                  0.0
 [6,]  1.25       0           0         0       0       0     0      0.0        0.1                  0.0
 [7,]  1.50       0           0         0       0       0     0      0.0        0.0                  0.0
 [8,]  1.75       0           0         0       0       0     0      0.0        0.2                  0.0
 [9,]  2.00       0           0         0       0       0     0      0.0        0.1                  0.0
      tsmr_uwr tsmr_simple_median tsmr_pwm tsmr_meta_fixed_simple tsmr_IVW_mre tsmr_IVW_fe tsmr_Egger_Boot
 [1,]        0                0.1      0.1                    0.1          0.1         0.1             0.0
 [2,]        0                0.0      0.0                    0.0          0.0         0.0             0.0
 [3,]        0                0.0      0.0                    0.0          0.0         0.0             0.0
 [4,]        0                0.0      0.0                    0.0          0.0         0.0             0.2
 [5,]        0                0.0      0.0                    0.0          0.0         0.0             0.0
 [6,]        0                0.0      0.0                    0.0          0.0         0.0             0.0
 [7,]        0                0.0      0.0                    0.0          0.0         0.0             0.0
 [8,]        0                0.0      0.0                    0.0          0.0         0.0             0.2
 [9,]        0                0.0      0.0                    0.0          0.0         0.0             0.1
      BDCDcML.S.DP BDMRcML.S.DP BDCD.Ratio.S BDCD.Egger.S
 [1,]            0            0            0            0
 [2,]            0            0            0            0
 [3,]            0            0            0            0
 [4,]            0            0            0            0
 [5,]            0            0            0            0
 [6,]            0            0            0            0
 [7,]            0            0            0            0
 [8,]            0            0            0            0
 [9,]            0            0            0            0
```

```
round(results$mat_total3, 2)

      betaX MRS.ivw MRS.wMedian MRS.Egger CDRatio CDEgger CDgls tsmr_IVW tsmr_Egger tsmr_weighted_median
 [1,]  0.00     0.8         0.9       1.0     0.9       1     1      0.8        1.0                  0.8
 [2,]  0.25     0.1         0.2       1.0     0.8       1     1      0.1        1.0                  0.2
 [3,]  0.50     0.0         0.1       1.0     0.3       1     1      0.2        1.0                  0.2
 [4,]  0.75     0.0         0.0       0.9     0.3       1     1      0.1        0.8                  0.2
 [5,]  1.00     0.0         0.0       0.6     0.5       1     1      0.1        0.5                  0.0
 [6,]  1.25     0.2         0.2       0.8     0.9       1     1      0.0        0.7                  0.0
 [7,]  1.50     0.8         0.8       0.8     1.0       1     1      0.2        0.6                  0.1
 [8,]  1.75     1.0         1.0       1.0     1.0       1     1      0.1        0.4                  0.0
 [9,]  2.00     1.0         1.0       1.0     1.0       1     1      0.4        0.6                  0.3
      tsmr_uwr tsmr_simple_median tsmr_pwm tsmr_meta_fixed_simple tsmr_IVW_mre tsmr_IVW_fe tsmr_Egger_Boot
 [1,]        1                0.8      0.8                    0.8          0.8         0.8             1.0
 [2,]        1                0.3      0.2                    0.0          0.1         0.0             1.0
 [3,]        1                0.2      0.2                    0.2          0.3         0.2             1.0
 [4,]        1                0.0      0.1                    0.2          0.2         0.2             0.8
 [5,]        1                0.0      0.0                    0.1          0.2         0.1             0.8
 [6,]        1                0.0      0.0                    0.0          0.1         0.0             0.8
 [7,]        1                0.0      0.1                    0.2          0.2         0.2             0.3
 [8,]        1                0.1      0.0                    0.1          0.3         0.1             0.2
 [9,]        1                0.0      0.2                    0.4          0.6         0.4             0.4
      BDCDcML.S.DP BDMRcML.S.DP BDCD.Ratio.S BDCD.Egger.S
 [1,]          0.8          1.0          0.7          0.9
 [2,]          0.9          1.0          0.9          1.0
 [3,]          0.4          0.7          0.6          1.0
 [4,]          0.3          0.5          0.8          1.0
 [5,]          0.3          0.3          0.3          1.0
 [6,]          0.3          0.4          0.5          1.0
 [7,]          0.1          0.1          0.1          1.0
 [8,]          0.0          0.0          0.1          1.0
 [9,]          0.0          0.0          0.0          1.0
```



## References
