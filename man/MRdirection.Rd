\name{MRdirection}
\alias{MRdirection}
\title{
MRDirection
}
\description{
Runs multiple MR approaches on simulated data for a user specified scenario and for each method returns the number of simulations that returned case 1 (X --> Y), case 2 (Y --> X), or an error.
}
\usage{
MRdirection(nSim = 10, nX = 1000, nY = 1000, MAF_GX = rep(0.5, 10), MAF_GY = rep(0.5, 10), gamma0 = 0, gammaGX = rep(0.2, 10), gammaGY = rep(0, 10), measurementError = F, delta0 = 0, deltaX = 1, varME = 1, beta0 = 0, betaX = seq(0, 2, by = 0.5), betaGX = rep(0, 10), betaGY = rep(0.2, 10), varY = 0.2, varX = 1, unmeasuredConfounding = F, meanU = 0, varU = 1, gammaU = 1, betaU = 1, eta0 = 0, etaGX = rep(0.2, 10), etaGY = rep(0, 10), snpEffectConfounder = F, long1 = F, long2 = F, kappa0 = 0, kappaY = 0.2, iota0 = 0, iotaX = 0.2, sig.level = 0.05, contX = TRUE, contY = TRUE, nSimTotal = 500, runMethods = "All", SEED = 1, path1, libPath = NULL, plot.pdf = F)
}
\arguments{
  \item{nSim}{
is the number of simulations that will be run.
}
  \item{nX}{
is the sample size of X that is being simulated.
}
  \item{nY}{
is the sample size of Y that is being simulated.
}
  \item{MAF_GX}{
is the vector of the minor allele frequencies of the SNPs GX.
}
  \item{MAF_GY}{
is the vector of the minor allele frequencies of the SNPs GY.
}
  \item{gamma0}{
is the intercept for the association of the SNPs with the exposure X.
}
  \item{gammaGX}{
is the vector of the slopes for the association of the SNPs GX with the exposure X.
}
  \item{gammaGY}{
is the vector of the slopes for the association of the SNPs GY with the exposure X.
}
  \item{measurementError}{
if measurementError=T, then the exposure is generated with measurement error.
}
  \item{delta0}{
is the intercept for the association of the exposure with measurement error.
}
  \item{deltaX}{
is the slope for the assosciation of the exposure with measurement error.
}
  \item{varME}{
is the variance of the exposure X with measurement error which is generated from a normal distribution.
}
  \item{beta0}{
is the intercept for the association of the exposure X with the outcome Y1.
}
  \item{betaX}{
is the slope for the association of the exposure X with the outcome Y1.
}
  \item{betaGX}{
is the vector of slopes for the association of the SNPs GX with the outcome Y1.
}
  \item{betaGY}{
is the vector of slopes for the association of the SNPs GY with the outcome Y1.
}
  \item{varY}{
is the variance of the outcome Y1 which is generated from a normal distribution.
}
  \item{varX}{
is the variance of the exposure X which is generated from a normal distribution.
}
  \item{unmeasuredConfounding}{
if unmeasuredConfounding=T, then an unmeasured confounder U of the exposure X outcome Y1 relationship is generated from a normal distribution with mean meanU and variance varU.
}
  \item{meanU}{
is the mean of the unmeasured confounder U.
}
  \item{varU}{
is the variance of the unmeasured confounder U.
}
  \item{gammaU}{
is the effect of the unmeasured confounder U on the exposure X.
}
  \item{betaU}{
is the effect of the unmeasured confounder U on the outcome Y1.
}
  \item{eta0}{
is the intercept for the association of the SNPs with the confounder.
}
  \item{etaGX}{
is the slope for the association of the SNPs GX with the confounder U.
}
  \item{etaGY}{
is the slope for the association of the SNPs GY with the confounder U.
}
  \item{snpEffectConfounder}{
if snpEffectConfounder=TRUE, SNPs have causal effect on confounder U.
}
  \item{long1}{
if long1=T, then X is generated from a longitudinal feedback loop such that X1 causes Y1 and Y1 causes X2. X2 and Y1 are used for phenotypes 1 (X) and 2 (Y1), respectively.
}
  \item{long2}{
if long2=T, then X and Y1 are generated from a longitudnial feedback loop such that X1 causes Y1, Y1 causes X2, and X2 causes Y2. X2 and Y2 are used for phenotypes 1 (X) and 2 (Y1), respectively.
}
  \item{kappa0}{
is the intercept for the association of Y1 with X2 (longitudinal feedback loop).
}
  \item{kappaY}{
is the slope for the association of Y1 with X2 (longitudinal feedback loop).
}
  \item{iota0}{
is the intercept for the association of X2 with Y2 (longitudinal feedback loop).
}
  \item{iotaX}{
is the slope for the association of X2 with Y2 (longitudinal feedback loop).
}
  \item{sig.level}{
is the significance level, default=0.05
}
  \item{contX}{
if contX=TRUE, X is generated from a normal distribution. If contX=FALSE, X is generated from a Bernoulli distribution.
}
  \item{contY}{
if contY=TRUE, Y1 is generated from a normal distirbution. If contY=FALSE, Y1 is generated from a Bernoulli distribution.
}
  \item{nSimTotal}{
is the total number of sims to be run for this scenario. It will be larger than nSim if the sims for the scenario are being split up and run (e.g., on multiple computers). Otherwise, it should be equal to nSim.
}
  \item{runMethods}{
is the methods to be run. If all methods should be run, runMethods="All" (this is the default). Other available options include: 
"MRS.ivw": will run MR Steiger (IVW) method from the TwoSampleMR package; 
"MRS.wMedian": will run MR Steiger (weighted median) method from the TwoSampleMR package; 
"MRS.Egger": will run MR Steiger (Egger) method from the TwoSampleMR package; 
"CDRatio": will run the CD-Ratio method from the MRCD package; 
"CDEgger": will run the CD-Egger method from the MRCD package; 
"CD-GLS": will run the CD-GLS method from the MRCD package; 
"BDCDcML.S.DP": runs the bidirectional CDcML method with data perturbation and screening from the BiDirectCausal pacakge; 
"BDMRcML.S.DP": runs the bidirectional MRcML method with data perturbation and screening from the BiDirectCausal pacakge; 
"BDCD.Ratio.S": runs the CD-Ratio method with screening from the BiDirectCausal package; 
"BDCD.Egger.S": runs the CD-Egger method with screening from the BiDirectCausal package; 
"tsmr_IVW": will run the bidirectional MR IVW method from the TwoSampleMR package; 
"tsmr_IVW_mre": will run the bidirectional MR IVW multiplicative random effects method from the TwoSampleMR package; 
"tsmr_IVW_fe": will run the bidirectional MR IVW fixed effects method from the TwoSampleMR package; 
"tsmr_simple_median": will run the bidirectional MR simple median method from the TwoSampleMR package; 
"tsmr_weighted_median": will run the bidirectional MR weighted median method from the TwoSampleMR package; 
"tsmr_pwm": : will run the bidirectional MR penalized weighted median method from the TwoSampleMR package; 
"tsmr_Egger": will run the bidirectional MR Egger method from the TwoSampleMR package; 
"tsmr_Egger_Boot": will run the bidirectional MR Egger Boot method from the TwoSampleMR package; 
"tsmr_uwr": will run the bidirectional MR unweighted ression method from the TwoSampleMR package.}
  \item{SEED}{
is the seed for the random number generator
}
  \item{path1}{
is the path to directory where .txt file output and plots should be saved.
}
  \item{libPath}{
is the path to directory packages are saved in if saved outside standard library directory. Leave as NULL if packages are saved in standard library directory.
}
  \item{plot.pdf}{
if TRUE, plots will be saved to the path1 directory.
}
}

\details{
This function outputs matrices of the number of simulations where case 1 (X --> Y), case 2(Y --> X), or an error is detected by using the MR methods specified by the user. These 3 matrices are also saved as .txt files to the path1 location.
This function also outputs plots of case 1 and case 2, which are saved to the the path1 location specified by the user.
}
\value{
\item{mat1}{matrix of the number of simulations that returned case 1 for each approach}
\item{mat2}{matrix of the number of simulations that returned case 2 for each approach}
\item{matE}{matrix of the number of simulations that returned an error for each approach}
}
\references{
G. Hemani, K. Tilling and G. Davey Smith, "Orienting the causal relationship between imprecisely measured traits using GWAS summary data," PLoS genetics, vol. 13, no. 11, 2017. 

H. Xue and W. Pan, "Inferring causal direction between two traits in the presence of horizontal pleiotropy with GWAS summary data," PLoS Genet, vol. 16, no. 11, 2020. 

H. Xue and W. Pan, "Robust inference of bi-directional causal relationships in presence of correlated pleiotropy with GWAS summary data," PLoS Genet, vol. 18, no. 5, 2022. 

G. Hemani, J. Zheng, B. Elsworth, K. Wade, D. Baird, V. Haberland, C. Laurin, S. Burgess, J. Bowden, R. Langdon, V. Tan, J. Yarmolinsky, H. Shihab, N. Timpson, D. Evans, C. Relton, R. Martin, G. Davey Smith, T. Gaunt, P. Haycock and T. M.-B. Collaboration., "The MR-Base platform supports systematic causal inference across the human phenome," eLife, vol. 7, 2018.
}
\author{
Sharon Lutz, Kirsten Voorhies
}


\examples{
tmpwd = getwd()
MRdirection(nSim=10, runMethods=c("All"), path1=tmpwd)
}

\section{Warning}{library(MRCD), library(TwoSampleMR), library(MRcML), and library(BiDirectCausal) are needed to run this function.}

\keyword{MR}
\keyword{MRdirection}
