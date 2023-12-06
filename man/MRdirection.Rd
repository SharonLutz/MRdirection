\name{MRdirection}
\alias{MRdirection}
\title{MRdirection}
\description{
Runs multiple MR approaches on simulated data under different scenarios specified by the user.
}
\usage{
MRdirection(nSim = 100, nX = 1000, nY = 1000, MAF_GX = rep(0.5, 10), MAF_GY = rep(0.5, 10), gamma0 = 0, gammaGX = rep(0.2, 10), gammaGY = rep(0, 10), measurementError = F, delta0 = 0, deltaX = 1, varME = 1, beta0 = 0, betaX = seq(0, 2, by = 0.5), betaGX = rep(0, 10), betaGY = rep(0.2, 10), varY = 0.2, varX = 1, unmeasuredConfounding = F, meanU = 0, varU = 1, gammaU = 1, betaU = 1, eta0 = 0, etaGX = rep(0.2, 10), etaGY = rep(0, 10), snpEffectConfounder = F, independentSNPs = T, long1 = F, long2 = F, kappa0 = 0, kappaY = 0.2, iota0 = 0, iotaX = 0.2, sig.level = 0.05, contX = TRUE, contY = TRUE, SEED = 1, runMethods = c("All"))
}

\arguments{
  \item{nSim}{
is the number of simulations that will be run
}
  \item{nX}{
is the sample size of X that is being simulated
}
  \item{nY}{
is the sample size of Y that is being simulated
}
  \item{MAF_GX}{
vector of the minor allele frequencies of the SNPs GX
}
  \item{MAF_GY}{
vector of the minor allele frequencies of the SNPs GY
}
  \item{gamma0}{
the intercept for the linear association of the SNPs with the exposure X
}
  \item{gammaGX}{
the vector of the slopes for the linear association of the SNPs GX with the exposure X
}
  \item{gammaGY}{
the vector of the slopes for the linear association of the SNPs GY with the exposure X
}
  \item{measurementError}{
if measurementError=T, then the exposure is generated with measurement error
}
  \item{delta0}{
the intercept for the linear association of the exposure with measurement error
}
  \item{deltaX}{
the slope for the linear assosciation of the exposure with measurement error
}
  \item{varME}{
the variance of the exposure X with measurement error which is generated from a normal distribution
}
  \item{beta0}{
the intercept for the linear association of the exposure X with the outcome Y1
}
  \item{betaX}{
the slope for the linear association of the exposure X with the outcome Y1
}
  \item{betaGX}{
the vector of slopes for the linear association of the SNPs GX with the outcome Y1
}
  \item{betaGY}{
the vector of slopes for the linear association of the SNPs GY with the outcome Y1
}
  \item{varY}{
the variance of the outcome Y1 which is generated from a normal distribution
}
  \item{varX}{
the variance of the exposure X which is generated from a normal distribution
}
  \item{unmeasuredConfounding}{
if unmeasuredConfounding=T, then an unmeasured confounder U of the exposure X outcome Y1 relationship is generated from a normal distribution with mean meanU and variance varU
}
  \item{meanU}{
the mean of the unmeasured confounder U
}
  \item{varU}{
the variance of the unmeasured confounder U
}
  \item{gammaU}{
the effect of the unmeasured confounder U on the exposure X
}
  \item{betaU}{
the effect of the unmeasured confounder U on the outcome Y1
}
  \item{eta0}{
the intercept for the linear association of the SNPs with the confounder
}
  \item{etaGX}{
the slope for the linear association of the SNPs GX with the confounder U
}
  \item{etaGY}{
the slope for the linear association of the SNPs GY with the confounder U
}
  \item{snpEffectConfounder}{
if snpEffectConfounder=TRUE, SNPs have causal effect on confounder U
}
  \item{independentSNPs}{
specifies whether the SNPs are independent.
}
  \item{long1}{
if long1=T, then X is generated from a longitudinal feedback loop such that X1 causes Y1 and Y1 causes X2. X2 and Y1 are used for phenotypes 1 (X) and 2 (Y1), respectively.
}
  \item{long2}{
if long2=T, then X and Y1 are generated from a longitudnial feedback loop such that X1 causes Y1, Y1 causes X2, and X2 causes Y2. X2 and Y2 are used for phenotypes 1 (X) and 2 (Y1), respectively.
}
  \item{kappa0}{
the intercept for the linear association of Y1 with X2 (longitudinal feedback loop)
}
  \item{kappaY}{
the slope for the linear association of Y1 with X2 (longitudinal feedback loop)
}
  \item{iota0}{
the intercept for the linear association of X2 with Y2 (longitudinal feedback loop)
}
\item{iotaX}{
the slope for the linear association of X2 with Y2 (longitudinal feedback loop)
}
  \item{sig.level}{
the significance level, default=0.05
}
  \item{contX}{
if contX=TRUE, X is generated from a normal distribution. If contX=FALSE, X is generated from a Bernoulli distribution.
}
  \item{contY}{
if contY=TRUE, Y1 is generated from a normal distirbution. If contY=FALSE, Y1 is generated from a Bernoulli distribution.
}
  \item{SEED}{
the seed for the random number generator
}
  \item{runMethods}{
the methods to test. Default is "All".
}
}
\details{
This function outputs matrices of the proportion of simulations where case 1 (X->Y), case 2(Y->X), or case 3 (inconclusive) is detected by using the MR methods specified by the user. These 3 matrices are also saved as .txt files to the working directory.
}
\value{

\item{mat_total1 }{matrix of the proportion of simulations that returned case 1 for each approach}
\item{mat_total2 }{matrix of the proportion of simulations that returned case 2 for each approach}
\item{mat_total3 }{matrix of the proportion of simulations that returned case 3 for each approach}

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

\section{Warning }{library(MRCD), library(TwoSampleMR), library(MRcML), and library(BiDirectCausal) are needed to run this function.}


\examples{
MRdirection(nSim=50, runMethods=c("All"))
}

\keyword{MR}
\keyword{MRdirection}
