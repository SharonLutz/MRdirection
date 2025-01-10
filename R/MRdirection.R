MRdirection <-
function(nSim=10, nX=1000, nY=1000,
           MAF_GX=rep(0.5,10), MAF_GY=rep(0.5,10),
           gamma0=0, gammaGX=rep(0.2,10), gammaGY=rep(0,10),
           measurementError=F, delta0=0, deltaX=1, varME=1,
           beta0=0, betaX=seq(0,2,by=0.5), betaGX=rep(0,10), betaGY=rep(0.2,10),
           varY=0.2, varX=1, unmeasuredConfounding=F, meanU=0, varU=1, gammaU=1, betaU=1, 
           eta0=0, etaGX=rep(0.2,10), etaGY=rep(0,10), snpEffectConfounder=F, 
           long1=F, long2=F, kappa0=0, kappaY=0.2, iota0=0,
           iotaX=0.2, sig.level=0.05, contX= TRUE, contY=TRUE, nSimTotal=500, runMethods="All",
           SEED=1, path1, libPath=NULL, plot.pdf=F){

    ################################################################################  
    # load libraries
    ################################################################################
    
    if(is.null(libPath)){
      library(MRCD)
      library(MRcML)
      library(BiDirectCausal)
      library(TwoSampleMR)
    }
    if(!is.null(libPath)){
      library(MRCD, lib.loc=libPath)
      library(MRcML, lib.loc=libPath)
      library(BiDirectCausal, lib.loc=libPath)
      library(TwoSampleMR, lib.loc=libPath)
    }
    
    ################################################################################
    # Simple Error checks
    ################################################################################
    if(nX<0 | nX==0 | floor(nX)!=ceiling(nX) ){stop("n must be an integer greater than 0")}
    if(nY<0 | nY==0 | floor(nY)!=ceiling(nY) ){stop("n must be an integer greater than 0")}
    if(nSim<0 | nSim==0 | floor(nSim)!=ceiling(nSim) ){stop("nSim must be an integer greater than 0")}
    if(length(sig.level)!=1| sig.level<0 | sig.level>1){stop("sig.level must be a single value greater than 0 and less than 1")}
    if(length(unique(betaX))<2){stop("betaX must be a vector with at least two values")}
    
    if(!varX>0){stop("varX must be greater than 0")}
    if(length(varX)!=1){stop("length(varX) must equal 1")}
    if(!varY>0){stop("varY must be greater than 0")}
    if(length(varY)!=1){stop("length(varY) must equal 1")}
    if(!varME>0){stop("varME must be greater than 0")}
    if(length(varME)!=1){stop("length(varME) must equal 1")}
    
    if(length(gamma0)!=1){stop("length(gamma0), the intercept, must equal 1")}
    if(length(delta0)!=1){stop("length(delta0), the intercept, must equal 1")}
    if(length(beta0)!=1){stop("length(beta0), the intercept, must equal 1")}
    if(length(eta0)!=1){stop("length(eta0), the intercept, must equal 1")}
    if(length(iota0)!=1){stop("length(iota0), the intercept, must equal 1")}
    if(length(kappa0)!=1){stop("length(kappa0), the intercept, must equal 1")}
    
    if(length(gammaGX)!=length(MAF_GX)){stop("length(MAF_GX) must equal length(gammaGX) for the same number of SNPs")}
    if(length(betaGX)!=length(MAF_GX)){stop("length(MAF_GX) must equal length(betaGX) for the same number of SNPs")}
    if(length(etaGX)!=length(MAF_GX)){stop("length(MAF_GX) must equal length(etaGX) for the same number of SNPs")}
    if(length(deltaX)!=1){stop("length(deltaX) must equal 1")}
    
    if(length(gammaGY)!=length(MAF_GY)){stop("length(MAF_GY) must equal length(gammaGY) for the same number of SNPs")}
    if(length(etaGY)!=length(MAF_GY)){stop("length(MAF_GY) must equal length(etaGY) for the same number of SNPs")}
    if(length(betaGY)!=length(MAF_GY)){stop("length(MAF_GY) must equal length(betaGY) for the same number of SNPs")}
    
    if(measurementError==T|measurementError=="T"|measurementError=="t"|measurementError=="True"|measurementError=="true"){measurementError<-T}  
    if(measurementError==F|measurementError=="F"|measurementError=="f"|measurementError=="False"|measurementError=="false"){measurementError<-F}
    if(measurementError!=T & measurementError!=F){stop("measurementError must equal True or False. Please note that R is case sensitive.")}
    
    if(unmeasuredConfounding==T|unmeasuredConfounding=="T"|unmeasuredConfounding=="t"|unmeasuredConfounding=="True"|unmeasuredConfounding=="true"){unmeasuredConfounding<-T}  
    if(unmeasuredConfounding==F|unmeasuredConfounding=="F"|unmeasuredConfounding=="f"|unmeasuredConfounding=="False"|unmeasuredConfounding=="false"){unmeasuredConfounding<-F}
    if(unmeasuredConfounding!=T & unmeasuredConfounding!=F){stop("unmeasuredConfounding must equal True or False. Please note that R is case sensitive.")}
    
    if(snpEffectConfounder==T|snpEffectConfounder=="T"|snpEffectConfounder=="t"|snpEffectConfounder=="True"|snpEffectConfounder=="true"){snpEffectConfounder<-T}  
    if(snpEffectConfounder==F|snpEffectConfounder=="F"|snpEffectConfounder=="f"|snpEffectConfounder=="False"|snpEffectConfounder=="false"){snpEffectConfounder<-F}
    if(snpEffectConfounder!=T & snpEffectConfounder!=F){stop("snpEffectConfounder must equal True or False. Please note that R is case sensitive.")}
    
    if(long1==T & long2==T){stop("long1 and long2 cannot both be true.")}  
    
    if(long1==T|long1=="T"|long1=="t"|long1=="True"|long1=="true"){long1<-T}  
    if(long1==F|long1=="F"|long1=="f"|long1=="False"|long1=="false"){long1<-F}
    if(long1!=T & long1!=F){stop("long1 must equal True or False. Please note that R is case sensitive.")}
    
    if(long2==T|long2=="T"|long2=="t"|long2=="True"|long2=="true"){long2<-T}  
    if(long2==F|long2=="F"|long2=="f"|long2=="False"|long2=="false"){long2<-F}
    if(long2!=T & long2!=F){stop("long2 must equal True or False. Please note that R is case sensitive.")}
    
    if(length(MAF_GX)==1){stop("The minimum number of SNPs required for X  is 2.")}
    if(length(MAF_GY)==1){stop("The minimum number of SNPsrequired for Y is 2.")}
    
    ################################################################################
    # Matrix to save Results
    ################################################################################
    
    AllMethods <- c("MRS.ivw","MRS.wMedian","MRS.Egger","CDRatio", "CDEgger","CDgls",
                    "BDCDcML.S.DP","BDMRcML.S.DP", "BDCD.Ratio.S", "BDCD.Egger.S",
                    "tsmr_IVW","tsmr_IVW_mre","tsmr_IVW_fe",
                    "tsmr_simple_median","tsmr_weighted_median","tsmr_pwm",
                    "tsmr_Egger","tsmr_Egger_Boot","tsmr_uwr")
                            
    MRSteigerMethods <- c("MRS.ivw","MRS.wMedian","MRS.Egger")
    
    CDMethods <- c("CDRatio","CDEgger","CDgls","BDCDcML.S.DP","BDMRcML.S.DP","BDCD.Ratio.S", "BDCD.Egger.S")
    
    BidirectionalMethods <- c("tsmr_IVW","tsmr_IVW_mre","tsmr_IVW_fe",
                            "tsmr_simple_median","tsmr_weighted_median","tsmr_pwm",
                            "tsmr_Egger","tsmr_Egger_Boot","tsmr_uwr")
    
    
    colnames.1 <- c("bX")
    if("All" %in% runMethods){colnames.1 <- c(colnames.1, AllMethods)}
    if("MRSteigerMethods" %in% runMethods){colnames.1 <- c(colnames.1, MRSteigerMethods)}
    if("CDMethods" %in% runMethods){colnames.1 <- c(colnames.1, CDMethods)}
    if("BidirectionalMethods" %in% runMethods){colnames.1 <- c(colnames.1, BidirectionalMethods)}

    if("MRS.ivw" %in% runMethods){colnames.1 <- c(colnames.1, "MRS.ivw")}
    if("MRS.wMedian" %in% runMethods){colnames.1 <- c(colnames.1, "MRS.wMedian")}
    if("MRS.Egger" %in% runMethods){colnames.1 <- c(colnames.1, "MRS.Egger")}
    
    if("CDRatio" %in% runMethods){colnames.1 <- c(colnames.1, "CDRatio")}
    if("CDEgger" %in% runMethods){colnames.1 <- c(colnames.1, "CDEgger")}
    if("CDgls" %in% runMethods){colnames.1 <- c(colnames.1, "CDgls")}
    
    if("tsmr_IVW" %in% runMethods){colnames.1 <- c(colnames.1, "tsmr_IVW")}
    if("tsmr_Egger" %in% runMethods){colnames.1 <- c(colnames.1, "tsmr_Egger")}
    
    if("tsmr_weighted_median" %in% runMethods){colnames.1 <- c(colnames.1, "tsmr_weighted_median")}
    if("tsmr_uwr" %in% runMethods){colnames.1 <- c(colnames.1, "tsmr_uwr")}
    if("tsmr_simple_median" %in% runMethods){colnames.1 <- c(colnames.1, "tsmr_simple_median")}
    if("tsmr_pwm" %in% runMethods){colnames.1 <- c(colnames.1, "tsmr_pwm")}
    
    if("tsmr_IVW_mre" %in% runMethods){colnames.1 <- c(colnames.1, "tsmr_IVW_mre")}
    if("tsmr_IVW_fe" %in% runMethods){colnames.1 <- c(colnames.1, "tsmr_IVW_fe")}
    if("tsmr_Egger_Boot" %in% runMethods){colnames.1 <- c(colnames.1, "tsmr_Egger_Boot")}
    
    if("BDCDcML.S.DP" %in% runMethods){colnames.1 <- c(colnames.1, "BDCDcML.S.DP")}
    if("BDMRcML.S.DP" %in% runMethods){colnames.1 <- c(colnames.1, "BDMRcML.S.DP")}
    if("BDCD.Ratio.S" %in% runMethods){colnames.1 <- c(colnames.1, "BDCD.Ratio.S")}
    if("BDCD.Egger.S" %in% runMethods){colnames.1 <- c(colnames.1, "BDCD.Egger.S")}
        
    colnames.methods <- unique(colnames.1)
    
    mat1 <- matrix(0,ncol=length(colnames.methods),nrow=length(betaX))
    mat2 <- matrix(0,ncol=length(colnames.methods),nrow=length(betaX))
    mat3 <- matrix(0,ncol=length(colnames.methods),nrow=length(betaX))
    matError <- matrix(0,ncol=length(colnames.methods),nrow=length(betaX))
    colnames.methods <- colnames.1
    
    colnames(mat1) <- colnames.methods
    colnames(mat2) <- colnames.methods
    colnames(mat3) <- colnames.methods
    colnames(matError) <- colnames.methods
    
    mat1[,"bX"] <- betaX
    mat2[,"bX"] <- betaX
    mat3[,"bX"] <- betaX
    matError[,"bX"] <- betaX    
    
    ################################################################################
    # cycle through the simulations
    ################################################################################  
    #cycle through the simulations
    for(ii in 1:nSim){
      printCut<-1
      if(floor(ii/printCut)==ceiling(ii/printCut)){print(paste(ii,"of",nSim,"simulations"))}
      ################################################################################
      # cycle through values of betaX
      ################################################################################
      #cycle through values of betaX
      for(bX in 1:length(betaX)){
        print(bX)
        seedPlus = (bX * nSimTotal)+ii
        set.seed(SEED+seedPlus-(nSimTotal+1))
        
        ################################################################################
        # generate data
        ################################################################################                            
        # Generate GX and GY (SNPs)
        nSNP_GX<-length(MAF_GX) # number of SNPs
        nSNP_GY<-length(MAF_GY)
        GX<-matrix(0,nrow=max(nX,nY),ncol=nSNP_GX) 
        GY<-matrix(0,nrow=max(nX,nY),ncol=nSNP_GY) 
        
        for(mm in 1:nSNP_GX){
          GX[,mm] <- rbinom(max(nX,nY),2,MAF_GX[mm])
        }
        
        for(mm in 1:nSNP_GY){
          GY[,mm] <- rbinom(max(nX,nY),2,MAF_GY[mm])
        }
        
        # Generate unmeasured confounder U
        if(snpEffectConfounder==F){# if FALSE, confounder U not affected by SNP GX or GY
          U<-rnorm(max(nX,nY),mean=meanU,sd=sqrt(varU))}
        if(snpEffectConfounder==T){ # if TRUE, SNP G has causal effect on confounder U
          U<-rnorm(max(nX,nY),mean=(eta0 + GX%*%etaGX + GY%*%etaGY),sd=sqrt(varU))}
        
        # Generate X (exposure/ intermediate phenotype)
        gammaGX<-matrix(gammaGX,nrow=length(gammaGX),ncol=1)
        gammaGY<-matrix(gammaGY,nrow=length(gammaGY),ncol=1)
        
        # Unmeasured confouding True or False
        if(unmeasuredConfounding==F & contX==T){
          Xtrue <- gamma0 + GX%*%gammaGX + GY%*%gammaGY + rnorm(max(nX,nY),mean=0,sd=sqrt(varX)) }
        if(unmeasuredConfounding==T & contX==T){
          Xtrue <- gamma0 + GX%*%gammaGX + GY%*%gammaGY+ gammaU*U+ rnorm(max(nX,nY),mean=0,sqrt(varX)) }
        
        if(unmeasuredConfounding==F & contX==F){
          muTrue <- gamma0 + GX%*%gammaGX + GY%*%gammaGY
          muTrue <- muTrue - mean(muTrue)
          Xtrue <- rbinom(max(nX,nY),1, exp(muTrue)/(1+exp(muTrue)))}
        if(unmeasuredConfounding==T & contX==F){
          muTrue <- gamma0 + GX%*%gammaGX + GY%*%gammaGY+ gammaU*U  
          muTrue <- muTrue - mean(muTrue)
          Xtrue <- rbinom(max(nX,nY),1, exp(muTrue)/(1+exp(muTrue)))}
        
        
        # Measurement error True or False
        if(measurementError==T & contX==T){
          X<-rnorm(max(nX,nY),(delta0+deltaX*Xtrue),sqrt(varME))}
        if(measurementError==T & contX==F){
          mu1<-delta0+deltaX*Xtrue
          mu1 <- mu1 - mean(mu1)
          X <- rbinom(max(nX,nY),1,exp(mu1)/(1+exp(mu1)))}
        if(measurementError==F){
          X<-Xtrue}
        
        # Generate Y1 (outcome)
        if(unmeasuredConfounding==F & contY==T){
          Y1 <- beta0 + betaX[bX]*Xtrue +GX%*%betaGX + GY%*%betaGY+rnorm(max(nX,nY),mean=0,sd=sqrt(varY))}
        if(unmeasuredConfounding==T & contY==T){
          Y1 <- beta0 + betaX[bX]*Xtrue +GX%*%betaGX + GY%*%betaGY+betaU*U+rnorm(max(nX,nY),mean=0,sd=sqrt(varY))}
        
        if(unmeasuredConfounding==F & contY==F){
          mu1 <- beta0 + betaX[bX]*Xtrue +GX%*%betaGX + GY%*%betaGY
          mu1 <- mu1 - mean(mu1)
          Y1 <- rbinom(max(nX,nY), 1, exp(mu1)/(1+exp(mu1)))}
        if(unmeasuredConfounding==T & contY==F){
          mu1 <- beta0 + betaX[bX]*Xtrue +GX%*%betaGX + GY%*%betaGY+betaU*U
          mu1 <- mu1 - mean(mu1)
          Y1 <- rbinom(max(nX,nY), 1, exp(mu1)/(1+exp(mu1)))}
        
        
        # Generate X2 if long1==T or long2==T, and Y2 if long2==T
        if(long1==T & contX==T){
          X2 <- rnorm(max(nX,nY),kappa0 + kappaY*Y1, 1)
          X <- X2
        }
        if(long1==T & contX==F){
          mu2 <-kappa0 + kappaY*Y1
          mu2 <- mu2 - mean(mu2)
          X2 <- rbinom(max(nX,nY),1,exp(mu2)/(1+exp(mu2)))
          X <- X2
        }
        
        if(long2==T & contX==T & contY == T){
          X2 <- rnorm(max(nX,nY), kappa0 + kappaY*Y1, 1)
          Y2 <- rnorm(max(nX,nY), iota0 + iotaX*X2, 1)
          X<- X2
          Y1 <- Y2
        }
        if(long2==T & contX==F & contY == T){
          mu2 <-kappa0 + kappaY*Y1
          mu2 <- mu2 - mean(mu2)
          X2 <- rbinom(max(nX,nY),1,exp(mu2)/(1+exp(mu2)))
          Y2 <- rnorm(max(nX,nY), iota0 + iotaX*X2, 1)
          X<- X2
          Y1 <- Y2
        }
        if(long2==T & contX==T & contY == F){
          X2 <- rnorm(max(nX,nY), kappa0 + kappaY*Y1, 1)
          mu3 <- iota0 + iotaX*X2
          mu3 <- mu3 - mean(mu3)
          Y2 <- rbinom(max(nX,nY),1,exp(mu3)/(1+exp(mu3)))
          X<- X2
          Y1 <- Y2
        }
        if(long2==T & contX==F & contY == F){
          mu2 <- kappa0 + kappaY*Y1
          mu2 <- mu2 - mean(mu2)
          X2 <- rbinom(max(nX,nY),1,exp(mu2)/(1+exp(mu2)))
          mu3 <- iota0 + iotaX*X2
          mu3 <- mu3 - mean(mu3)
          Y2 <- rbinom(max(nX,nY),1,exp(mu3)/(1+exp(mu3)))
          X<- X2
          Y1 <- Y2
        }
        
        GX1 <- GX
        GX2 <- GX
        GY1 <- GY
        GY2 <- GY
        
        # Differing sample sizes
        if(nX!=nY){
          sampX <- sample(c(1:max(nX,nY)),size=nX)
          sampY <- sample(c(1:max(nX,nY)),size=nY)
          
          GX1 <- GX[sampX,]
          GY1 <- GY[sampX,]
          X <- X[sampX]
          
          GX2 <- GX[sampY,]
          GY2 <- GY[sampY,]
          Y1 <- Y1[sampY]
        }
        
        #################################################################################################################################################
        # More than 1 SNP based on code from mr_steiger function in from TwoSampleMR package, Peters et al 2014 codeANM package, MRCD package
        # https://github.com/MRCIEU/TwoSampleMR
        # http://web.math.ku.dk/~peters/code.html
        # https://github.com/xue-hr/MRCD
        #################################################################################################################################################

        # BIDIRECTIONAL MR
        # GX --> X --> Y
        betaGX_X <- NULL
        seGX_X <- NULL
        betaGX_Y <- NULL
        seGX_Y <- NULL
        
        if(contX==T){
          # save beta/se for each snp GX -> Y
          # Get matrices
          nM<-nX
          nSNP <- ncol(GX)
          Xmat<-GX1 # dim nrow=nM ncol=nSNP
          Ymat<-matrix(X,nrow=nM,ncol=1)
          ####
          Xbar<- matrix(colSums(Xmat)/nM,ncol=nSNP,nrow=nM,byrow=T)
          # numerator sum of (x-mean(x))*(y-mean(y))
          Beta1N<- t(Xmat- Xbar)%*% (Ymat-mean(Ymat))
          #denominator sum of (x-mean(x))^2
          Beta1D<- colSums((Xmat- Xbar )^2)
          #beta1 vector
          betaHat<- Beta1N/ Beta1D
          #beta0 vector
          betaO<-matrix(mean(Ymat),nrow=nSNP,ncol=1)- betaHat* matrix(colSums(Xmat)/nM,ncol=1,nrow=nSNP)
          # SE for beta
          Yhat<-(matrix(betaO,nrow=nM,ncol=nSNP,byrow=T)+ matrix(betaHat,nrow=nM,ncol=nSNP,byrow=T) * Xmat  )
          sigma2<- colSums(( matrix(Ymat,ncol=nSNP,nrow=nM)- Yhat )^2)/(nM-2)
          seBeta<-sqrt( sigma2/Beta1D)
          #####
          betaGX_X <- betaHat
          seGX_X <- seBeta
          
        }
        
        # save beta/se for each snp GX -> Y
        # Get matrices
        if(contY==T){
          nM<-nY
          nSNP <- ncol(GX2)
          Xmat<-GX2 # dim nrow=nM ncol=nSNP
          Ymat<-matrix(Y1,nrow=nM,ncol=1)
          ####
          Xbar<- matrix(colSums(Xmat)/nM,ncol=nSNP,nrow=nM,byrow=T)
          # numerator sum of (x-mean(x))*(y-mean(y))
          Beta1N<- t(Xmat- Xbar)%*% (Ymat-mean(Ymat))
          #denominator sum of (x-mean(x))^2
          Beta1D<- colSums((Xmat- Xbar )^2)
          #beta1 vector
          betaHat<- Beta1N/ Beta1D
          #beta0 vector
          betaO<-matrix(mean(Ymat),nrow=nSNP,ncol=1)- betaHat* matrix(colSums(Xmat)/nM,ncol=1,nrow=nSNP)
          # SE for beta
          Yhat<-(matrix(betaO,nrow=nM,ncol=nSNP,byrow=T)+ matrix(betaHat,nrow=nM,ncol=nSNP,byrow=T) * Xmat  )
          sigma2<- colSums(( matrix(Ymat,ncol=nSNP,nrow=nM)- Yhat )^2)/(nM-2)
          seBeta<-sqrt( sigma2/Beta1D)
          #####
          betaGX_Y <- betaHat
          seGX_Y <- seBeta
        }
        
        if(contX==F){
          for(i in 1:ncol(GX1)){
            betaGX_X <- c(betaGX_X, summary(glm(X~GX1[,i], family="binomial"))$coef[2,1])
            seGX_X <- c(seGX_X, summary(glm(X~GX1[,i], family="binomial"))$coef[2,2])
          }
        }
        
        if(contY==F){
          for(i in 1:ncol(GX1)){
            betaGX_Y <- c(betaGX_Y, summary(glm(Y1~GX2[,i], family="binomial"))$coef[2,1])
            seGX_Y <- c(seGX_Y, summary(glm(Y1~GX2[,i], family="binomial"))$coef[2,2])
          }
        }
        
        # GY --> Y --> X
        betaGY_Y <- NULL
        seGY_Y <- NULL          
        betaGY_X <- NULL
        seGY_X <- NULL
        
        if(contX==T){
          # save beta/se for each snp GX -> Y
          # Get matrices
          nM<-nX
          nSNP <- ncol(GY1)
          Xmat<-GY1 # dim nrow=nM ncol=nSNP
          Ymat<-matrix(X,nrow=nM,ncol=1)
          ####
          Xbar<- matrix(colSums(Xmat)/nM,ncol=nSNP,nrow=nM,byrow=T)
          # numerator sum of (x-mean(x))*(y-mean(y))
          Beta1N<- t(Xmat- Xbar)%*% (Ymat-mean(Ymat))
          #denominator sum of (x-mean(x))^2
          Beta1D<- colSums((Xmat- Xbar )^2)
          #beta1 vector
          betaHat<- Beta1N/ Beta1D
          #beta0 vector
          betaO<-matrix(mean(Ymat),nrow=nSNP,ncol=1)- betaHat* matrix(colSums(Xmat)/nM,ncol=1,nrow=nSNP)
          # SE for beta
          Yhat<-(matrix(betaO,nrow=nM,ncol=nSNP,byrow=T)+ matrix(betaHat,nrow=nM,ncol=nSNP,byrow=T) * Xmat  )
          sigma2<- colSums(( matrix(Ymat,ncol=nSNP,nrow=nM)- Yhat )^2)/(nM-2)
          seBeta<-sqrt( sigma2/Beta1D)
          #####
          betaGY_X <- betaHat
          seGY_X <- seBeta
        }
        
        if(contY==T){
          nM<-nY
          nSNP <- ncol(GY2)
          Xmat<-GY2 # dim nrow=nM ncol=nSNP
          Ymat<-matrix(Y1,nrow=nM,ncol=1)
          ####
          Xbar<- matrix(colSums(Xmat)/nM,ncol=nSNP,nrow=nM,byrow=T)
          # numerator sum of (x-mean(x))*(y-mean(y))
          Beta1N<- t(Xmat- Xbar)%*% (Ymat-mean(Ymat))
          #denominator sum of (x-mean(x))^2
          Beta1D<- colSums((Xmat- Xbar )^2)
          #beta1 vector
          betaHat<- Beta1N/ Beta1D
          #beta0 vector
          betaO<-matrix(mean(Ymat),nrow=nSNP,ncol=1)- betaHat* matrix(colSums(Xmat)/nM,ncol=1,nrow=nSNP)
          # SE for beta
          Yhat<-(matrix(betaO,nrow=nM,ncol=nSNP,byrow=T)+ matrix(betaHat,nrow=nM,ncol=nSNP,byrow=T) * Xmat  )
          sigma2<- colSums(( matrix(Ymat,ncol=nSNP,nrow=nM)- Yhat )^2)/(nM-2)
          seBeta<-sqrt( sigma2/Beta1D)
          #####
          betaGY_Y <- betaHat
          seGY_Y <- seBeta
        }
        
        if(contX==F){
          for(i in 1:ncol(GY1)){
            betaGY_X <- c(betaGY_X, summary(glm(X~GY1[,i], family="binomial"))$coef[2,1])
            seGY_X <- c(seGY_X, summary(glm(X~GY1[,i], family="binomial"))$coef[2,2])
          }
        }
        
        if(contY==F){
          for(i in 1:ncol(GY2)){
            betaGY_Y <- c(betaGY_Y, summary(glm(Y1~GY2[,i], family="binomial"))$coef[2,1])
            seGY_Y <- c(seGY_Y, summary(glm(Y1~GY2[,i], family="binomial"))$coef[2,2])
          }
        }
        
        betaGX_X_1<-betaGX_X[abs(betaGX_X)>1e-8 & abs(betaGX_Y)>1e-8]
        seGX_X_1<-seGX_X[abs(betaGX_X)>1e-8 & abs(betaGX_Y)>1e-8]
        
        # To remove betas that round to zero
        betaGX_Y_1<-betaGX_Y[abs(betaGX_Y)>1e-8 & abs(betaGX_X)>1e-8]
        seGX_Y_1<-seGX_Y[abs(betaGX_Y)>1e-8 & abs(betaGX_X)>1e-8]
        
        # To remove betas that round to zero
        betaGY_X_1<-betaGY_X[abs(betaGY_X)>1e-8 & abs(betaGY_Y)>1e-8]
        seGY_X_1<-seGY_X[abs(betaGY_X)>1e-8 & abs(betaGY_Y)>1e-8]
        
        # To remove betas that round to zero
        betaGY_Y_1<-betaGY_Y[abs(betaGY_Y)>1e-8 & abs(betaGY_X)>1e-8]
        seGY_Y_1<-seGY_Y[abs(betaGY_Y)>1e-8 & abs(betaGY_X)>1e-8]
        
        betaGX_X <- betaGX_X_1
        seGX_X <- seGX_X_1
        betaGX_Y <- betaGX_Y_1
        seGX_Y <- seGX_Y_1
        
        betaGY_X <- betaGY_X_1
        seGY_X <- seGY_X_1
        betaGY_Y <- betaGY_Y_1
        seGY_Y <- seGY_Y_1
        
        # define b_X, b_Y, se_X, se_Y for CD screening analyses
        b_X= c(betaGX_X, betaGY_X)
        b_Y=c(betaGX_Y, betaGY_Y)
        se_X= c(seGX_X, seGY_X)
        se_Y= c(seGX_Y, seGY_Y)
        
       ####################################################################################
       # start bidirectional MR approaches (from TwoSampleMR package) 
       ####################################################################################
        # Run bidirectional MR with TwoSampleMR::mr_ivw()
        if( ("tsmr_IVW" %in% colnames.methods) | ("MRS.ivw" %in% colnames.methods) ){
            tsmr.IVW.GX <- TwoSampleMR::mr_ivw(b_exp=betaGX_X, se_exp = seGX_X, b_out=betaGX_Y, se_out = seGX_Y, parameters = default_parameters())
            tsmr.IVW.GY <- TwoSampleMR::mr_ivw(b_exp=betaGY_Y, se_exp = seGY_Y, b_out=betaGY_X, se_out = seGY_X, parameters = default_parameters())
            # TwoSampleMR mr_ivw cases
            if( ("tsmr_IVW" %in% colnames.methods) ){
                
              if(tsmr.IVW.GX$pval <= sig.level & tsmr.IVW.GY$pval>sig.level){mat1[bX, "tsmr_IVW"] <- mat1[bX, "tsmr_IVW"] + 1
              }else if(tsmr.IVW.GX$pval > sig.level & tsmr.IVW.GY$pval<=sig.level){mat2[bX, "tsmr_IVW"] <- mat2[bX, "tsmr_IVW"] + 1
              }
            }
        }
        
        # Run bidirectional MR with TwoSampleMR::mr_egger_regression()
        if( ("tsmr_Egger" %in% colnames.methods) | ("MRS.Egger" %in% colnames.methods)){
            tsmr.Egger.GX<- TwoSampleMR::mr_egger_regression(b_exp=betaGX_X, se_exp = seGX_X, b_out=betaGX_Y, se_out = seGX_Y, parameters = default_parameters())
            tsmr.Egger.GY<- TwoSampleMR::mr_egger_regression(b_exp=betaGY_Y, se_exp = seGY_Y, b_out=betaGY_X, se_out = seGY_X, parameters = default_parameters())
            # TwoSampleMR mr_egger cases
            if("tsmr_Egger" %in% colnames.methods){
              if(tsmr.Egger.GX$pval <= sig.level & tsmr.Egger.GY$pval>sig.level){mat1[bX, "tsmr_Egger"] <- mat1[bX, "tsmr_Egger"] + 1
              }else if(tsmr.Egger.GX$pval > sig.level & tsmr.Egger.GY$pval<=sig.level){mat2[bX, "tsmr_Egger"] <- mat2[bX, "tsmr_Egger"] + 1
              }
            }
        }
        
        # Run bidirectional MR with TwoSampleMR::mr_egger_regression_bootstrap()
        if("tsmr_Egger_Boot" %in% colnames.methods){
            tsmr.Egger_boot.GX<- TwoSampleMR::mr_egger_regression_bootstrap(b_exp=betaGX_X, se_exp = seGX_X, b_out=betaGX_Y, se_out = seGX_Y, parameters = default_parameters())
            tsmr.Egger_boot.GY<- TwoSampleMR::mr_egger_regression_bootstrap(b_exp=betaGY_Y, se_exp = seGY_Y, b_out=betaGY_X, se_out = seGY_X, parameters = default_parameters())
            # TwoSampleMR mr_egger cases
            if(tsmr.Egger_boot.GX$pval <= sig.level & tsmr.Egger_boot.GY$pval>sig.level){mat1[bX, "tsmr_Egger_Boot"] <- mat1[bX, "tsmr_Egger_Boot"] + 1
            }else if(tsmr.Egger_boot.GX$pval > sig.level & tsmr.Egger_boot.GY$pval<=sig.level){mat2[bX, "tsmr_Egger_Boot"] <- mat2[bX, "tsmr_Egger_Boot"] + 1
            }
        }
        
        # Run bidirectional MR with TwoSampleMR::mr_ivw_fe() (Inverse variance weighted regression - fixed effects)
        if("tsmr_IVW_fe" %in% colnames.methods){
            tsmr.IVW_fe.GX <- TwoSampleMR::mr_ivw_fe(b_exp=betaGX_X, se_exp = seGX_X, b_out=betaGX_Y, se_out = seGX_Y, parameters = default_parameters())
            tsmr.IVW_fe.GY <- TwoSampleMR::mr_ivw_fe(b_exp=betaGY_Y, se_exp = seGY_Y, b_out=betaGY_X, se_out = seGY_X, parameters = default_parameters())
            # TwoSampleMR mr_ivw_fe cases
            if(tsmr.IVW_fe.GX$pval <= sig.level & tsmr.IVW_fe.GY$pval>sig.level){mat1[bX, "tsmr_IVW_fe"] <- mat1[bX, "tsmr_IVW_fe"] + 1
            }else if(tsmr.IVW_fe.GX$pval > sig.level & tsmr.IVW_fe.GY$pval<=sig.level){mat2[bX, "tsmr_IVW_fe"] <- mat2[bX, "tsmr_IVW_fe"] + 1
            }
        }
        
        # Run bidirectional MR with TwoSampleMR::mr_ivw_mre() (Inverse variance weighted regression - multiplicative random effects model)
        if("tsmr_IVW_mre" %in% colnames.methods){
            tsmr.IVW_mre.GX <- TwoSampleMR::mr_ivw_mre(b_exp=betaGX_X, se_exp = seGX_X, b_out=betaGX_Y, se_out = seGX_Y, parameters = default_parameters())
            tsmr.IVW_mre.GY <- TwoSampleMR::mr_ivw_mre(b_exp=betaGY_Y, se_exp = seGY_Y, b_out=betaGY_X, se_out = seGY_X, parameters = default_parameters())
            # TwoSampleMR mr_ivw_mre cases
            if(tsmr.IVW_mre.GX$pval <= sig.level & tsmr.IVW_mre.GY$pval>sig.level){mat1[bX, "tsmr_IVW_mre"] <- mat1[bX, "tsmr_IVW_mre"] + 1
            }else if(tsmr.IVW_mre.GX$pval > sig.level & tsmr.IVW_mre.GY$pval<=sig.level){mat2[bX, "tsmr_IVW_mre"] <- mat2[bX, "tsmr_IVW_mre"] + 1
            }
        }
        
        
        # Run bidirectional MR with TwoSampleMR::mr_penalized_weighted_median() (penalized weighted median MR)
        if("tsmr_pwm" %in% colnames.methods){
            tsmr.penalized_weighted_median.GX <- TwoSampleMR::mr_penalised_weighted_median(b_exp=betaGX_X_1, se_exp = seGX_X_1, b_out=betaGX_Y_1, se_out = seGX_Y_1, parameters = default_parameters())
            tsmr.penalized_weighted_median.GY <- TwoSampleMR::mr_penalised_weighted_median(b_exp=betaGY_Y_1, se_exp = seGY_Y_1, b_out=betaGY_X_1, se_out = seGY_X_1, parameters = default_parameters())
            # TwoSampleMR mr_penalized_weighted_median cases
            if(tsmr.penalized_weighted_median.GX$pval <= sig.level & tsmr.penalized_weighted_median.GY$pval>sig.level){mat1[bX, "tsmr_pwm"] <- mat1[bX, "tsmr_pwm"] + 1
            }else if(tsmr.penalized_weighted_median.GX$pval > sig.level & tsmr.penalized_weighted_median.GY$pval<=sig.level){mat2[bX, "tsmr_pwm"] <- mat2[bX, "tsmr_pwm"] + 1
            }
        }
        
        
        # Run bidirectional MR with TwoSampleMR::mr_simple_median
        if("tsmr_simple_median" %in% colnames.methods){
            tsmr.simple_median.GX <- TwoSampleMR::mr_simple_median(b_exp=betaGX_X, se_exp = seGX_X, b_out=betaGX_Y, se_out = seGX_Y, parameters = default_parameters())
            tsmr.simple_median.GY <- TwoSampleMR::mr_simple_median(b_exp=betaGY_Y, se_exp = seGY_Y, b_out=betaGY_X, se_out = seGY_X, parameters = default_parameters())
            # TwoSampleMR mr_simple_median cases
            if(tsmr.simple_median.GX$pval <= sig.level & tsmr.simple_median.GY$pval>sig.level){mat1[bX, "tsmr_simple_median"] <- mat1[bX, "tsmr_simple_median"] + 1
            }else if(tsmr.simple_median.GX$pval > sig.level & tsmr.simple_median.GY$pval<=sig.level){mat2[bX, "tsmr_simple_median"] <- mat2[bX, "tsmr_simple_median"] + 1
            }
        }
        
        
        # Run bidirectional MR with TwoSampleMR::mr_uwr (Unweighted regression)
        if("tsmr_uwr" %in% colnames.methods){
            tsmr.uwr.GX <- TwoSampleMR::mr_uwr(b_exp=betaGX_X, se_exp = seGX_X, b_out=betaGX_Y, se_out = seGX_Y, parameters = default_parameters())
            tsmr.uwr.GY <- TwoSampleMR::mr_uwr(b_exp=betaGY_Y, se_exp = seGY_Y, b_out=betaGY_X, se_out = seGY_X, parameters = default_parameters())
            # TwoSampleMR mr_uwr cases
            if(tsmr.uwr.GX$pval <= sig.level & tsmr.uwr.GY$pval>sig.level){mat1[bX, "tsmr_uwr"] <- mat1[bX, "tsmr_uwr"] + 1
            }else if(tsmr.uwr.GX$pval > sig.level & tsmr.uwr.GY$pval<=sig.level){mat2[bX, "tsmr_uwr"] <- mat2[bX, "tsmr_uwr"] + 1
            }
        }
        
        # Run bidirectional MR with TwoSampleMR::mr_weighted_median
        if( ("tsmr_weighted_median" %in% colnames.methods) | ("MRS.wMedian" %in% colnames.methods) ){
            tsmr.weighted_median.GX <- TwoSampleMR::mr_weighted_median(b_exp=betaGX_X, se_exp = seGX_X, b_out=betaGX_Y, se_out = seGX_Y, parameters = default_parameters())
            tsmr.weighted_median.GY <- TwoSampleMR::mr_weighted_median(b_exp=betaGY_Y, se_exp = seGY_Y, b_out=betaGY_X, se_out = seGY_X, parameters = default_parameters())
            # TwoSampleMR mr_weighted_median cases
            if("tsmr_weighted_median" %in% colnames.methods){
              if(tsmr.weighted_median.GX$pval <= sig.level & tsmr.weighted_median.GY$pval>sig.level){mat1[bX, "tsmr_weighted_median"] <- mat1[bX, "tsmr_weighted_median"] + 1
              }else if(tsmr.weighted_median.GX$pval > sig.level & tsmr.weighted_median.GY$pval<=sig.level){mat2[bX, "tsmr_weighted_median"] <- mat2[bX, "tsmr_weighted_median"] + 1
              }
            }
        }
        
        ####################################################################################
        # Start CD approaches (MRCD package)
        ####################################################################################
        # CD-Ratio
        
        if(("CDRatio" %in% colnames.methods)  | ("CDEgger" %in% colnames.methods) | ("CDgls" %in% colnames.methods)){

        
            tmp.mat <- matrix(NA,ncol=12,nrow=length(betaGX_X_1))
            tmp.df <- as.data.frame(tmp.mat)
            colnames(tmp.df) <- c("chr", "pos", "rsid", "A1", "A2", "beta_T1",
                                  "se_T1", "N_T1", "beta_T2", "se_T2", "N_T2", "loci")
            tmp.df$rsid<-paste0("SNP", c(1:length(betaGX_X_1)))
            tmp.df$beta_T1 <- betaGX_X_1
            tmp.df$se_T1 <- seGX_X_1
            tmp.df$N_T1 <- nX
            tmp.df$beta_T2 <- betaGX_Y_1
            tmp.df$se_T2 <- seGX_Y_1
            tmp.df$N_T2 <- nY
            
            pruned1 <- list(sig_part=tmp.df)
            mtry <- try(cd3 <- CD_3_methods_Independent(pruned1$sig_part))
            
            if(inherits(mtry, "try-error")){
              if("CDRatio" %in% colnames.methods){matError[bX, "CDRatio"] <- matError[bX, "CDRatio"] + 1}
              if("CDEgger" %in% colnames.methods){matError[bX, "CDEgger"] <- matError[bX, "CDEgger"] + 1}
              if("CDgls" %in% colnames.methods){matError[bX, "CDgls"] <- matError[bX, "CDgls"] + 1}
            }
            
            if(!inherits(mtry, "try-error")){
              ratio.YX <- cd3$CD_Ratio_result$T1toT2
              ratio.XY <- cd3$CD_Ratio_result$T2toT1
              
              egger.YX<- cd3$CD_Egger_result$T1toT2
              egger.XY <- cd3$CD_Egger_result$T2toT1
              
              gls.YX <- cd3$CD_GLS_result$T1toT2
              gls.XY <- cd3$CD_GLS_result$T2toT1
              
              # Confidence intervals for K - needed for use in decision rules
              # CD-ratio CIs
              lowerCIyx.cdRatio <- ratio.YX["K"] - qnorm(1-sig.level/2)*ratio.YX["se(K)"]
              upperCIyx.cdRatio <- ratio.YX["K"] + qnorm(1-sig.level/2)*ratio.YX["se(K)"]
              
              lowerCIxy.cdRatio <- ratio.XY["K"] - qnorm(1-sig.level/2)*ratio.XY["se(K)"]
              upperCIxy.cdRatio <- ratio.XY["K"] + qnorm(1-sig.level/2)*ratio.XY["se(K)"]
              
              # CD-Egger CIs
              lowerCIyx.cdEgger <- egger.YX["K"] - qnorm(1-sig.level/2)*egger.YX["se(K)"]
              upperCIyx.cdEgger <- egger.YX["K"] + qnorm(1-sig.level/2)*egger.YX["se(K)"]
              
              lowerCIxy.cdEgger <- egger.XY["K"] - qnorm(1-sig.level/2)*egger.XY["se(K)"]
              upperCIxy.cdEgger <- egger.XY["K"] + qnorm(1-sig.level/2)*egger.XY["se(K)"]
              
              # CD-GLS CIs
              lowerCIyx.cdGLS <- gls.YX["K"] - qnorm(1-sig.level/2)*gls.YX["se(K)"]
              upperCIyx.cdGLS <- gls.YX["K"] + qnorm(1-sig.level/2)*gls.YX["se(K)"]
              
              lowerCIxy.cdGLS <- gls.XY["K"] - qnorm(1-sig.level/2)*gls.XY["se(K)"]
              upperCIxy.cdGLS <- gls.XY["K"] + qnorm(1-sig.level/2)*gls.XY["se(K)"]
            
              # case 1: X->Y if CIxy not inside [-1,1] and CIyx inside [-1,0) OR (0,1]
              # case 2: X<-Y if CIyx not inside [-1,1] and CIxy inside [-1,0) OR (0,1]
              # case 3: neither otherwise
              # CD-Ratio decisions
              if(upperCIxy.cdRatio < lowerCIxy.cdRatio){
                stop("LowerCIxy > UpperCIxy")}
              if(upperCIyx.cdRatio < lowerCIyx.cdRatio){
                stop("LowerCIyx > UpperCIyx")}
              
              if(upperCIxy.cdEgger < lowerCIxy.cdEgger){
                stop("LowerCIxyEgger > UpperCIxyEgger")}
              if(upperCIyx.cdEgger < lowerCIyx.cdEgger){
                stop("LowerCIyxEgger > UpperCIyxEgger")}
              
              if(upperCIxy.cdGLS < lowerCIxy.cdGLS){
                stop("LowerCIxy > UpperCIxy")}
              if(upperCIyx.cdGLS < lowerCIyx.cdGLS){
                stop("LowerCIyx > UpperCIyx")}
              
              if("CDRatio" %in% colnames.methods){
                  if((upperCIxy.cdRatio < (-1) | lowerCIxy.cdRatio>1) &
                     ((lowerCIyx.cdRatio>=(-1) & upperCIyx.cdRatio<0) | (lowerCIyx.cdRatio>(0) & upperCIyx.cdRatio<=1)) ){
                    mat1[bX, "CDRatio"] <- mat1[bX, "CDRatio"] + 1}
                  if((upperCIyx.cdRatio < (-1) | lowerCIyx.cdRatio>1) &
                     ((lowerCIxy.cdRatio>=(-1) & upperCIxy.cdRatio< 0) | (lowerCIxy.cdRatio>(0) & upperCIxy.cdRatio<= 1)) ){
                    mat2[bX, "CDRatio"] <- mat2[bX, "CDRatio"] + 1}
              }
            
              # CD-Egger decisions
              if("CDEgger" %in% colnames.methods){
                  if((upperCIxy.cdEgger < (-1)  | lowerCIxy.cdEgger>1) &
                     ((lowerCIyx.cdEgger>=(-1) & upperCIyx.cdEgger<0) | (lowerCIyx.cdEgger>(0) & upperCIyx.cdEgger<=1))){
                    mat1[bX, "CDEgger"] <- mat1[bX, "CDEgger"] + 1}
                  if((upperCIyx.cdEgger < (-1) | lowerCIyx.cdEgger>1) &
                     ((lowerCIxy.cdEgger>=(-1)  & upperCIxy.cdEgger< 0) | (lowerCIxy.cdEgger>(0)  & upperCIxy.cdEgger<=1))){
                    mat2[bX, "CDEgger"] <- mat2[bX, "CDEgger"] + 1}
              }
              
              # CD-GLS decisions
              if("CDgls" %in% colnames.methods){
                  if((upperCIxy.cdGLS < (-1) | lowerCIxy.cdGLS>1) &
                     ((lowerCIyx.cdGLS>= (-1) & upperCIyx.cdGLS<0) | (lowerCIyx.cdGLS> (0) & upperCIyx.cdGLS<=1)) ){
                    mat1[bX, "CDgls"] <- mat1[bX, "CDgls"] + 1}
                  if((upperCIyx.cdGLS < (-1) | lowerCIyx.cdGLS>1) &
                     ((lowerCIxy.cdGLS>= (-1) & upperCIxy.cdGLS<0) | (lowerCIxy.cdGLS> (0) & upperCIxy.cdGLS<=1)) ){
                    mat2[bX, "CDgls"] <- mat2[bX, "CDgls"] + 1}
              }
            }
       }
        
        ################################################################################
        # Start BiDirectCausal methods (from https://github.com/xue-hr/BiDirectCausal/)
        # "BDCDcML.S.DP","BDMRcML.S.DP", "BDCD.Ratio.S", "BDCD.Egger.S"
        ################################################################################
        
        if( ("BDCDcML.S.DP" %in% colnames.methods) | ("BDMRcML.S.DP" %in% colnames.methods) ){
        
            # Methods error out if no snps after screening - check that SNPs remain
            sig.cutoff1 = 0.05/(ncol(GX1)+ncol(GY1))
            # taken from https://github.com/xue-hr/BiDirectCausal -- .R files
            # Check if indices or empty
            pvalue.X1 = pnorm(-abs(c(betaGX_X, betaGY_X)/c(seGX_X, seGY_X))) * 2
            pvalue.Y1 = pnorm(-abs(c(betaGX_Y, betaGY_Y)/c(seGX_Y, seGY_Y))) * 2
            ind_X1 = which(pvalue.X1 < (sig.cutoff1))
            ind_Y1 = which(pvalue.Y1 < (sig.cutoff1))
            cor_X1 = c(betaGX_X, betaGY_X)/sqrt(c(betaGX_X, betaGY_X)^2 + (nX - 2) * c(seGX_X, seGY_X)^2)
            cor_Y1 = c(betaGX_Y, betaGY_Y)/sqrt(c(betaGX_Y, betaGY_Y)^2 + (nY - 2) * c(seGX_Y, seGY_Y)^2)
            se_corX1 = sqrt((1 - cor_X1^2)^2/nX)
            se_corY1 = sqrt((1 - cor_Y1^2)^2/nY)
            intersect.ind.X.Y1 = intersect(ind_X1, ind_Y1)
            
            ind_X_new1 = setdiff(ind_X1, intersect.ind.X.Y1[(abs(cor_X1)[intersect.ind.X.Y1]) < (abs(cor_Y1)[intersect.ind.X.Y1])])
            ind_Y_new1 = setdiff(ind_Y1, intersect.ind.X.Y1[(abs(cor_X1)[intersect.ind.X.Y1]) > (abs(cor_Y1)[intersect.ind.X.Y1])])
            
            ## Check first! If true, case 3!
            if(identical(ind_X_new1, integer(0)) | identical(ind_Y_new1, integer(0)) ){
                if("BDCDcML.S.DP" %in% colnames.methods){
                    mat3[bX,"BDCDcML.S.DP"] <- mat3[bX,"BDCDcML.S.DP"]+1
                    matError[bX, "BDCDcML.S.DP"] <- matError[bX, "BDCDcML.S.DP"]+1
                }
                if( "BDMRcML.S.DP" %in% colnames.methods ){
                    mat3[bX,"BDMRcML.S.DP"] <- mat3[bX,"BDMRcML.S.DP"]+1
                    matError[bX, "BDMRcML.S.DP"] <- matError[bX, "BDMRcML.S.DP"]+1
                }
              
            }else{
              ###############################################################################
              # Apply bi-directional CDcML methods: .S.DP
              ###############################################################################
              
              if("BDCDcML.S.DP" %in% colnames.methods){
                  CDcML.out=BiDirCDcML(b_X= c(betaGX_X, betaGY_X),
                                       b_Y=c(betaGX_Y, betaGY_Y),
                                       se_X= c(seGX_X, seGY_X),
                                       se_Y= c(seGX_Y, seGY_Y),
                                       n_X=nX,
                                       n_Y=nY,
                                       sig.cutoff = sig.cutoff1,
                                       num_pert = 100)
                  
                  
                  # CIs - CDcML.S.DP
                  ci.low.XY.S.DP <- CDcML.out$XtoY.est.S.DP - (qnorm(0.975)*CDcML.out$XtoY.se.S.DP)
                  ci.upp.XY.S.DP <- CDcML.out$XtoY.est.S.DP + (qnorm(0.975)*CDcML.out$XtoY.se.S.DP)
                  
                  ci.low.YX.S.DP <- CDcML.out$YtoX.est.S.DP - (qnorm(0.975)*CDcML.out$YtoX.se.S.DP)
                  ci.upp.YX.S.DP <- CDcML.out$YtoX.est.S.DP + (qnorm(0.975)*CDcML.out$YtoX.se.S.DP)
                  
                  if(ci.upp.XY.S.DP < ci.low.XY.S.DP){
                    stop("CDcML: LowerCIxy > UpperCIxy")}
                  if(ci.upp.YX.S.DP < ci.low.YX.S.DP){
                    stop("CDcML: LowerCIyx > UpperCIyx")}
                  
                  # decision rules with CIs  --- .S.DP
                  if((ci.low.XY.S.DP>0 & ci.upp.XY.S.DP<1)| (ci.low.XY.S.DP>-1 & ci.upp.XY.S.DP<0)){
                    if((ci.low.YX.S.DP>0 & ci.upp.YX.S.DP<1) | (ci.low.YX.S.DP>-1 & ci.upp.YX.S.DP<0)){mat3[bX,"BDCDcML.S.DP"] <- mat3[bX,"BDCDcML.S.DP"]+1
                    }else{mat1[bX,"BDCDcML.S.DP"] <- mat1[bX,"BDCDcML.S.DP"]+1}
                  }else if((ci.low.YX.S.DP>0 & ci.upp.YX.S.DP<1) | (ci.low.YX.S.DP>-1 & ci.upp.YX.S.DP<0)){
                    if((ci.low.XY.S.DP>0 & ci.upp.XY.S.DP<1)| (ci.low.XY.S.DP>-1 & ci.upp.XY.S.DP<0)){mat3[bX,"BDCDcML.S.DP"] <- mat3[bX,"BDCDcML.S.DP"]+1
                    }else{mat2[bX,"BDCDcML.S.DP"] <- mat2[bX,"BDCDcML.S.DP"]+1}
                  }else{mat3[bX,"BDCDcML.S.DP"] <- mat3[bX,"BDCDcML.S.DP"]+1}
             }
             
             if( "BDMRcML.S.DP" %in% colnames.methods ){
                MRcML.out <- BiDirMRcML(b_X=c(betaGX_X, betaGY_X),
                                          b_Y=c(betaGX_Y, betaGY_Y),
                                          se_X=c(seGX_X, seGY_X),
                                          se_Y=c(seGX_Y, seGY_Y),
                                          n_X=nX,
                                          n_Y=nY,
                                          sig.cutoff = sig.cutoff1,
                                          num_pert = 100)
                  
                  
                  p.XY.BDMRcML.S.DP<- pnorm(-abs(MRcML.out$XtoY.est.S.DP / MRcML.out$XtoY.se.S.DP))*2
                  p.YX.BDMRcML.S.DP<- pnorm(-abs(MRcML.out$YtoX.est.S.DP / MRcML.out$YtoX.se.S.DP))*2
                  
                  # decision rules
                  if(p.XY.BDMRcML.S.DP <= 0.05 & p.YX.BDMRcML.S.DP > 0.05){mat1[bX,"BDMRcML.S.DP"] <- mat1[bX,"BDMRcML.S.DP"]+1}#return case 1
                  if(p.XY.BDMRcML.S.DP > 0.05 & p.YX.BDMRcML.S.DP <= 0.05){mat2[bX,"BDMRcML.S.DP"] <- mat2[bX,"BDMRcML.S.DP"]+1}# return case 2
                  if((p.XY.BDMRcML.S.DP >0.05 & p.YX.BDMRcML.S.DP > 0.05) | (p.XY.BDMRcML.S.DP <=0.05 & p.YX.BDMRcML.S.DP <= 0.05)){mat3[bX,"BDMRcML.S.DP"] <- mat3[bX,"BDMRcML.S.DP"]+1}# return case 3
              }
                  
            }
        }
        
        # CD-Ratio screening, CD-Egger screening
        if( ("BDCD.Ratio.S" %in% colnames.methods) | ("BDCD.Egger.S" %in% colnames.methods)){
            # methods error out if only 0-1 SNPs included after screening, before running BiDirCDMethods run this check
            
            pvalue.X = pnorm(-abs(b_X/se_X))*2
            pvalue.Y = pnorm(-abs(b_Y/se_Y))*2
            
            sig.cutoff1=0.05/(ncol(GX1)+ncol(GY1))
            
            ind_X = which(pvalue.X<(sig.cutoff1))
            ind_Y = which(pvalue.Y<(sig.cutoff1))
            
            cor_X = b_X / sqrt(b_X^2 + (nX-2)*se_X^2)
            cor_Y = b_Y / sqrt(b_Y^2 + (nY-2)*se_Y^2)
            
            se_corX = sqrt((1-cor_X^2)^2/nX)
            se_corY = sqrt((1-cor_Y^2)^2/nY)
            
            intersect.ind.X.Y = intersect(ind_X,ind_Y)
            ind_X_new = setdiff(ind_X,
                                intersect.ind.X.Y[(abs(cor_X)[intersect.ind.X.Y])<
                                                    (abs(cor_Y)[intersect.ind.X.Y])])
            ind_Y_new = setdiff(ind_Y,
                                intersect.ind.X.Y[(abs(cor_X)[intersect.ind.X.Y])>
                                                    (abs(cor_Y)[intersect.ind.X.Y])])
            
            if(identical(ind_X, integer(0)) | identical(ind_Y, integer(0)) | identical(ind_X_new, integer(0)) | identical(ind_Y_new, integer(0)) |
               length(ind_X_new) == 1 | length(ind_Y_new)==1){
              
              if("BDCD.Ratio.S" %in% colnames.methods){
                mat3[bX,"BDCD.Ratio.S"] <- mat3[bX,"BDCD.Ratio.S"]+1
                matError[bX, "BDCD.Ratio.S"] <- matError[bX, "BDCD.Ratio.S"]+1
                
              }
              if("BDCD.Egger.S" %in% colnames.methods){
                mat3[bX,"BDCD.Egger.S"] <- mat3[bX,"BDCD.Egger.S"]+1
                matError[bX, "BDCD.Egger.S"] <- matError[bX, "BDCD.Egger.S"]+1
              }
                
            }else{
              
              CD.Out <- BiDirCDMethod(b_X = c(betaGX_X, betaGY_X),
                                      b_Y = c(betaGX_Y, betaGY_Y),
                                      se_X = c(seGX_X, seGY_X),
                                      se_Y = c(seGX_Y, seGY_Y),
                                      n_X = nX,
                                      n_Y = nY,
                                      sig.cutoff = 0.05/(ncol(GX1)+ncol(GY1)), random.seed = 0)
              
              K.ratio.XY.S = CD.Out$CDRatio.XtoY.est.S
              
              seK.ratio.XY.S = CD.Out$CDRatio.XtoY.se.S
              
              K.egger.XY.S = CD.Out$CDEgger.XtoY.est.S
              
              seK.egger.XY.S =CD.Out$CDEgger.XtoY.se.S
              
              K.ratio.YX.S = CD.Out$CDRatio.YtoX.est.S
              
              seK.ratio.YX.S= CD.Out$CDRatio.YtoX.se.S
              
              K.egger.YX.S = CD.Out$CDEgger.YtoX.est.S
              
              seK.egger.YX.S = CD.Out$CDEgger.YtoX.se.S
              
              # create 1-alpha CIs for K_XY and K_YX
              ci.low.ratioXY.S <- K.ratio.XY.S - (qnorm(0.975)*seK.ratio.XY.S)
              ci.upp.ratioXY.S <- K.ratio.XY.S + (qnorm(0.975)*seK.ratio.XY.S)
              
              ci.low.ratioYX.S <- K.ratio.YX.S - (qnorm(0.975)*seK.ratio.YX.S)
              ci.upp.ratioYX.S <- K.ratio.YX.S + (qnorm(0.975)*seK.ratio.YX.S)
              
              ci.low.eggerXY.S <- K.egger.XY.S - (qnorm(0.975)*seK.egger.XY.S)
              ci.upp.eggerXY.S <- K.egger.XY.S + (qnorm(0.975)*seK.egger.XY.S)
              
              ci.low.eggerYX.S <- K.egger.YX.S - (qnorm(0.975)*seK.egger.YX.S)
              ci.upp.eggerYX.S <- K.egger.YX.S + (qnorm(0.975)*seK.egger.YX.S)
              
              if(ci.upp.ratioXY.S < ci.low.ratioXY.S){
                stop("CD-Ratio: LowerCIxy > UpperCIxy")}
              if(ci.upp.ratioYX.S < ci.low.ratioYX.S){
                stop("CD-Ratio: LowerCIyx > UpperCIyx")}
              
              if(ci.upp.eggerXY.S < ci.low.eggerXY.S){
                stop("CD-Egger: LowerCIxy > UpperCIxy")}
              if(ci.upp.eggerYX.S < ci.low.eggerYX.S){
                stop("CD-Egger: LowerCIyx > UpperCIyx")}
              
              #  DECISION RULES CD Ratio - .S
              
              if("BDCD.Ratio.S" %in% colnames.methods){
                if((ci.low.ratioXY.S>0 & ci.upp.ratioXY.S<1)| (ci.low.ratioXY.S>-1 & ci.upp.ratioXY.S<0)){
                  if((ci.low.ratioYX.S>0 & ci.upp.ratioYX.S<1) | (ci.low.ratioYX.S>-1 & ci.upp.ratioYX.S<0)){ mat3[bX,"BDCD.Ratio.S"] <- mat3[bX,"BDCD.Ratio.S"]+1
                  }else{mat1[bX,"BDCD.Ratio.S"] <- mat1[bX,"BDCD.Ratio.S"]+1}
                }else if((ci.low.ratioYX.S>0 & ci.upp.ratioYX.S<1) | (ci.low.ratioYX.S>-1 & ci.upp.ratioYX.S<0)){
                  if((ci.low.ratioXY.S>0 & ci.upp.ratioXY.S<1)| (ci.low.ratioXY.S>-1 & ci.upp.ratioXY.S<0)){ mat3[bX,"BDCD.Ratio.S"] <- mat3[bX,"BDCD.Ratio.S"]+1
                  }else{mat2[bX,"BDCD.Ratio.S"] <- mat2[bX,"BDCD.Ratio.S"]+1}
                }else{ mat3[bX,"BDCD.Ratio.S"] <- mat3[bX,"BDCD.Ratio.S"]+1}
              }
              
              #  DECISION RULES CD Egger
              if("BDCD.Egger.S" %in% colnames.methods){
                
                if((ci.low.eggerXY.S>0 & ci.upp.eggerXY.S<1)| (ci.low.eggerXY.S>-1 & ci.upp.eggerXY.S<0)){
                  if((ci.low.eggerYX.S>0 & ci.upp.eggerYX.S<1) | (ci.low.eggerYX.S>-1 & ci.upp.eggerYX.S<0)){mat3[bX,"BDCD.Egger.S"] <- mat3[bX,"BDCD.Egger.S"]+1
                  }else{mat1[bX,"BDCD.Egger.S"] <- mat1[bX,"BDCD.Egger.S"]+1}
                }else if((ci.low.eggerYX.S>0 & ci.upp.eggerYX.S<1) | (ci.low.eggerYX.S>-1 & ci.upp.eggerYX.S<0)){
                  if((ci.low.eggerXY.S>0 & ci.upp.eggerXY.S<1)| (ci.low.eggerXY.S>-1 & ci.upp.eggerXY.S<0)){mat3[bX,"BDCD.Egger.S"] <- mat3[bX,"BDCD.Egger.S"]+1
                  }else{mat2[bX,"BDCD.Egger.S"] <- mat2[bX,"BDCD.Egger.S"]+1}
                }else{mat3[bX,"BDCD.Egger.S"] <- mat3[bX,"BDCD.Egger.S"]+1}
              }
              
            }
         }
        
        ################################################################################
        # Start MR Steiger - Get values and p-values from Y~G and X~G
        ################################################################################
        if( ("MRS.wMedian" %in% colnames.methods) | ("MRS.ivw" %in% colnames.methods) | ("MRS.Egger" %in% colnames.methods) ){
        
            #Vector of p-values of SNP-exposure
            p_exp<-rep(0,nSNP_GX)
            for(pe in 1:nSNP_GX){
              if(contX==T){p_exp[pe]<-summary(lm(X~GX1[,pe]))$coef[2,4]}
              if(contX==F){p_exp[pe]<-summary(glm(X~GX1[,pe], family="binomial"))$coef[2,4]}
            }
            
            #Vector of p-values of SNP-outcome
            p_out<-rep(0,nSNP_GX)
            for(po in 1:nSNP_GX){
              if(contY==T){p_out[po]<-summary(lm(Y1~GX2[,po]))$coef[2,4]}
              if(contY==F){p_out[po]<-summary(glm(Y1~GX2[,po], family="binomial"))$coef[2,4]}
              
            }
            
            #Sample sizes for p_exp, p_out
            n_exp<-nX
            n_out<-nY
            
            #Vector of absolute correlations for SNP-exposure
            r_exp<-rep(0,nSNP_GX)
            
            for(re in 1:nSNP_GX){
              if(contX==T){r_exp[re]<-abs(cor(GX1[,re],X))}
              if(contX==F){
                beta <- summary(glm(X~GX1[,re], family="binomial"))$coef[2,1]
                maf <- mean(GX1[,re], na.rm=T)/2
                ncases <- length(X[X==1])
                ncontrols <- length(X[X==0])
                prev<- (ncases)/(ncases+ncontrols)
                
                r_exp[re] <- get_r_from_lor(lor=beta, af=maf, ncase=ncases, ncontrol=ncontrols, prevalence = prev, model="logit")}
            }
            
            #Vector of absolute correlations for SNP-outcome
            r_out<-rep(0,nSNP_GX)
            
            for(ro in 1:nSNP_GX){
              if(contY==T){r_out[ro]<-abs(cor(GX2[,ro],Y1))}
              if(contY==F){
                beta <- summary(glm(Y1~GX2[,ro], family="binomial"))$coef[2,1]
                maf <- mean(GX2[,re], na.rm=T)/2
                ncases <- length(Y1[Y1==1])
                ncontrols <- length(Y1[Y1==0])
                prev<- (ncases)/(ncases+ncontrols)
                
                r_out[ro] <- get_r_from_lor(lor=beta, af=maf, ncase=ncases, ncontrol=ncontrols, prevalence = prev, model="logit")}
            }
            
            #r_xxo Measurememt precision of exposure, default 1
            #r_yyo Measurement precision of outcome, default 1
            
            
            ################################################################################
            # MR Steiger
            ################################################################################
            
            #A statistical test for whether the assumption that exposure causes outcome is valid
            mtryMRS <- try(mrs<-mr_steiger(p_exp, p_out, n_exp, n_out, r_exp, r_out, r_xxo = 1, r_yyo = 1))
            
            if( ("MRS.wMedian" %in% colnames.methods) ){
              # matError if MR Steiger errors out
              if(inherits(mtryMRS, "try-error") | is.na(mtryMRS$steiger_test)){
                matError[bX, "MRS.wMedian"] <- matError[bX, "MRS.wMedian"]+1
              }
              # return cases if no error returned
              if(!inherits(mtryMRS, "try-error") & !is.na(mtryMRS$steiger_test)){
                # case 1: X->Y1 if steiger_test<alpha and weighted median pval<alpha and correct causal direction == T
                # case 2: X<-Y1 if correct causal direction == F & weighted median pval < alpha & steiger_test<alpha
                if(mrs$correct_causal_direction==TRUE & mrs$steiger_test <= sig.level & tsmr.weighted_median.GX$pval <= sig.level){
                  mat1[bX,"MRS.wMedian"]<-mat1[bX,"MRS.wMedian"]+1}
                if(mrs$correct_causal_direction==FALSE & mrs$steiger_test <= sig.level & tsmr.weighted_median.GX$pval <= sig.level){
                  mat2[bX,"MRS.wMedian"]<-mat2[bX,"MRS.wMedian"]+1}
              }
            }
            
            if( ("MRS.ivw" %in% colnames.methods)){
              # matError if MR Steiger errors out
              if(inherits(mtryMRS, "try-error") | is.na(mtryMRS$steiger_test)){
                matError[bX, "MRS.ivw"] <- matError[bX, "MRS.ivw"]+1
              }
              if(!inherits(mtryMRS, "try-error") & !is.na(mtryMRS$steiger_test)){
                # case 1: X->Y1 if steiger_test<alpha and MR IVW pval<alpha and correct causal direction == T
                # case 2: X<-Y1 if correct causal direction == F & MR IVW pval < alpha & steiger_test<alpha
                if(mrs$correct_causal_direction==TRUE & mrs$steiger_test <= sig.level & tsmr.IVW.GX$pval <= sig.level){
                  mat1[bX,"MRS.ivw"]<-mat1[bX,"MRS.ivw"]+1}
                if(mrs$correct_causal_direction==FALSE & mrs$steiger_test <= sig.level & tsmr.IVW.GX$pval <= sig.level){
                  mat2[bX,"MRS.ivw"]<-mat2[bX,"MRS.ivw"]+1}
              }
           }
            
            if( ("MRS.Egger" %in% colnames.methods)){
              # matError if MR Steiger errors out
              if(inherits(mtryMRS, "try-error") | is.na(mtryMRS$steiger_test)){
                matError[bX, "MRS.Egger"] <- matError[bX, "MRS.Egger"]+1
              }
              if(!inherits(mtryMRS, "try-error") & !is.na(mtryMRS$steiger_test)){
                # case 1: X->Y1 if steiger_test<alpha and MR Egger pval<alpha and correct causal direction == T
                # case 2: X<-Y1 if correct causal direction == F & MR Egger pval < alpha & steiger_test<alpha
                if(mrs$correct_causal_direction==TRUE & mrs$steiger_test <= sig.level & tsmr.Egger.GX$pval <= sig.level){
                  mat1[bX,"MRS.Egger"]<-mat1[bX,"MRS.Egger"]+1}
                if(mrs$correct_causal_direction==FALSE & mrs$steiger_test <= sig.level & tsmr.Egger.GX$pval <= sig.level){
                  mat2[bX,"MRS.Egger"]<-mat2[bX,"MRS.Egger"]+1}
              }
            }
        
        }
        
      ################################################################################
        # end loops
        ################################################################################    
      }#beta loop
    }#sim loop
    
    pl=F
    if(any(gammaGY!=0) | any(betaGX!=0)){pl=T}
    
    ##############################
    # start plots
    ##############################
    mat_total1 <- as.data.frame(cbind(mat1[,1],(mat1[,-1]/nSim)))
    mat_total2 <- as.data.frame(cbind(mat2[,1],(mat2[,-1]/nSim)))
    mat_totalE <- as.data.frame(cbind(matError[,1],(matError[,-1]/nSim)))
    
    colnames(mat_total1) <- colnames(mat1)
    colnames(mat_total2) <- colnames(mat2)
    colnames(mat_totalE) <- colnames(matError)
    
    
    if(plot.pdf){
        #########################################
        # case 1 plot
        #########################################
        pdf(paste(path1,"Plot_Case1nSim",nSim,"seed",SEED ,"nX",nX,"nY",nY,"nSNPX", length(MAF_GX),"nSNPY", length(MAF_GY), "contX", contX, "contY",contY,"u",unmeasuredConfounding, "me",measurementError,"p",pl, "l1", long1, "l2",long2,"dx", deltaX, "bU", betaU, "gU",gammaU,"bX",betaGX[1]*100, "gY", gammaGY[1]*100,"eGX", etaGX[1]*100, "eGY", etaGY[1]*100,".pdf", sep = ""))
        
        plot(-2,-2,xlim=c(min(betaX),max(betaX)),ylim=c(0,1.1),main="",xlab=expression(beta[X]),ylab="Proportion of Simulations")
        
        if( ("MRS.ivw" %in% colnames.methods)){lines(betaX,mat_total1[,"MRS.ivw"],col="blue4",pch=1,lty=2,type="b", lwd=2.4)}
        if( ("MRS.wMedian" %in% colnames.methods)){lines(betaX,mat_total1[,"MRS.wMedian"],col="steelblue1",pch=3,lty=3,type="b", lwd=2.4)}
        if( ("MRS.Egger" %in% colnames.methods)){lines(betaX,mat_total1[,"MRS.Egger"],col="darkslategray2",pch=4,lty=4,type="b", lwd=2.4)}
        
        # CD methods - pch=2; PURPLES!!! -- dark purple, bright purple, medium purple
        if( ("CDRatio" %in% colnames.methods)){lines(betaX,mat_total1[,"CDRatio"],col="darkorange2",pch=1,lty=1,type="b", lwd=2.4)}
        if( ("CDEgger" %in% colnames.methods)){lines(betaX,mat_total1[,"CDEgger"],col="darkorange1",pch=3,lty=2,type="b", lwd=2.4)}
        if( ("CDgls" %in% colnames.methods)){lines(betaX,mat_total1[,"CDgls"],col="orange2",pch=4,lty=3,type="b", lwd=2.4)}
        
        # Bi-directional methods - reds/pinks for TSMR, reds/browns for MR
        # IVW methods - pch=4; REDS
        if( ("tsmr_IVW" %in% colnames.methods)){lines(betaX,mat_total1[,"tsmr_IVW"],col="plum2",pch=3,lty=1,type="b", lwd=2.4)}
        if( ("tsmr_IVW_mre" %in% colnames.methods)){lines(betaX,mat_total1[,"tsmr_IVW_mre"],col="mediumpurple3",pch=4,lty=4,type="b", lwd=2.4)}
        if( ("tsmr_IVW_fe" %in% colnames.methods)){lines(betaX,mat_total1[,"tsmr_IVW_fe"],col="mediumorchid1",pch=6,lty=1,type="b", lwd=2.4)}
        ##########
        # Egger methods - pch = 0
        if( ("tsmr_Egger" %in% colnames.methods)){lines(betaX,mat_total1[,"tsmr_Egger"],col="deeppink2",pch=1,lty=1,type="b", lwd=2.4)}
        if( ("tsmr_Egger_Boot" %in% colnames.methods)){lines(betaX,mat_total1[,"tsmr_Egger_Boot"],col="firebrick1",pch=4,lty=2,type="b", lwd=2.4)}
        
        # median methods - pch=6 & 8
        if( ("tsmr_weighted_median" %in% colnames.methods)){lines(betaX,mat_total1[,"tsmr_weighted_median"],col="darkgreen",pch=4,lty=3,type="b", lwd=2.4)}
        if( ("tsmr_simple_median" %in% colnames.methods)){lines(betaX,mat_total1[,"tsmr_simple_median"],col="chartreuse2",pch=6,lty=4,type="b", lwd=2.4)}
        if( ("tsmr_pwm" %in% colnames.methods)){lines(betaX,mat_total1[,"tsmr_pwm"],col="lightgreen",pch=0,lty=5,type="b", lwd=2.4)}
        
        # other tsmr methods - pch=7
        if( ("tsmr_uwr" %in% colnames.methods)){lines(betaX,mat_total1[,"tsmr_uwr"],col="black",pch=1,lty=1,type="b", lwd=2.4)}
        
        if( ("BDCDcML.S.DP" %in% colnames.methods)){lines(betaX,mat_total1[,"BDCDcML.S.DP"],col="goldenrod3",pch=2,lty=2,type="b", lwd=2.4)}
        if( ("BDMRcML.S.DP" %in% colnames.methods)){lines(betaX,mat_total1[,"BDMRcML.S.DP"],col="goldenrod1",pch=3,lty=1,type="b", lwd=2.4)}
        if( ("BDCD.Ratio.S" %in% colnames.methods)){lines(betaX,mat_total1[,"BDCD.Ratio.S"],col="yellow4",pch=6,lty=1,type="b", lwd=2.4)}
        if( ("BDCD.Egger.S" %in% colnames.methods)){lines(betaX,mat_total1[,"BDCD.Egger.S"],col="yellow3",pch=0,lty=4,type="b", lwd=2.4)}
        
        dev.off()
        
        #########################################
        # case 2 plot
        #########################################
        pdf(paste(path1,"Plot_Case2nSim",nSim,"seed",SEED ,"nX",nX,"nY",nY,"nSNPX", length(MAF_GX),"nSNPY", length(MAF_GY), "contX", contX, "contY",contY,"u",unmeasuredConfounding, "me",measurementError,"p",pl, "l1", long1, "l2",long2,"dx", deltaX, "bU", betaU, "gU",gammaU,"bX",betaGX[1]*100, "gY", gammaGY[1]*100,"eGX", etaGX[1]*100, "eGY", etaGY[1]*100,".pdf", sep = ""))
        plot(-2,-2,xlim=c(min(betaX),max(betaX)),ylim=c(0,1.1),main="",xlab=expression(beta[X]),ylab="Proportion of Simulations")
        ###########
        
        if( ("MRS.ivw" %in% colnames.methods)){lines(betaX,mat_total2[,"MRS.ivw"],col="blue4",pch=1,lty=2,type="b", lwd=2.4)}
        if( ("MRS.wMedian" %in% colnames.methods)){lines(betaX,mat_total2[,"MRS.wMedian"],col="steelblue1",pch=3,lty=3,type="b", lwd=2.4)}
        if( ("MRS.Egger" %in% colnames.methods)){lines(betaX,mat_total2[,"MRS.Egger"],col="darkslategray2",pch=4,lty=4,type="b", lwd=2.4)}
        
        # CD methods - pch=2; PURPLES!!! -- dark purple, bright purple, medium purple
        if( ("CDRatio" %in% colnames.methods)){lines(betaX,mat_total2[,"CDRatio"],col="darkorange2",pch=1,lty=1,type="b", lwd=2.4)}
        if( ("CDEgger" %in% colnames.methods)){lines(betaX,mat_total2[,"CDEgger"],col="darkorange1",pch=3,lty=2,type="b", lwd=2.4)}
        if( ("CDgls" %in% colnames.methods)){lines(betaX,mat_total2[,"CDgls"],col="orange2",pch=4,lty=3,type="b", lwd=2.4) }
        
        # Bi-directional methods - reds/pinks for TSMR, reds/browns for MR
        # IVW methods - pch=4; REDS
        if( ("tsmr_IVW" %in% colnames.methods)){lines(betaX,mat_total2[,"tsmr_IVW"],col="plum2",pch=3,lty=1,type="b", lwd=2.4)}
        if( ("tsmr_IVW_mre" %in% colnames.methods)){lines(betaX,mat_total2[,"tsmr_IVW_mre"],col="mediumpurple3",pch=4,lty=4,type="b", lwd=2.4)}
        if( ("tsmr_IVW_fe" %in% colnames.methods)){lines(betaX,mat_total2[,"tsmr_IVW_fe"],col="mediumorchid1",pch=6,lty=1,type="b", lwd=2.4)}
        ##########
        # Egger methods - pch = 0
        if( ("tsmr_Egger" %in% colnames.methods)){lines(betaX,mat_total2[,"tsmr_Egger"],col="deeppink2",pch=1,lty=1,type="b", lwd=2.4)}
        if( ("tsmr_Egger_Boot" %in% colnames.methods)){lines(betaX,mat_total2[,"tsmr_Egger_Boot"],col="firebrick1",pch=4,lty=2,type="b", lwd=2.4)}
        
        # median methods - pch=6 & 8
        if( ("tsmr_weighted_median" %in% colnames.methods)){lines(betaX,mat_total2[,"tsmr_weighted_median"],col="darkgreen",pch=4,lty=3,type="b", lwd=2.4)}
        if( ("tsmr_simple_median" %in% colnames.methods)){lines(betaX,mat_total2[,"tsmr_simple_median"],col="chartreuse2",pch=6,lty=4,type="b", lwd=2.4)}
        if( ("tsmr_pwm" %in% colnames.methods)){lines(betaX,mat_total2[,"tsmr_pwm"],col="lightgreen",pch=0,lty=5,type="b", lwd=2.4)}
        
        # other tsmr methods - pch=7
        if( ("tsmr_uwr" %in% colnames.methods)){lines(betaX,mat_total2[,"tsmr_uwr"],col="black",pch=1,lty=1,type="b", lwd=2.4)}
        
        if( ("BDCDcML.S.DP" %in% colnames.methods)){lines(betaX,mat_total2[,"BDCDcML.S.DP"],col="goldenrod3",pch=2,lty=2,type="b", lwd=2.4)}
        if( ("BDMRcML.S.DP" %in% colnames.methods)){lines(betaX,mat_total2[,"BDMRcML.S.DP"],col="goldenrod1",pch=3,lty=1,type="b", lwd=2.4)}
          
        if( ("BDCD.Ratio.S" %in% colnames.methods)){lines(betaX,mat_total2[,"BDCD.Ratio.S"],col="yellow4",pch=6,lty=1,type="b", lwd=2.4)}
        if( ("BDCD.Egger.S" %in% colnames.methods)){lines(betaX,mat_total2[,"BDCD.Egger.S"],col="yellow3",pch=0,lty=4,type="b", lwd=2.4)}
        
        dev.off()
        
        #########################################
        # legend
        #########################################
        pdf(paste(path1,"Legend_nSim",nSim,"seed",SEED ,"nX",nX,"nY",nY,"nSNPX", length(MAF_GX),"nSNPY", length(MAF_GY), "contX", contX, "contY",contY,"u",unmeasuredConfounding, "me",measurementError,"p",pl, "l1", long1, "l2",long2,"dx", deltaX, "bU", betaU, "gU",gammaU,"bX",betaGX[1]*100, "gY", gammaGY[1]*100,"eGX", etaGX[1]*100, "eGY", etaGY[1]*100,".pdf", sep = ""))
        plot(NULL, xlim=0:1, ylim=0:1, ylab='', xlab='', xaxt='n', yaxt='n', bty='n')
        legend("topleft", xpd=T,legend = c("MR Steiger IVW", "MR Steiger weighted median", "MR Steiger Egger",
                                           "CD-Ratio", "CD-Egger", "CD-GLS",
                                           "CD-cML (screening, data perturbation)", "MR-cML (screening, data perturbation)",
                                           "CD-Ratio (screening)", "CD-Egger (screening)",
                                           "Bi-MR IVW", "Bi-MR IVW MRE", "Bi-MR IVW FE",
                                           "Bi-MR Egger", "Bi-MR Egger Boot",
                                           "Bi-MR weighted median", "Bi-MR simple median", "Bi-MR penalized weighted median",
                                           "Bi-MR unweighted regression"),
               pch = c(1,3,4,
                       1,3,4,
                       2,3,6,0,
                       3,4,6,
                       1,4,
                       4,6,0,
                       1),
               
               lty = c(2,3,4,1,2,3,2,1,1,4,1,4,1,1,2,3,4,5,1),
               col = c("blue4","steelblue1","darkslategray2",
                       "darkorange2","darkorange1","orange2",
                       "goldenrod3","goldenrod1","yellow4","yellow3",
                       "plum2","mediumpurple3","mediumorchid1",
                       "deeppink2","firebrick1",
                       "darkgreen","chartreuse2","lightgreen","black"),
               cex=0.6,
               lwd=1.4)
        dev.off()
    }
    
    # save matrices
    write.table(mat_total1,file=paste0(path1,"matCase1nSim", nSim,"seed",SEED ,"nX",nX,"nY",nY,"nSNPX", length(MAF_GX),"nSNPY", length(MAF_GY), "contX", contX, "contY",contY,"u",unmeasuredConfounding, "me",measurementError,"p",pl, "l1", long1, "l2",long2,"dx", deltaX, "bU", betaU, "gU",gammaU,"bX",betaGX[1]*100, "gY", gammaGY[1]*100,"eGX", etaGX[1]*100, "eGY", etaGY[1]*100,".txt"), quote=F, row.names = F)
    write.table(mat_total2,file=paste0(path1,"matCase2nSim", nSim,"seed",SEED,"nX",nX,"nY",nY,"nSNPX", length(MAF_GX),"nSNPY", length(MAF_GY),"contX", contX, "contY",contY,"u",unmeasuredConfounding, "me",measurementError,"p",pl, "l1", long1, "l2",long2,"dx", deltaX,  "bU", betaU, "gU",gammaU,"bX",betaGX[1]*100, "gY", gammaGY[1]*100,"eGX", etaGX[1]*100, "eGY", etaGY[1]*100 ,".txt"), quote=F, row.names = F)
    write.table(mat_totalE,file=paste0(path1,"matCaseErrornSim", nSim,"seed",SEED,"nX",nX,"nY",nY,"nSNPX", length(MAF_GX),"nSNPY", length(MAF_GY),"contX", contX, "contY",contY,"u",unmeasuredConfounding, "me",measurementError,"p",pl, "l1", long1, "l2",long2,"dx", deltaX,  "bU", betaU, "gU",gammaU,"bX",betaGX[1]*100, "gY", gammaGY[1]*100,"eGX", etaGX[1]*100, "eGY", etaGY[1]*100 ,".txt"), quote=F, row.names = F)
    
    return(list("mat1"=mat_total1, "mat2"=mat_total2, "matE"=mat_totalE))
  }
