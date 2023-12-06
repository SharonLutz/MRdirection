MRdirection <-
function(nSim=100, nX=1000, nY=1000, 
           MAF_GX=rep(0.5,10), MAF_GY=rep(0.5,10),
           gamma0=0, gammaGX=rep(0.2,10), gammaGY=rep(0,10),
           measurementError=F, delta0=0, deltaX=1, varME=1,
           beta0=0, betaX=seq(0,2,by=0.5), betaGX=rep(0,10), betaGY=rep(0.2,10),
           varY=0.2, varX=1, 
           unmeasuredConfounding=F, meanU=0, varU=1, gammaU=1, betaU=1, 
           eta0=0, etaGX=rep(0.2,10), etaGY=rep(0,10), snpEffectConfounder=F, 
           independentSNPs=T, long1=F, long2=F, kappa0=0, kappaY=0.2, iota0=0,
           iotaX=0.2, sig.level=0.05, contX= TRUE, contY=TRUE,
           SEED=1, runMethods=c("All")){
    ################################################################################   
    # load libraries
    ################################################################################
    
    # load library for CD-Ratio, CD-GLS, and CD-Egger tests
    library(MRCD)
    library(TwoSampleMR)
    library(MRcML)
    library(BiDirectCausal)
    
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
    
    if(independentSNPs==F & length(MAF_GX) > ncol(MRCD::pruned$loci_bed)){stop("Number of non-independent SNPs must be 22 or less.")}
    
    ################################################################################
    # Matrix to save Results
    ################################################################################
    #save results for type 1 error rate betaX=0 and power betaX>0
    
    AllMethods <- c("MRS.ivw","MRS.wMedian","MRS.Egger","CDRatio", "CDEgger","CDgls",
                    "tsmr_IVW","tsmr_Egger",
                    "tsmr_weighted_median","tsmr_uwr","tsmr_simple_median","tsmr_pwm", 
                    "tsmr_meta_fixed_simple","tsmr_IVW_mre","tsmr_IVW_fe", "tsmr_Egger_Boot", 
                    "BDCDcML.S.DP","BDMRcML.S.DP","BDCD.Ratio.S", "BDCD.Egger.S")

    MRSteigerMethods <- c("MRS.ivw","MRS.wMedian","MRS.Egger")
    
    CDMethods <- c("CDRatio","CDEgger","CDgls","BDCDcML.S.DP","BDMRcML.S.DP","BDCD.Ratio.S", "BDCD.Egger.S")
    
    BidirectionalMethods <- c("tsmr_IVW","MR_IVW","tsmr_Egger",
                              "tsmr_weighted_median","tsmr_uwr","tsmr_simple_median","tsmr_pwm", 
                              "tsmr_meta_fixed_simple","tsmr_IVW_mre","tsmr_IVW_fe", "tsmr_Egger_Boot")
    
    
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
    
    if("tsmr_meta_fixed_simple" %in% runMethods){colnames.1 <- c(colnames.1, "tsmr_meta_fixed_simple")}
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
    colnames.methods <- colnames.1
    
    colnames(mat1) <- colnames.methods
    colnames(mat2) <- colnames.methods
    colnames(mat3) <- colnames.methods
    mat1[,"bX"] <- betaX    
    mat2[,"bX"] <- betaX    
    mat3[,"bX"] <- betaX    
    ################################################################################
    # cycle through the simulations
    ################################################################################  
    #cycle through the simulations
    for(ii in 1:nSim){
      printCut<-1
      if(floor(ii/printCut)==ceiling(ii/printCut)){print(paste(ii,"of",nSim,"simulations"))}
      set.seed(SEED+ii-1)
      ################################################################################
      # cycle through values of betaX
      ################################################################################ 
      #cycle through values of betaX
      for(bX in 1:length(betaX)){
        print(paste(bX,"of", length(betaX),"betas"))
        ################################################################################
        # generate data
        ################################################################################                            
        # Generate GX and GY (SNPs)
        nSNP_GX<-length(MAF_GX)
        nSNP_GY<-length(MAF_GY)
        GX<-matrix(0,nrow=max(nX,nY),ncol=nSNP_GX) 
        GY<-matrix(0,nrow=max(nX,nY),ncol=nSNP_GY) 
        
        for(mm in 1:nSNP_GX){
          GX[,mm] <- rbinom(max(nX,nY),2,MAF_GX[mm])
        }
        
        for(mm in 1:nSNP_GY){
          GY[,mm] <- rbinom(max(nX,nY),2,MAF_GY[mm])
        }
        
        # refPan is from MRCD - originally from "European samples of 489 unrelated individuals in the 1000 Genomes project"
        if(independentSNPs==F | nSNP_GX==1){
          refPan <- MRCD::pruned$loci_bed[,1:nSNP_GX]} # reference panel - use number of SNPs nSNP_GX
        
        # Generate unmeasured confounder U
        if(snpEffectConfounder==F){# if FALSE, confounder U not affected by SNP GX or GY
          U<-rnorm(max(nX,nY),mean=meanU,sd=sqrt(varU))
          }
        if(snpEffectConfounder==T){ # if TRUE, SNPs GX or GY have causal effect on confounder U
          U<-rnorm(max(nX,nY),mean=(eta0 + GX%*%etaGX + GY%*%etaGY),sd=sqrt(varU))
          }
        
        # Generate X (exposure/ intermediate phenotype)
        gammaGX<-matrix(gammaGX,nrow=length(gammaGX),ncol=1)
        gammaGY<-matrix(gammaGY,nrow=length(gammaGY),ncol=1)
        
        # Unmeasured confouding True or False
        if(unmeasuredConfounding==F & contX==T){
          Xtrue <- gamma0 + GX%*%gammaGX + GY%*%gammaGY + rnorm(max(nX,nY),mean=0,sd=sqrt(varX)) 
          }
        if(unmeasuredConfounding==T & contX==T){
          Xtrue <- gamma0 + GX%*%gammaGX + GY%*%gammaGY+ gammaU*U+ rnorm(max(nX,nY),mean=0,sqrt(varX)) 
          }
        if(unmeasuredConfounding==F & contX==F){
          muTrue <- gamma0 + GX%*%gammaGX + GY%*%gammaGY 
          muTrue <- muTrue - mean(muTrue)
          Xtrue <- rbinom(max(nX,nY),1, exp(muTrue)/(1+exp(muTrue)))
          }
        if(unmeasuredConfounding==T & contX==F){
          muTrue <- gamma0 + GX%*%gammaGX + GY%*%gammaGY+ gammaU*U  
          muTrue <- muTrue - mean(muTrue)
          Xtrue <- rbinom(max(nX,nY),1, exp(muTrue)/(1+exp(muTrue)))
          }
        
        
        # Measurment error True or False
        if(measurementError==T & contX==T){X<-rnorm(max(nX,nY),(delta0+deltaX*Xtrue),sqrt(varME))}
        if(measurementError==T & contX==F){
          mu1<-delta0+deltaX*Xtrue
          mu1 <- mu1 - mean(mu1)
          X <- rbinom(max(nX,nY),1,exp(mu1)/(1+exp(mu1)))
          }
        if(measurementError==F ){X<-Xtrue}
        
        # Generate Y1 (outcome)
        if(unmeasuredConfounding==F & contY==T){ 
          Y1 <- beta0 + betaX[bX]*Xtrue +GX%*%betaGX + GY%*%betaGY+rnorm(max(nX,nY),mean=0,sd=sqrt(varY))
          }
        if(unmeasuredConfounding==T & contY==T){ 
          Y1 <- beta0 + betaX[bX]*Xtrue +GX%*%betaGX + GY%*%betaGY+betaU*U+rnorm(max(nX,nY),mean=0,sd=sqrt(varY))
          }
        
        if(unmeasuredConfounding==F & contY==F){ 
          mu1 <- beta0 + betaX[bX]*Xtrue + GX%*%betaGX + GY%*%betaGY
          mu1 <- mu1 - mean(mu1)
          Y1 <- rbinom(max(nX,nY), 1, exp(mu1)/(1+exp(mu1)))
          }
        if(unmeasuredConfounding==T & contY==F){ 
          mu1 <- beta0 + betaX[bX]*Xtrue +GX%*%betaGX + GY%*%betaGY+betaU*U
          mu1 <- mu1 - mean(mu1)
          Y1 <- rbinom(max(nX,nY), 1, exp(mu1)/(1+exp(mu1)))
          }

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
        
        
        
        ################################################################################################################################################################
        # More than 1 SNP using functions from: TwoSampleMR package, Peters et al 2014 codeANM package, MRCD package, BidirectionalMR GitHub, MendelianRandomization
        # https://github.com/MRCIEU/TwoSampleMR
        # http://web.math.ku.dk/~peters/code.html
        # https://github.com/xue-hr/MRCD
        # https://github.com/saili0103/BidirectionalMR
        # https://cran.r-project.org/web/packages/MendelianRandomization/MendelianRandomization.pdf
        ################################################################################################################################################################
        #more than 1 SNP
        
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
        
      
        # Run bidirectional MR with TwoSampleMR::mr_ivw() 
        if("tsmr_IVW" %in% colnames.methods){
          tsmr.IVW.GX <- TwoSampleMR::mr_ivw(b_exp=betaGX_X, se_exp = seGX_X, b_out=betaGX_Y, se_out = seGX_Y, parameters = default_parameters())
          tsmr.IVW.GY <- TwoSampleMR::mr_ivw(b_exp=betaGY_Y, se_exp = seGY_Y, b_out=betaGY_X, se_out = seGY_X, parameters = default_parameters())
          # TwoSampleMR mr_ivw cases
          if(tsmr.IVW.GX$pval <= sig.level & tsmr.IVW.GY$pval>sig.level){
            mat1[bX, "tsmr_IVW"] <- mat1[bX, "tsmr_IVW"] + 1
          }else if(tsmr.IVW.GX$pval > sig.level & tsmr.IVW.GY$pval<=sig.level){
            mat2[bX, "tsmr_IVW"] <- mat2[bX, "tsmr_IVW"] + 1
          }
        }
        
        # Run bidirectional MR with TwoSampleMR::mr_egger_regression()
        if("tsmr_Egger" %in% colnames.methods){
          tsmr.Egger.GX<- TwoSampleMR::mr_egger_regression(b_exp=betaGX_X, se_exp = seGX_X, b_out=betaGX_Y, se_out = seGX_Y, parameters = default_parameters())
          tsmr.Egger.GY<- TwoSampleMR::mr_egger_regression(b_exp=betaGY_Y, se_exp = seGY_Y, b_out=betaGY_X, se_out = seGY_X, parameters = default_parameters())
          # TwoSampleMR mr_egger cases
          if(tsmr.Egger.GX$pval <= sig.level & tsmr.Egger.GY$pval>sig.level){mat1[bX, "tsmr_Egger"] <- mat1[bX, "tsmr_Egger"] + 1
          }else if(tsmr.Egger.GX$pval > sig.level & tsmr.Egger.GY$pval<=sig.level){mat2[bX, "tsmr_Egger"] <- mat2[bX, "tsmr_Egger"] + 1
          } 
        }
        
        # Run bidirectional MR with TwoSampleMR::mr_egger_regression()
        if("tsmr_Egger_Boot" %in% colnames.methods){
          tsmr.Egger_boot.GX<- TwoSampleMR::mr_egger_regression_bootstrap(b_exp=betaGX_X, se_exp = seGX_X, b_out=betaGX_Y, se_out = seGX_Y, parameters = default_parameters())
          tsmr.Egger_boot.GY<- TwoSampleMR::mr_egger_regression_bootstrap(b_exp=betaGY_Y, se_exp = seGY_Y, b_out=betaGY_X, se_out = seGY_X, parameters = default_parameters())
          # TwoSampleMR mr_egger cases
          if(tsmr.Egger_boot.GX$pval <= sig.level & tsmr.Egger_boot.GY$pval>sig.level){
            mat1[bX, "tsmr_Egger_Boot"] <- mat1[bX, "tsmr_Egger_Boot"] + 1
          }else if(tsmr.Egger_boot.GX$pval > sig.level & tsmr.Egger_boot.GY$pval<=sig.level){
            mat2[bX, "tsmr_Egger_Boot"] <- mat2[bX, "tsmr_Egger_Boot"] + 1
          }       
        }
        
        # Run bidirectional MR with TwoSampleMR::mr_ivw_fe() (Inverse variance weighted regression - fixed effects)
        if("tsmr_IVW_fe" %in% colnames.methods){
          tsmr.IVW_fe.GX <- TwoSampleMR::mr_ivw_fe(b_exp=betaGX_X, se_exp = seGX_X, b_out=betaGX_Y, se_out = seGX_Y, parameters = default_parameters())
          tsmr.IVW_fe.GY <- TwoSampleMR::mr_ivw_fe(b_exp=betaGY_Y, se_exp = seGY_Y, b_out=betaGY_X, se_out = seGY_X, parameters = default_parameters())
          # TwoSampleMR mr_ivw_fe cases
          if(tsmr.IVW_fe.GX$pval <= sig.level & tsmr.IVW_fe.GY$pval>sig.level){
            mat1[bX, "tsmr_IVW_fe"] <- mat1[bX, "tsmr_IVW_fe"] + 1
          }else if(tsmr.IVW_fe.GX$pval > sig.level & tsmr.IVW_fe.GY$pval<=sig.level){
            mat2[bX, "tsmr_IVW_fe"] <- mat2[bX, "tsmr_IVW_fe"] + 1
          }
        }
        
        # Run bidirectional MR with TwoSampleMR::mr_ivw_mre() (Inverse variance weighted regression - multiplicative random effects model)
        if("tsmr_IVW_mre" %in% colnames.methods){
          tsmr.IVW_mre.GX <- TwoSampleMR::mr_ivw_mre(b_exp=betaGX_X, se_exp = seGX_X, b_out=betaGX_Y, se_out = seGX_Y, parameters = default_parameters())
          tsmr.IVW_mre.GY <- TwoSampleMR::mr_ivw_mre(b_exp=betaGY_Y, se_exp = seGY_Y, b_out=betaGY_X, se_out = seGY_X, parameters = default_parameters())
          # TwoSampleMR mr_ivw_mre cases
          if(tsmr.IVW_mre.GX$pval <= sig.level & tsmr.IVW_mre.GY$pval>sig.level){
            mat1[bX, "tsmr_IVW_mre"] <- mat1[bX, "tsmr_IVW_mre"] + 1
          }else if(tsmr.IVW_mre.GX$pval > sig.level & tsmr.IVW_mre.GY$pval<=sig.level){
            mat2[bX, "tsmr_IVW_mre"] <- mat2[bX, "tsmr_IVW_mre"] + 1
          }
        }
        
        
        # Run bidirectional MR with TwoSampleMR::mr_meta_fixed_simple() (simple standard error)
        if("tsmr_meta_fixed_simple" %in% colnames.methods){
          tsmr.meta_fixed_simple.GX <- TwoSampleMR::mr_meta_fixed_simple(b_exp=betaGX_X, se_exp = seGX_X, b_out=betaGX_Y, se_out = seGX_Y, parameters = default_parameters())
          tsmr.meta_fixed_simple.GY <- TwoSampleMR::mr_meta_fixed_simple(b_exp=betaGY_Y, se_exp = seGY_Y, b_out=betaGY_X, se_out = seGY_X, parameters = default_parameters())
          # TwoSampleMR mr_meta_fixed_simple cases
          if(tsmr.meta_fixed_simple.GX$pval <= sig.level & tsmr.meta_fixed_simple.GY$pval>sig.level){
            mat1[bX, "tsmr_meta_fixed_simple"] <- mat1[bX, "tsmr_meta_fixed_simple"] + 1
          }else if(tsmr.meta_fixed_simple.GX$pval > sig.level & tsmr.meta_fixed_simple.GY$pval<=sig.level){
            mat2[bX, "tsmr_meta_fixed_simple"] <- mat2[bX, "tsmr_meta_fixed_simple"] + 1
          }
        }
        
        
        # Run bidirectional MR with TwoSampleMR::mr_penalized_weighted_median() (penalized weighted median MR)
        if("tsmr_pwm" %in% colnames.methods){
          tsmr.penalized_weighted_median.GX <- TwoSampleMR::mr_penalised_weighted_median(b_exp=betaGX_X_1, se_exp = seGX_X_1, b_out=betaGX_Y_1, se_out = seGX_Y_1, parameters = default_parameters())
          tsmr.penalized_weighted_median.GY <- TwoSampleMR::mr_penalised_weighted_median(b_exp=betaGY_Y_1, se_exp = seGY_Y_1, b_out=betaGY_X_1, se_out = seGY_X_1, parameters = default_parameters())
          # TwoSampleMR mr_penalized_weighted_median cases
          if(tsmr.penalized_weighted_median.GX$pval <= sig.level & tsmr.penalized_weighted_median.GY$pval>sig.level){
            mat1[bX, "tsmr_pwm"] <- mat1[bX, "tsmr_pwm"] + 1
          }else if(tsmr.penalized_weighted_median.GX$pval > sig.level & tsmr.penalized_weighted_median.GY$pval<=sig.level){
            mat2[bX, "tsmr_pwm"] <- mat2[bX, "tsmr_pwm"] + 1
          }
        }
        
        # Run bidirectional MR with TwoSampleMR::mr_simple_median
        if("tsmr_simple_median" %in% colnames.methods){
          tsmr.simple_median.GX <- TwoSampleMR::mr_simple_median(b_exp=betaGX_X, se_exp = seGX_X, b_out=betaGX_Y, se_out = seGX_Y, parameters = default_parameters())
          tsmr.simple_median.GY <- TwoSampleMR::mr_simple_median(b_exp=betaGY_Y, se_exp = seGY_Y, b_out=betaGY_X, se_out = seGY_X, parameters = default_parameters())
          # TwoSampleMR mr_simple_median cases
          if(tsmr.simple_median.GX$pval <= sig.level & tsmr.simple_median.GY$pval>sig.level){
            mat1[bX, "tsmr_simple_median"] <- mat1[bX, "tsmr_simple_median"] + 1
          }else if(tsmr.simple_median.GX$pval > sig.level & tsmr.simple_median.GY$pval<=sig.level){
            mat2[bX, "tsmr_simple_median"] <- mat2[bX, "tsmr_simple_median"] + 1
          }
        }
        
        # Run bidirectional MR with TwoSampleMR::mr_uwr (Unweighted regression)
        if("tsmr_uwr" %in% colnames.methods){
          tsmr.uwr.GX <- TwoSampleMR::mr_uwr(b_exp=betaGX_X, se_exp = seGX_X, b_out=betaGX_Y, se_out = seGX_Y, parameters = default_parameters())
          tsmr.uwr.GY <- TwoSampleMR::mr_uwr(b_exp=betaGY_Y, se_exp = seGY_Y, b_out=betaGY_X, se_out = seGY_X, parameters = default_parameters())
          # TwoSampleMR mr_uwr cases
          if(tsmr.uwr.GX$pval <= sig.level & tsmr.uwr.GY$pval>sig.level){
            mat1[bX, "tsmr_uwr"] <- mat1[bX, "tsmr_uwr"] + 1
          }else if(tsmr.uwr.GX$pval > sig.level & tsmr.uwr.GY$pval<=sig.level){
            mat2[bX, "tsmr_uwr"] <- mat2[bX, "tsmr_uwr"] + 1
          }
        }
        
        # Run bidirectional MR with TwoSampleMR::mr_weighted_median
        if("tsmr_weighted_median" %in% colnames.methods){
          tsmr.weighted_median.GX <- TwoSampleMR::mr_weighted_median(b_exp=betaGX_X, se_exp = seGX_X, b_out=betaGX_Y, se_out = seGX_Y, parameters = default_parameters())
          tsmr.weighted_median.GY <- TwoSampleMR::mr_weighted_median(b_exp=betaGY_Y, se_exp = seGY_Y, b_out=betaGY_X, se_out = seGY_X, parameters = default_parameters())
          # TwoSampleMR mr_weighted_median cases
          if(tsmr.weighted_median.GX$pval <= sig.level & tsmr.weighted_median.GY$pval>sig.level){
            mat1[bX, "tsmr_weighted_median"] <- mat1[bX, "tsmr_weighted_median"] + 1
          }else if(tsmr.weighted_median.GX$pval > sig.level & tsmr.weighted_median.GY$pval<=sig.level){
            mat2[bX, "tsmr_weighted_median"] <- mat2[bX, "tsmr_weighted_median"] + 1
          }
        }
        
    if(("CDRatio" %in% colnames.methods)  | ("CDEgger" %in% colnames.methods) | ("CDgls" %in% colnames.methods)){
        # CD-Ratio (1-Sample, multiple SNPs):
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
        
        if(independentSNPs==T){
          pruned1 <- list(sig_part=tmp.df)
          cd3 <- CD_3_methods_Independent(pruned1$sig_part)
        }
        
        if(independentSNPs==F){
          pruned1 <- list(sig_part=tmp.df, loci_bed=refPan)
          cd3 <- CD_3_methods(pruned1)
        }
        
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
          stop("LowerCIxy > UpperCIxy")
          }
        if(upperCIyx.cdRatio < lowerCIyx.cdRatio){
          stop("LowerCIyx > UpperCIyx")
          }
        
        if(upperCIxy.cdEgger < lowerCIxy.cdEgger){
          stop("LowerCIxyEgger > UpperCIxyEgger")
          }
        if(upperCIyx.cdEgger < lowerCIyx.cdEgger){
          stop("LowerCIyxEgger > UpperCIyxEgger")
          }
        
        if(upperCIxy.cdGLS < lowerCIxy.cdGLS){
          stop("LowerCIxy > UpperCIxy")
          }
        if(upperCIyx.cdGLS < lowerCIyx.cdGLS){
          stop("LowerCIyx > UpperCIyx")
          }
    }
      
        
    if("CDRatio" %in% colnames.methods){    
        if((upperCIxy.cdRatio < (-1) | lowerCIxy.cdRatio>1) &
           ((lowerCIyx.cdRatio>=(-1) & upperCIyx.cdRatio<0) | (lowerCIyx.cdRatio>(0) & upperCIyx.cdRatio<=1)) ){
          mat1[bX, "CDRatio"] <- mat1[bX, "CDRatio"] + 1
          }
        if((upperCIyx.cdRatio < (-1) | lowerCIyx.cdRatio>1) &
           ((lowerCIxy.cdRatio>=(-1) & upperCIxy.cdRatio< 0) | (lowerCIxy.cdRatio>(0) & upperCIxy.cdRatio<= 1)) ){
          mat2[bX, "CDRatio"] <- mat2[bX, "CDRatio"] + 1
          }
    }
        
    if("CDEgger" %in% colnames.methods){    
        # CD-Egger decisions
        if((upperCIxy.cdEgger < (-1)  | lowerCIxy.cdEgger>1) & 
           ((lowerCIyx.cdEgger>=(-1) & upperCIyx.cdEgger<0) | (lowerCIyx.cdEgger>(0) & upperCIyx.cdEgger<=1))){
          mat1[bX, "CDEgger"] <- mat1[bX, "CDEgger"] + 1
          }
        if((upperCIyx.cdEgger < (-1) | lowerCIyx.cdEgger>1) & 
           ((lowerCIxy.cdEgger>=(-1)  & upperCIxy.cdEgger< 0) | (lowerCIxy.cdEgger>(0)  & upperCIxy.cdEgger<=1))){
          mat2[bX, "CDEgger"] <- mat2[bX, "CDEgger"] + 1
          }
    }   
     
    if("CDgls" %in% colnames.methods){   
        # CD-GLS decisions
        if((upperCIxy.cdGLS < (-1) | lowerCIxy.cdGLS>1) & 
           ((lowerCIyx.cdGLS>= (-1) & upperCIyx.cdGLS<0) | (lowerCIyx.cdGLS> (0) & upperCIyx.cdGLS<=1)) ){
          mat1[bX, "CDgls"] <- mat1[bX, "CDgls"] + 1
          }
        if((upperCIyx.cdGLS < (-1) | lowerCIyx.cdGLS>1) & 
           ((lowerCIxy.cdGLS>= (-1) & upperCIxy.cdGLS<0) | (lowerCIxy.cdGLS> (0) & upperCIxy.cdGLS<=1)) ){
          mat2[bX, "CDgls"] <- mat2[bX, "CDgls"] + 1
          }
    }
        


        
        ################################################################################
        # MR Steiger
        ################################################################################ 
        #Vector of p-values of SNP-exposure
        
        
      if(("MRS.ivw" %in% colnames.methods) | ("MRS.Egger" %in% colnames.methods) | ("MRS.wMedian" %in% colnames.methods)){
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
            ncases <- length(X==1)
            ncontrols <- length(X==0)
            prev<- (ncases)/(ncases+ncontrols)
            
            r_exp[re] <- get_r_from_lor(lor=beta, af=maf, ncase=ncases, ncontrol=ncontrols, prevalence = prev, model="logit")
            }
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
            
            r_out[ro] <- get_r_from_lor(lor=beta, af=maf, ncase=ncases, ncontrol=ncontrols, prevalence = prev, model="logit")
            }
        }
        
        #r_xxo Measurememt precision of exposure, default 1
        #r_yyo Measurement precision of outcome, default 1
        
        #A statistical test for whether the assumption that exposure causes outcome is valid
        mrs<-mr_steiger(p_exp, p_out, n_exp, n_out, r_exp, r_out, r_xxo = 1, r_yyo = 1)
      }
        
        # case 1: X->Y1 if steiger_test<alpha and weighted median pval<alpha and correct causal direction == T
        # case 2: X<-Y1 if correct causal direction == F & weighted median pval < alpha & steiger_test<alpha
        # case 3: neither if pSteiger>alpha or pMR>alpha 
        
        if(("MRS.wMedian" %in% colnames.methods) & ("tsmr_weighted_median" %in% colnames.methods)){
        if(mrs$correct_causal_direction==TRUE & mrs$steiger_test <= sig.level & tsmr.weighted_median.GX$pval <= sig.level){
          mat1[bX,"MRS.wMedian"]<-mat1[bX,"MRS.wMedian"]+1
          }
        if(mrs$correct_causal_direction==FALSE & mrs$steiger_test <= sig.level & tsmr.weighted_median.GX$pval <= sig.level){
          mat2[bX,"MRS.wMedian"]<-mat2[bX,"MRS.wMedian"]+1
          }
        }
        
        if(("MRS.wMedian" %in% colnames.methods) & !("tsmr_weighted_median" %in% colnames.methods)){
          tsmr.weighted_median.GX <- TwoSampleMR::mr_weighted_median(b_exp=betaGX_X, se_exp = seGX_X, b_out=betaGX_Y, se_out = seGX_Y, parameters = default_parameters())
          
          if(mrs$correct_causal_direction==TRUE & mrs$steiger_test <= sig.level & tsmr.weighted_median.GX$pval <= sig.level){
            mat1[bX,"MRS.wMedian"]<-mat1[bX,"MRS.wMedian"]+1}
          if(mrs$correct_causal_direction==FALSE & mrs$steiger_test <= sig.level & tsmr.weighted_median.GX$pval <= sig.level){
            mat2[bX,"MRS.wMedian"]<-mat2[bX,"MRS.wMedian"]+1}
        }
          
        
        # case 1: X->Y1 if steiger_test<alpha and MR IVW pval<alpha and correct causal direction == T
        # case 2: X<-Y1 if correct causal direction == F & MR IVW pval < alpha & steiger_test<alpha
        # case 3: neither if pSteiger>alpha or MR Egger pval >alpha 
        if(("MRS.ivw" %in% colnames.methods) & ("tsmr_IVW" %in% colnames.methods)){
          if(mrs$correct_causal_direction==TRUE & mrs$steiger_test <= sig.level & tsmr.IVW.GX$pval <= sig.level){
            mat1[bX,"MRS.ivw"]<-mat1[bX,"MRS.ivw"]+1}
          if(mrs$correct_causal_direction==FALSE & mrs$steiger_test <= sig.level & tsmr.IVW.GX$pval <= sig.level){
            mat2[bX,"MRS.ivw"]<-mat2[bX,"MRS.ivw"]+1}
        }
        
        if(("MRS.ivw" %in% colnames.methods) & !("tsmr_IVW" %in% colnames.methods)){
          tsmr.IVW.GX <- TwoSampleMR::mr_ivw(b_exp=betaGX_X, se_exp = seGX_X, b_out=betaGX_Y, se_out = seGX_Y, parameters = default_parameters())
        
          if(mrs$correct_causal_direction==TRUE & mrs$steiger_test <= sig.level & tsmr.IVW.GX$pval <= sig.level){
            mat1[bX,"MRS.ivw"]<-mat1[bX,"MRS.ivw"]+1}
          if(mrs$correct_causal_direction==FALSE & mrs$steiger_test <= sig.level & tsmr.IVW.GX$pval <= sig.level){
            mat2[bX,"MRS.ivw"]<-mat2[bX,"MRS.ivw"]+1}
        }
          
        
        # case 1: X->Y1 if steiger_test<alpha and MR Egger pval<alpha and correct causal direction == T
        # case 2: X<-Y1 if correct causal direction == F & MR Egger pval < alpha & steiger_test<alpha
        # case 3: neither if pSteiger>alpha or MR Egger pval >alpha 
        if(("MRS.Egger" %in% colnames.methods) & ("tsmr_Egger" %in% colnames.methods)){
          if(mrs$correct_causal_direction==TRUE & mrs$steiger_test <= sig.level & tsmr.Egger.GX$pval <= sig.level){
            mat1[bX,"MRS.Egger"]<-mat1[bX,"MRS.Egger"]+1}
          if(mrs$correct_causal_direction==FALSE & mrs$steiger_test <= sig.level & tsmr.Egger.GX$pval <= sig.level){
            mat2[bX,"MRS.Egger"]<-mat2[bX,"MRS.Egger"]+1}
        }
        
        if(("MRS.Egger" %in% colnames.methods) & !("tsmr_Egger" %in% colnames.methods)){
          tsmr.Egger.GX<- TwoSampleMR::mr_egger_regression(b_exp=betaGX_X, se_exp = seGX_X, b_out=betaGX_Y, se_out = seGX_Y, parameters = default_parameters())
          if(mrs$correct_causal_direction==TRUE & mrs$steiger_test <= sig.level & tsmr.Egger.GX$pval <= sig.level){
            mat1[bX,"MRS.Egger"]<-mat1[bX,"MRS.Egger"]+1}
          if(mrs$correct_causal_direction==FALSE & mrs$steiger_test <= sig.level & tsmr.Egger.GX$pval <= sig.level){
            mat2[bX,"MRS.Egger"]<-mat2[bX,"MRS.Egger"]+1}
        }
        
        
        
        ################################################################################
        # BiDirectCausal methods
        ################################################################################

        if(("BDCDcML.S.DP" %in% colnames.methods) | ("BDMRcML.S.DP" %in% colnames.methods) |
           ("BDCD.Ratio.S" %in% colnames.methods) | ("BDCD.Egger.S" %in% colnames.methods)){
          
          sig.cutoff1 = 0.05/(ncol(GX1)+ncol(GY1))
          
          # taken from https://github.com/xue-hr/BiDirectCausal -- .R files
          # Check if indices or empty
          pvalue.X1 = pnorm(-abs(c(betaGX_X, betaGY_X)/c(seGX_X, seGY_X))) * 2
          pvalue.Y1 = pnorm(-abs(c(betaGY_Y, betaGX_Y)/c(seGY_Y, seGX_Y))) * 2
          ind_X1 = which(pvalue.X1 < (sig.cutoff1))
          ind_Y1 = which(pvalue.Y1 < (sig.cutoff1))
          cor_X1 = c(betaGX_X, betaGY_X)/sqrt(c(betaGX_X, betaGY_X)^2 + (nX - 2) * c(seGX_X, seGY_X)^2)
          cor_Y1 = c(betaGY_Y, betaGX_Y)/sqrt(c(betaGY_Y, betaGX_Y)^2 + (nY - 2) * c(seGY_Y, seGX_Y)^2)
          se_corX1 = sqrt((1 - cor_X1^2)^2/nX)
          se_corY1 = sqrt((1 - cor_Y1^2)^2/nY)
          intersect.ind.X.Y1 = intersect(ind_X1, ind_Y1)
          ind_X_new1 = setdiff(ind_X1, intersect.ind.X.Y1[(abs(cor_X1)[intersect.ind.X.Y1]) <
                                                            (abs(cor_Y1)[intersect.ind.X.Y1])])
          ind_Y_new1 = setdiff(ind_Y1, intersect.ind.X.Y1[(abs(cor_X1)[intersect.ind.X.Y1]) >
                                                            (abs(cor_Y1)[intersect.ind.X.Y1])])
          # Check first! If true, return case 3!
          if(identical(ind_X_new1, integer(0)) | identical(ind_Y_new1, integer(0)) ){
            if("BDCDcML.S.DP" %in% colnames.methods){mat3[bX,"BDCDcML.S.DP"] <- mat3[bX,"BDCDcML.S.DP"]+1}
            if("BDMRcML.S.DP" %in% colnames.methods){mat3[bX,"BDMRcML.S.DP"] <- mat3[bX,"BDMRcML.S.DP"]+1}
            
            if("BDCD.Ratio.S" %in% colnames.methods){mat3[bX,"BDCD.Ratio.S"] <- mat3[bX,"BDCD.Ratio.S"]+1}
            if("BDCD.Egger.S" %in% colnames.methods){mat3[bX,"BDCD.Egger.S"] <- mat3[bX,"BDCD.Egger.S"]+1}
          }else{
            
            if("BDCDcML.S.DP" %in% colnames.methods){
              CDcMLxy=mr_cML_DP(b_exp = cor_X1[ind_X_new1],
                                b_out = cor_Y1[ind_X_new1],
                                se_exp = se_corX1[ind_X_new1],
                                se_out = se_corY1[ind_X_new1],
                                n = min(nX,nY),
                                num_pert = 100)
              
              CDcMLyx=mr_cML_DP(b_exp = cor_Y1[ind_Y_new1],
                                b_out = cor_X1[ind_Y_new1],
                                se_exp = se_corY1[ind_Y_new1],
                                se_out = se_corX1[ind_Y_new1],
                                n = min(nX,nY),
                                num_pert = 100)
              
              # CIs - CDcML.S.DP
              ci.low.XY.S.DP <- CDcMLxy$BIC_DP_theta - (qnorm(0.975)*CDcMLxy$BIC_DP_se)
              ci.upp.XY.S.DP <- CDcMLxy$BIC_DP_theta + (qnorm(0.975)*CDcMLxy$BIC_DP_se)
              
              ci.low.YX.S.DP <- CDcMLyx$BIC_DP_theta - (qnorm(0.975)*CDcMLyx$BIC_DP_se)
              ci.upp.YX.S.DP <- CDcMLyx$BIC_DP_theta + (qnorm(0.975)*CDcMLyx$BIC_DP_se)
              
              # decision rules with CIs  --- .S.DP
              if((ci.low.XY.S.DP>0 & ci.upp.XY.S.DP<1)| (ci.low.XY.S.DP>-1 & ci.upp.XY.S.DP<0)){
                if((ci.low.YX.S.DP>0 & ci.upp.YX.S.DP<1) | (ci.low.YX.S.DP>-1 & ci.upp.YX.S.DP<0)){mat3[bX,"BDCDcML.S.DP"] <- mat3[bX,"BDCDcML.S.DP"]+1
                }else{mat1[bX,"BDCDcML.S.DP"] <- mat1[bX,"BDCDcML.S.DP"]+1}
              }else if((ci.low.YX.S.DP>0 & ci.upp.YX.S.DP<1) & (ci.low.YX.S.DP>-1 & ci.upp.YX.S.DP<0)){
                if((ci.low.XY.S.DP>0 & ci.upp.XY.S.DP<1)| (ci.low.XY.S.DP>-1 & ci.upp.XY.S.DP<0)){mat3[bX,"BDCDcML.S.DP"] <- mat3[bX,"BDCDcML.S.DP"]+1
                }else{mat2[bX,"BDCDcML.S.DP"] <- mat2[bX,"BDCDcML.S.DP"]+1}
              }else{mat3[bX,"BDCDcML.S.DP"] <- mat3[bX,"BDCDcML.S.DP"]+1}
            }
            
            
            ###############################################################################
            # Apply bi-directional MRcML methods: .S.DP
            ###############################################################################
            if("BDMRcML.S.DP" %in% colnames.methods){
              MRcMLxy <- mr_cML_DP(b_exp = c(betaGX_X, betaGY_X)[ind_X_new1],
                                   b_out = c(betaGY_Y, betaGX_Y)[ind_X_new1],
                                   se_exp = c(seGX_X, seGY_X)[ind_X_new1],
                                   se_out = c(seGY_Y, seGX_Y)[ind_X_new1],
                                   n = min(nX,nY),
                                   num_pert = 100)
              MRcMLyx <- mr_cML_DP(b_exp = c(betaGY_Y, betaGX_Y)[ind_Y_new1],
                                   b_out = c(betaGX_X, betaGY_X)[ind_Y_new1],
                                   se_exp = c(seGY_Y, seGX_Y)[ind_Y_new1],
                                   se_out = c(seGX_X, seGY_X)[ind_Y_new1],
                                   n = min(nX,nY),
                                   num_pert =100)
              
              
              p.XY.BDMRcML.S.DP<- pnorm(-abs(MRcMLxy$BIC_DP_theta / MRcMLxy$BIC_DP_se))*2
              p.YX.BDMRcML.S.DP<- pnorm(-abs(MRcMLyx$BIC_DP_theta / MRcMLyx$BIC_DP_se))*2
              
              # decision rules
              if(p.XY.BDMRcML.S.DP <= 0.05 & p.YX.BDMRcML.S.DP > 0.05){mat1[bX,"BDMRcML.S.DP"] <- mat1[bX,"BDMRcML.S.DP"]+1}#return case 1
              if(p.XY.BDMRcML.S.DP > 0.05 & p.YX.BDMRcML.S.DP <= 0.05){mat2[bX,"BDMRcML.S.DP"] <- mat2[bX,"BDMRcML.S.DP"]+1}# return case 2
              if((p.XY.BDMRcML.S.DP >0.05 & p.YX.BDMRcML.S.DP > 0.05) |(p.XY.BDMRcML.S.DP <=0.05 & p.YX.BDMRcML.S.DP <= 0.05)){mat3[bX,"BDMRcML.S.DP"] <- mat3[bX,"BDMRcML.S.DP"]+1}# return case 3
            }
            
            
            ########################################
            
            # BiDirCDMethod
        if(("BDCD.Ratio.S" %in% colnames.methods) | ("BDCD.Egger.S" %in% colnames.methods)){    
            b_X = c(betaGX_X, betaGY_X)
            b_Y = c(betaGY_Y,betaGX_Y)
            se_X = c(seGX_X, seGY_X)
            se_Y = c(seGY_Y, seGX_Y)
            n_X = nX
            n_Y = nY
            sig.cutoff = 0.05/(ncol(GX1) + ncol(GY1))
            
            pvalue.X = pnorm(-abs(b_X/se_X))*2
            pvalue.Y = pnorm(-abs(b_Y/se_Y))*2
            
            ind_X = which(pvalue.X<(sig.cutoff))
            ind_Y = which(pvalue.Y<(sig.cutoff))
            
            cor_X = b_X / sqrt(b_X^2 + (n_X-2)*se_X^2)
            cor_Y = b_Y / sqrt(b_Y^2 + (n_Y-2)*se_Y^2)
            
            se_corX = sqrt((1-cor_X^2)^2/n_X)
            se_corY = sqrt((1-cor_Y^2)^2/n_Y)
            
            # With Screening
            intersect.ind.X.Y = intersect(ind_X,ind_Y)
            ind_X_new = setdiff(ind_X,
                                intersect.ind.X.Y[(abs(cor_X)[intersect.ind.X.Y])<
                                                    (abs(cor_Y)[intersect.ind.X.Y])])
            ind_Y_new = setdiff(ind_Y,
                                intersect.ind.X.Y[(abs(cor_X)[intersect.ind.X.Y])>
                                                    (abs(cor_Y)[intersect.ind.X.Y])])
            
            # CDMethods
            b_exp = b_X[ind_X_new]
            b_out = b_Y[ind_X_new]
            se_exp = se_X[ind_X_new]
            se_out = se_Y[ind_X_new]
            n_exp = n_X
            n_out = n_Y
            
            r_exp = b_exp / sqrt(b_exp^2 + (n_exp-2)*se_exp^2)
            r_out = b_out / sqrt(b_out^2 + (n_out-2)*se_out^2)
            
            ### calculate covariance matrices of r_exp and r_out
            
            p = length(r_exp)
            SNP = matrix(rnorm(p*p*100),p*100)
            SNP = scale(SNP)
            #
            rho_T1 = matrix(0, ncol = (p+1), nrow = (p+1))
            rho_T1[1:p,1:p] = diag(p)
            rho_T1[1:p,p+1] = r_exp
            rho_T1[p+1,1:p] = r_exp
            rho_T1[p+1,p+1] = 1
            #
            rho_T2 = matrix(0, ncol = (p+1), nrow = (p+1))
            rho_T2[1:p,1:p] = diag(p)
            rho_T2[1:p,p+1] = r_out
            rho_T2[p+1,1:p] = r_out
            rho_T2[p+1,p+1] = 1
            #
            #V_T1 = calculate_asymptotic_variance(SNP,rho_T1)
            #V_T2 = calculate_asymptotic_variance(SNP,rho_T2)
            V_T1 = diag((1-r_exp^2)^2)
            V_T2 = diag((1-r_out^2)^2)
            
            sig_part = data.frame(chr = 1, pos = 1, rsid = "a", A1 = "A", A2 = "G",
                                  beta_T1 = b_exp, se_T1 = se_exp, N_T1 = n_exp,
                                  beta_T2 = b_out, se_T2 = se_out, N_T2 = n_out,
                                  loci = 1)
            
            # CD_Ratio_Independent_New
            p = nrow(sig_part)
            SNP = matrix(rnorm(p * p * 100), p * 100)
            SNP = scale(SNP)
            N_T1 = sig_part[, 8]
            N_T2 = sig_part[, 11]
            T1_T = sig_part[, 6]/sig_part[, 7]
            T1_r = T1_T/sqrt(N_T1 - 2 + T1_T^2)
            T2_T = sig_part[, 9]/sig_part[, 10]
            T2_r = T2_T/sqrt(N_T2 - 2 + T2_T^2)
            rho_T1 = matrix(0, ncol = (p + 1), nrow = (p + 1))
            rho_T1[1:p, 1:p] = cor(SNP)
            rho_T1[1:p, p + 1] = T1_r
            rho_T1[p + 1, 1:p] = T1_r
            rho_T1[p + 1, p + 1] = 1
            rho_T2 = matrix(0, ncol = (p + 1), nrow = (p + 1))
            rho_T2[1:p, 1:p] = cor(SNP)
            rho_T2[1:p, p + 1] = T2_r
            rho_T2[p + 1, 1:p] = T2_r
            rho_T2[p + 1, p + 1] = 1
            if (is.null(V_T1)) {
              V_T1 = calculate_asymptotic_variance(SNP, rho_T1)
            }
            if (is.null(V_T2)) {
              V_T2 = calculate_asymptotic_variance(SNP, rho_T2)
            }
            
            if((dim(V_T1)[1]==0 & dim(V_T1)[2]==0) | (dim(V_T2)[1]==0 & dim(V_T2)[2]==0) ){
              
              if("BDCD.Ratio.S" %in% colnames.methods){mat3[bX,"BDCD.Ratio.S"] <- mat3[bX,"BDCD.Ratio.S"]+1}
              if("BDCD.Egger.S" %in% colnames.methods){mat3[bX,"BDCD.Egger.S"] <- mat3[bX,"BDCD.Egger.S"]+1}
              
            }else{    
            
              # CDMethods
              b_exp = b_Y[ind_Y_new]
              b_out = b_X[ind_Y_new]
              se_exp = se_Y[ind_Y_new]
              se_out = se_X[ind_Y_new]
              n_exp = n_Y
              n_out = n_X
              
              r_exp = b_exp / sqrt(b_exp^2 + (n_exp-2)*se_exp^2)
              r_out = b_out / sqrt(b_out^2 + (n_out-2)*se_out^2)
              
              ### calculate covariance matrices of r_exp and r_out
              
              p = length(r_exp)
              SNP = matrix(rnorm(p*p*100),p*100)
              SNP = scale(SNP)
              #
              rho_T1 = matrix(0, ncol = (p+1), nrow = (p+1))
              rho_T1[1:p,1:p] = diag(p)
              rho_T1[1:p,p+1] = r_exp
              rho_T1[p+1,1:p] = r_exp
              rho_T1[p+1,p+1] = 1
              #
              rho_T2 = matrix(0, ncol = (p+1), nrow = (p+1))
              rho_T2[1:p,1:p] = diag(p)
              rho_T2[1:p,p+1] = r_out
              rho_T2[p+1,1:p] = r_out
              rho_T2[p+1,p+1] = 1
              #
              #V_T1 = calculate_asymptotic_variance(SNP,rho_T1)
              #V_T2 = calculate_asymptotic_variance(SNP,rho_T2)
              V_T1 = diag((1-r_exp^2)^2)
              V_T2 = diag((1-r_out^2)^2)
              
              sig_part = data.frame(chr = 1, pos = 1, rsid = "a", A1 = "A", A2 = "G",
                                    beta_T1 = b_exp, se_T1 = se_exp, N_T1 = n_exp,
                                    beta_T2 = b_out, se_T2 = se_out, N_T2 = n_out,
                                    loci = 1)
              
              
              # CD_Ratio_Independent_New
              p = nrow(sig_part)
              SNP = matrix(rnorm(p * p * 100), p * 100)
              SNP = scale(SNP)
              N_T1 = sig_part[, 8]
              N_T2 = sig_part[, 11]
              T1_T = sig_part[, 6]/sig_part[, 7]
              T1_r = T1_T/sqrt(N_T1 - 2 + T1_T^2)
              T2_T = sig_part[, 9]/sig_part[, 10]
              T2_r = T2_T/sqrt(N_T2 - 2 + T2_T^2)
              rho_T1 = matrix(0, ncol = (p + 1), nrow = (p + 1))
              rho_T1[1:p, 1:p] = cor(SNP)
              rho_T1[1:p, p + 1] = T1_r
              rho_T1[p + 1, 1:p] = T1_r
              rho_T1[p + 1, p + 1] = 1
              rho_T2 = matrix(0, ncol = (p + 1), nrow = (p + 1))
              rho_T2[1:p, 1:p] = cor(SNP)
              rho_T2[1:p, p + 1] = T2_r
              rho_T2[p + 1, 1:p] = T2_r
              rho_T2[p + 1, p + 1] = 1
              if (is.null(V_T1)) {
                V_T1 = calculate_asymptotic_variance(SNP, rho_T1)
              }
              if (is.null(V_T2)) {
                V_T2 = calculate_asymptotic_variance(SNP, rho_T2)
              }
              
              if((dim(V_T1)[1]==0 & dim(V_T1)[2]==0) | (dim(V_T2)[1]==0 & dim(V_T2)[2]==0) ){
                
                if("BDCD.Ratio.S" %in% colnames.methods){mat3[bX,"BDCD.Ratio.S"] <- mat3[bX,"BDCD.Ratio.S"]+1}
                if("BDCD.Egger.S" %in% colnames.methods){mat3[bX,"BDCD.Egger.S"] <- mat3[bX,"BDCD.Egger.S"]+1}
                
              }else{
                
            
                # load functions needed
                # taken from: https://github.com/xue-hr/BiDirectCausal/blob/main/R/CDMethods.R
                CDMethods <- function(b_exp,b_out,
                                      se_exp,se_out,
                                      n_exp,n_out)
                {
                  r_exp = b_exp / sqrt(b_exp^2 + (n_exp-2)*se_exp^2)
                  r_out = b_out / sqrt(b_out^2 + (n_out-2)*se_out^2)
                  
                  ### calculate covariance matrices of r_exp and r_out
                  
                  p = length(r_exp)
                  SNP = matrix(rnorm(p*p*100),p*100)
                  SNP = scale(SNP)
                  #
                  rho_T1 = matrix(0, ncol = (p+1), nrow = (p+1))
                  rho_T1[1:p,1:p] = diag(p)
                  rho_T1[1:p,p+1] = r_exp
                  rho_T1[p+1,1:p] = r_exp
                  rho_T1[p+1,p+1] = 1
                  #
                  rho_T2 = matrix(0, ncol = (p+1), nrow = (p+1))
                  rho_T2[1:p,1:p] = diag(p)
                  rho_T2[1:p,p+1] = r_out
                  rho_T2[p+1,1:p] = r_out
                  rho_T2[p+1,p+1] = 1
                  #
                  #V_T1 = calculate_asymptotic_variance(SNP,rho_T1)
                  #V_T2 = calculate_asymptotic_variance(SNP,rho_T2)
                  V_T1 = diag((1-r_exp^2)^2)
                  V_T2 = diag((1-r_out^2)^2)
                  
                  sig_part = data.frame(chr = 1, pos = 1, rsid = "a", A1 = "A", A2 = "G",
                                        beta_T1 = b_exp, se_T1 = se_exp, N_T1 = n_exp,
                                        beta_T2 = b_out, se_T2 = se_out, N_T2 = n_out,
                                        loci = 1)
                  
                  CD_Ratio_result =
                    CD_Ratio_Independent_New(sig_part,V_T1 = V_T1,V_T2 = V_T2)
                  
                  CD_Egger_result =
                    CD_Egger_Independent(sig_part,num_iteration = 20,
                                         V_T1 = V_T1,V_T2 = V_T2)
                  
                  
                  return(list(CD_Ratio_result = CD_Ratio_result,
                              CD_Egger_result = CD_Egger_result)
                  )
                }
                
                # taken from: https://github.com/xue-hr/BiDirectCausal/blob/main/R/CD_Ratio_Independent_New.R
                CD_Ratio_Independent_New <- function (sig_part, V_T1 = NULL, V_T2 = NULL)
                {
                  p = nrow(sig_part)
                  SNP = matrix(rnorm(p * p * 100), p * 100)
                  SNP = scale(SNP)
                  N_T1 = sig_part[, 8]
                  N_T2 = sig_part[, 11]
                  T1_T = sig_part[, 6]/sig_part[, 7]
                  T1_r = T1_T/sqrt(N_T1 - 2 + T1_T^2)
                  T2_T = sig_part[, 9]/sig_part[, 10]
                  T2_r = T2_T/sqrt(N_T2 - 2 + T2_T^2)
                  rho_T1 = matrix(0, ncol = (p + 1), nrow = (p + 1))
                  rho_T1[1:p, 1:p] = cor(SNP)
                  rho_T1[1:p, p + 1] = T1_r
                  rho_T1[p + 1, 1:p] = T1_r
                  rho_T1[p + 1, p + 1] = 1
                  rho_T2 = matrix(0, ncol = (p + 1), nrow = (p + 1))
                  rho_T2[1:p, 1:p] = cor(SNP)
                  rho_T2[1:p, p + 1] = T2_r
                  rho_T2[p + 1, 1:p] = T2_r
                  rho_T2[p + 1, p + 1] = 1
                  if (is.null(V_T1)) {
                    V_T1 = calculate_asymptotic_variance(SNP, rho_T1)
                  }
                  if (is.null(V_T2)) {
                    V_T2 = calculate_asymptotic_variance(SNP, rho_T2)
                  }
                  jacobian = cbind(diag(1/T1_r), -diag(T2_r/T1_r^2))
                  combined_V = rbind(cbind(V_T2, matrix(0, ncol = p, nrow = p))/mean(N_T2),
                                     cbind(matrix(0, ncol = p, nrow = p), V_T1)/mean(N_T1))
                  V = jacobian %*% combined_V %*% t(jacobian)
                  inv_V = solve(V, tol = 0)
                  est_vec = T2_r/T1_r
                  gls_est = sum(inv_V %*% est_vec)/sum(inv_V)
                  gls_var = 1/sum(inv_V)
                  T1toT2 = c(gls_est, sqrt(gls_var))
                  Q_T1toT2 = (est_vec - gls_est) %*% inv_V %*% (est_vec - gls_est)
                  #jacobian = cbind(diag(1/T2_r), -diag(T1_r/T2_r^2))
                  #combined_V = rbind(cbind(V_T1, matrix(0, ncol = p, nrow = p))/mean(N_T1),
                  #                   cbind(matrix(0, ncol = p, nrow = p), V_T2)/mean(N_T2))
                  #V = jacobian %*% combined_V %*% t(jacobian)
                  #inv_V = solve(V, tol = 0)
                  #est_vec = T1_r/T2_r
                  #gls_est = sum(inv_V %*% est_vec)/sum(inv_V)
                  #gls_var = 1/sum(inv_V)
                  #T2toT1 = c(gls_est, sqrt(gls_var))
                  #Q_T2toT1 = (est_vec - gls_est) %*% inv_V %*% (est_vec - gls_est)
                  names(T1toT2) = c("K", "se(K)")
                  #names(T2toT1) = c("K", "se(K)")
                  return(list(T1toT2 = T1toT2, #T2toT1 = T2toT1,
                              Q_T1toT2 = Q_T1toT2
                              #Q_T2toT1 = Q_T2toT1
                  ))
                }
                
                
                # CD Egger and CD GLS - X -> Y
                CD_XtoY_S<- CDMethods(b_exp = c(betaGX_X, betaGY_X)[ind_X_new],
                                      b_out = c(betaGY_Y, betaGX_Y)[ind_X_new],
                                      se_exp = c(seGX_X, seGY_X)[ind_X_new],
                                      se_out = c(seGY_Y, seGX_Y)[ind_X_new],
                                      n_exp = nX,n_out = nY)
                # CD Egger and CD GLS - Y -> X
                CD_YtoX_S<-CDMethods(b_exp = c(betaGY_Y, betaGX_Y)[ind_Y_new],
                                     b_out = c(betaGX_X, betaGY_X)[ind_Y_new],
                                     se_exp = c(seGY_Y, seGX_Y)[ind_Y_new],
                                     se_out = c(seGX_X, seGY_X)[ind_Y_new],
                                     n_exp = nY,n_out = nX)
                
                K.ratio.XY.S = CD_XtoY_S$CD_Ratio_result$T1toT2["K"]
                seK.ratio.XY.S = CD_XtoY_S$CD_Ratio_result$T1toT2["se(K)"]
                K.egger.XY.S = CD_XtoY_S$CD_Egger_result$T1toT2["K"]
                seK.egger.XY.S = CD_XtoY_S$CD_Egger_result$T1toT2["se(K)"]
                K.ratio.YX.S = CD_YtoX_S$CD_Ratio_result$T1toT2["K"]
                seK.ratio.YX.S = CD_YtoX_S$CD_Ratio_result$T1toT2["se(K)"]
                K.egger.YX.S = CD_YtoX_S$CD_Egger_result$T1toT2["K"]
                seK.egger.YX.S = CD_YtoX_S$CD_Egger_result$T1toT2["se(K)"]
                # create 1-alpha CIs for K_XY and K_YX
                # if CI_XY completely within (-1,0) or (0,1), conclude X --> Y
                # if CI_YX completely within (-1,0) or (0,1) conclude Y--> X
                # if both true -- CASE 3
                ci.low.ratioXY.S <- K.ratio.XY.S - (qnorm(0.975)*seK.ratio.XY.S)
                ci.upp.ratioXY.S <- K.ratio.XY.S + (qnorm(0.975)*seK.ratio.XY.S)
                
                ci.low.ratioYX.S <- K.ratio.YX.S - (qnorm(0.975)*seK.ratio.YX.S)
                ci.upp.ratioYX.S <- K.ratio.YX.S + (qnorm(0.975)*seK.ratio.YX.S)
                
                ci.low.eggerXY.S <- K.egger.XY.S - (qnorm(0.975)*seK.egger.XY.S)
                ci.upp.eggerXY.S <- K.egger.XY.S + (qnorm(0.975)*seK.egger.XY.S)
                
                ci.low.eggerYX.S <- K.egger.YX.S - (qnorm(0.975)*seK.egger.YX.S)
                ci.upp.eggerYX.S <- K.egger.YX.S + (qnorm(0.975)*seK.egger.YX.S)
                
                
                #  DECISION RULES CD Ratio - .S
                if("BDCD.Ratio.S" %in% colnames.methods){
                  if((ci.low.ratioXY.S>0 & ci.upp.ratioXY.S<1)| (ci.low.ratioXY.S>-1 & ci.upp.ratioXY.S<0)){
                    if((ci.low.ratioYX.S>0 & ci.upp.ratioYX.S<1) | (ci.low.ratioYX.S>-1 & ci.upp.ratioYX.S<0)){ mat3[bX,"BDCD.Ratio.S"] <- mat3[bX,"BDCD.Ratio.S"]+1
                    }else{mat1[bX,"BDCD.Ratio.S"] <- mat1[bX,"BDCD.Ratio.S"]+1}
                  }else if((ci.low.ratioYX.S>0 & ci.upp.ratioYX.S<1) & (ci.low.ratioYX.S>-1 & ci.upp.ratioYX.S<0)){
                    if((ci.low.ratioXY.S>0 & ci.upp.ratioXY.S<1)| (ci.low.ratioXY.S>-1 & ci.upp.ratioXY.S<0)){ mat3[bX,"BDCD.Ratio.S"] <- mat3[bX,"BDCD.Ratio.S"]+1
                    }else{mat2[bX,"BDCD.Ratio.S"] <- mat2[bX,"BDCD.Ratio.S"]+1}
                  }else{ mat3[bX,"BDCD.Ratio.S"] <- mat3[bX,"BDCD.Ratio.S"]+1}
                }
                
                #  DECISION RULES CD Egger
                if("BDCD.Egger.S" %in% colnames.methods){
                  if((ci.low.eggerXY.S>0 & ci.upp.eggerXY.S<1)| (ci.low.eggerXY.S>-1 & ci.upp.eggerXY.S<0)){
                    if((ci.low.eggerYX.S>0 & ci.upp.eggerYX.S<1) | (ci.low.eggerYX.S>-1 & ci.upp.eggerYX.S<0)){mat3[bX,"BDCD.Egger.S"] <- mat3[bX,"BDCD.Egger.S"]+1
                    }else{mat1[bX,"BDCD.Egger.S"] <- mat1[bX,"BDCD.Egger.S"]+1}
                  }else if((ci.low.eggerYX.S>0 & ci.upp.eggerYX.S<1) & (ci.low.eggerYX.S>-1 & ci.upp.eggerYX.S<0)){
                    if((ci.low.eggerXY.S>0 & ci.upp.eggerXY.S<1)| (ci.low.eggerXY.S>-1 & ci.upp.eggerXY.S<0)){mat3[bX,"BDCD.Egger.S"] <- mat3[bX,"BDCD.Egger.S"]+1
                    }else{mat2[bX,"BDCD.Egger.S"] <- mat2[bX,"BDCD.Egger.S"]+1}
                  }else{mat3[bX,"BDCD.Egger.S"] <- mat3[bX,"BDCD.Egger.S"]+1}
                }
                
              } # end third else loop    
            } # end second else loop
        } # end first if loop
            
              
              
              
              
              
              
              
              
              
          }      
          

        }         
        ###############################################################
        
        ################################################################################
        # end loops
        ################################################################################    
      }#beta loop
    }#sim loop
    
    ################################################################################
    # save results
    ################################################################################    
    mat_total1 <- cbind(mat1[,1],(mat1[,-1]/nSim))
    mat_total2 <- cbind(mat2[,1],(mat2[,-1]/nSim))
    mat_total3 <- cbind(mat3[,1],(1-(mat_total1[,-1]+mat_total2[,-1])))
    
    if(ncol(mat1)>2){
    colnames(mat_total1)[1] <- c("betaX")
    colnames(mat_total2)[1] <- c("betaX")
    colnames(mat_total3)[1] <- c("betaX")
    }
    if(ncol(mat1) == 2){
      colnames(mat_total1) <- c("betaX", colnames.methods[2])
      colnames(mat_total2) <- c("betaX", colnames.methods[2])
      colnames(mat_total3) <- c("betaX", colnames.methods[2])
    }
    
    pl=F
    if(any(gammaGY!=0) | any(betaGX!=0)){pl=T}
    
    write.table(mat_total1,file=paste0("matCase1_sims", nSim,"seed", SEED,"nX",nX,"nY",nY,"nSNPX", length(MAF_GX),"nSNPY", length(MAF_GY), "contX", contX, "contY",contY,"u",unmeasuredConfounding, "me",measurementError,"p",pl, "l1", long1, "l2",long2,"dx", deltaX, "bU", betaU, "gU",gammaU,"bX",betaGX[1]*100, "gY", gammaGY[1]*100,"eGX", etaGX[1]*100, "eGY", etaGY[1]*100,".txt"), quote=F, row.names = F)
    write.table(mat_total2,file=paste0("matCase2_sims", nSim,"seed", SEED,"nX",nX,"nY",nY,"nSNPX", length(MAF_GX),"nSNPY", length(MAF_GY),"contX", contX, "contY",contY,"u",unmeasuredConfounding, "me",measurementError,"p",pl, "l1", long1, "l2",long2,"dx", deltaX,  "bU", betaU, "gU",gammaU,"bX",betaGX[1]*100, "gY", gammaGY[1]*100,"eGX", etaGX[1]*100, "eGY", etaGY[1]*100 ,".txt"), quote=F, row.names = F)
    write.table(mat_total3,file=paste0("matCase3_sims", nSim,"seed", SEED,"nX",nX, "nY",nY,"nSNPX", length(MAF_GX),"nSNPY", length(MAF_GY),"contX", contX, "contY",contY,"u",unmeasuredConfounding, "me",measurementError,"p",pl, "l1", long1, "l2",long2,"dx", deltaX,  "bU", betaU, "gU",gammaU,"bX",betaGX[1]*100, "gY", gammaGY[1]*100,"eGX", etaGX[1]*100, "eGY", etaGY[1]*100 ,".txt"), quote=F, row.names = F)
    
    return(list("mat_total1"=mat_total1, "mat_total2"=mat_total2, "mat_total3"=mat_total3))
  }
