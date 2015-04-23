#This function performs per-metabolite batch normalization using a mixture model with batch-specific thresholds and run order correction if desired

mxtrnorm <- function(ynames,batch="Batch",mxtrModel=NULL,batchTvals=NULL,correctSampleType=TRUE,sampleType=NULL,nNA=5,minProp=0.2,method="BFGS",cData,data){
  
  #Make sure batches are coded as factors with the same levels
  stopifnot(all(levels(cData[,batch])==levels(data[,batch])))

 
  
  #If not specified, indicate default mxtrModel 
  if (is.null(mxtrModel)) mxtrModel <- as.formula(paste("~",batch,"|",batch))
  mxtrModel2 <- as.Formula(mxtrModel)

  #By default, when there are control samples of different types,
  #normalized values are returned with a location shift correction for
  #control sample type.
  #If correctSampleType=FALSE, then normalized values will be returned without
  #this location correct.  In this case, sampleType must be specified and must be
  #included in at least one portion of the mxtrModel formula
  #Also make sure sampleType is coded as a factor with the same levels
  
  if(!correctSampleType){
    stopifnot(!is.null(sampleType))
    stopifnot(sampleType %in% all.vars(mxtrModel2))
    stopifnot(all(levels(cData[,sampleType])==levels(data[,sampleType])))

    #get variable names that will be associated with controlSampleType after correction
    sampleTypeVarNames <- paste(sampleType,levels(data[,sampleType])[-1],sep="")
  }
  
  #Define data frames containing response variables in control and observed data
  obsYc <- yvals(cData,ynames)
  obsY <- yvals(data,ynames)

  #Define the design matrices for the discrete portion of the model
  xVarsc<-xdesign(cData,mxtrModel2)
  xVars<-xdesign(data,mxtrModel2)
  stopifnot(dim(xVarsc)[2]==dim(xVars)[2])
  
  #Define the design matrices for the continuous portion of the model 
  zVarsc<-zdesign(cData,mxtrModel2)
  zVars<-zdesign(data,mxtrModel2)
  stopifnot(dim(zVarsc)[2]==dim(zVars)[2])
  

  #Determine if any of the observations have missing covariate values
  #Stop normalization if any values are missing
  stopifnot(all(apply(xVarsc,1,function(x){all(!is.na(x))})))
  stopifnot(all(apply(xVars,1,function(x){all(!is.na(x))})))
  stopifnot(all(apply(zVarsc,1,function(x){all(!is.na(x))})))
  stopifnot(all(apply(zVars,1,function(x){all(!is.na(x))})))
  


  
  #Find batch specific minima if batchTvals is not defined
  if(is.null(batchTvals)){
    cMin <- apply(obsYc,1,FUN=min,na.rm=TRUE)
    oMin <- apply(obsY,1,FUN=min,na.rm=TRUE)
    cMinB <- tapply(cMin,INDEX=cData[,batch],FUN=min,na.rm=TRUE)
    oMinB <- tapply(oMin,INDEX=data[,batch],FUN=min,na.rm=TRUE)
    batchTvals <- pmin(cMinB,oMinB)
  }

  #Create Tvals to assign batch-specific threshold to individuals
  Tvals <- batchTvals[cData[,batch]]

  #Estimate parameters from control data to use for normalization
  normParams <- mxtrmod(ynames=ynames,mxtrModel=mxtrModel,Tvals=Tvals,
                        nNA=nNA,minProp=minProp,method=method,
                        data=cData)
  zIntI <- which(colnames(normParams)=="zInt")
  if(length(zIntI)==0) zIntI <- NULL
  normParamsZ <- normParams[,c(zIntI,grep("z_",colnames(normParams)))]
  rownames(normParamsZ) <- normParams[,".id"]

  #Calculate correction values for each sample for all observed metabolites 
  #Don't include correction for (Intercept)
  if ("(Intercept)" %in% colnames(zVarsc)){
    zVarsc[,"(Intercept)"] <- rep(0,dim(zVarsc)[1])
    zVars[,"(Intercept)"] <- rep(0,dim(zVars)[1])
  }
  if(!correctSampleType){
    zVarsc[,sampleTypeVarNames] <- 0
    zVars[,sampleTypeVarNames] <- 0
  }
     
  corrValc <- zVarsc %*% t(normParamsZ)
  corrVal <- zVars %*% t(normParamsZ)

  #Subtract correction values from observed data
  ctlNorm <- obsYc[,colnames(corrValc)] - corrValc 
  obsNorm <- obsY[,colnames(corrVal)] - corrVal 

  ans <- list(normParamsZ=normParamsZ,ctlNorm=ctlNorm,obsNorm=obsNorm)
  ans
}


