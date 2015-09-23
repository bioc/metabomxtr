#This function performs per-metabolite batch normalization using a mixture model with batch-specific thresholds and run order correction if desired

mixnorm <- function(ynames,batch="Batch",mxtrModel=NULL,batchTvals=NULL,removeCorrection=NULL,nNA=5,minProp=0.2,method="BFGS",cData,data){
  
  #If not specified, indicate default mxtrModel 
  if (is.null(mxtrModel)) mxtrModel <- as.formula(paste("~",batch,"|",batch))
  mxtrModel2 <- as.Formula(mxtrModel)

  #confirm that all variables included in the mixture model are in both data and cData, in compatible forms 
  model.vars<-all.vars(mxtrModel2)
  lapply(c(model.vars,removeCorrection),function(x){
		#make sure all variables are in both data sets
		if (!x %in% colnames(cData)){
			stop(paste(x,"is not in cData."))
		} 
		if (!x %in% colnames(data)){
			stop(paste(x,"is not in data."))
		}
		#make sure the variable types are compatible in both data sets 
		cdata.type<-class(cData[,x])
		data.type<-class(data[,x])
		if (cdata.type!=data.type){
			stop(paste(x,"is of type",cdata.type,"in cData but type",data.type,"in data."))
		}		
		#make sure levels of factor/character variables are compatible
		if (cdata.type=="character" & data.type=="character"){
			cdata.factor<-as.factor(cData[,x])
			data.factor<-as.factor(data[,x])
			if (!all(levels(cdata.factor)==levels(data.factor))){
				stop(paste(x,"does not have the same levels in cData and data."))
			}
		}
		if (cdata.type=="factor" & data.type=="factor"){
			if (!all(levels(cData[,x])==levels(data[,x]))){
				stop(paste(x,"does not have the same levels in cData and data."))
			}
		}
		#make sure batch is coded as a factor 
		if (x==batch){
			if (cdata.type!="factor" | data.type!="factor"){
				stop(paste("Batch variable",x,"must be coded as a factor."))
			}
		}
		#make sure there are no missing values for the covariates in the model 
		if (!x %in% ynames){
			if (any(is.na(cData[,x]))){
				stop(paste(x,"is missing values in cData. Missing values for model covariates are not permitted."))
			}
			if (any(is.na(data[,x]))){
				stop(paste(x,"is missing values in data. Missing values for model covariates are not permitted."))
			}
		}
	})
  
  #Define data frames containing response variables in control and observed data
  obsYc <- yvals(cData,ynames)
  obsY <- yvals(data,ynames)
  
  #Define the design matrices for the continuous portion of the model 
  zVarsc<-zdesign(cData,mxtrModel2)
  zVars<-zdesign(data,mxtrModel2)
  stopifnot(dim(zVarsc)[2]==dim(zVars)[2])
   
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

  #By default, the effects of all variables included in the mixture model 
  #will be subtracted from the non-normalized data. However, if desired, 
  #variables specified via removeCorrection can be included in the mixture 
  #model as covariates, but their estimated effects will not be subtracted 
  #from the raw data. 
  if(!is.null(removeCorrection)){
	#get names of variables to remove 
	remove.these<-unlist(lapply(removeCorrection,function(x){
							#make sure removeCorrection variables are in the mixture model
							if (! x %in% model.vars){
								stop(paste(x,"is specified in removeCorrection but is not included in the mixture model."))
							}
							#then identify the names of the variables in the mixture model output normParams
							varname<-paste(x,levels(as.factor(data[,x]))[-1],sep="")
							return(varname)
					}))

	#then set effects to zero 
	zVarsc[,remove.these] <- 0
	zVars[,remove.these] <- 0
  }
   
  #Determine correction values   
  corrValc <- zVarsc %*% t(normParamsZ)
  corrVal <- zVars %*% t(normParamsZ)

  #Subtract correction values from observed data
  ctlNorm <- obsYc[,colnames(corrValc)] - corrValc 
  obsNorm <- obsY[,colnames(corrVal)] - corrVal 

  ans <- list(normParamsZ=normParamsZ,ctlNorm=ctlNorm,obsNorm=obsNorm)
  ans
}


