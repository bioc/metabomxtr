#This function performs per-metabolite batch normalization using a mixture model with batch-specific thresholds and run order correction if desired

mixnorm <- function(ynames, batch="Batch", mxtrModel=NULL, cData, data, batchTvals=NULL, removeCorrection=NULL, nNA=5,minProp=0.2,method="BFGS"){

		#If not specified, indicate default mxtrModel 
		if (is.null(mxtrModel)) mxtrModel <- as.formula(paste("~",batch,"|",batch))
		mxtrModelF <- as.Formula(mxtrModel)
		
		#in the input data are expression sets, convert to data frames 
		if (class(cData)=="ExpressionSet"){
		
			cData<-data.frame(cbind(t(exprs(cData)), pData(cData)))
		
		}
		
		if (class(data)=="ExpressionSet"){
		
			data<-data.frame(cbind(t(exprs(data)), pData(data)))
		
		}

		#confirm that all variables included in the mixture model are in both data and cData, in compatible forms 
		model.vars<-all.vars(mxtrModelF)
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
  
	    #create data frame containing response variable in control and observed data 
	    obsYc<- yvals(cData,ynames)
	    obsY <- yvals(data,ynames)
	   
	    #Define the design matrices for the continuous portion of the model 
	    zVarsc<-zdesign(cData,mxtrModelF)
	    zVars<-zdesign(data,mxtrModelF)
	    stopifnot(dim(zVarsc)[2]==dim(zVars)[2])
	  
	    #store order of columns in the full, control data set 
	    full.col.order<-colnames(zVarsc)
	    if ("(Intercept)" %in% full.col.order){
	  
			final.col.order<-c("zInt", paste0("z_", full.col.order[-1]))
		
	    } else {
	  
			final.col.order<-paste0("z_", full.col.order[-1])
	  
	    }
	   
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
	    names(Tvals)<-paste("batch", batch, cData[ , batch])

	    #Estimate parameters from control data to use for normalization
	    normParams <- mxtrmod(ynames=ynames,mxtrModel=mxtrModel,Tvals=Tvals,
							nNA=nNA,minProp=minProp,method=method,
							data=cData)
	    zIntI <- which(colnames(normParams)=="zInt")
	    if(length(zIntI)==0) zIntI <- NULL
	    normParamsZ <- normParams[,c(zIntI,grep("z_",colnames(normParams)))]
	    rownames(normParamsZ) <- normParams[,".id"]
	  
	    #Also, store model convergence information
	    model.conv<-normParams[ ,c(".id", "conv")]
	  
	    #add warning about missing data to model convergence data frame, if necessary 
	    if ("missing_levels_warning" %in% colnames(normParams)){
	  
			model.conv$missing_levels_warning<-normParams$missing_levels_warning
			model.conv$Warning<-ifelse(is.na(model.conv$missing_levels_warning), NA, "Effects of missing levels could not be estimated and were not included in normalization.")
	  
	    }

	    #Calculate correction values for each sample for all observed metabolites 
	    #Don't include correction for (Intercept)
	    norm.param.colnames<-colnames(normParamsZ)
	    normParamsZCopy<-normParamsZ
	    if ("zInt" %in% norm.param.colnames){ 
		
		  normParamsZCopy$zInt<-0
		
	    } 
	  
	    #also, if missing levels of categorical variables exist, 
	    #add parameter estimates equal to 0 to the normParamsZ data frame
		#and do not account for that factor level in normalization
	    missing.params<-final.col.order[! final.col.order %in% colnames(normParamsZ)]
	    na.vals<-any(apply(normParamsZ, 1, function(x) any(is.na(x))))
	    if (length(missing.params)>0){
	  
		  normParamsZCopy[ , missing.params]<-0
		  normParamsZCopy<-normParamsZCopy[ , colnames(normParamsZCopy)[match(final.col.order, colnames(normParamsZCopy))]]
	  
	    } else if (na.vals){
	  
		  for (i in 1:nrow(normParamsZCopy)){
		
			normParamsZCopy[i, ][is.na(normParamsZCopy[i, ])]<-0
		
		  }
	  
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
								param.name<-paste0("z_", varname)
								return(param.name)
						}))

		  #then set effects to zero 
		  normParamsZCopy[ , remove.these] <- 0
		
	    }
	   
	    #make sure normalization parameters are in the same order as the columns in the 
	    #obsY and obsYc datasets
	    metab.order<-colnames(obsY)
	    norm.params.order<-match(metab.order, row.names(normParamsZ))
	    normParamsZ<-normParamsZ[ norm.params.order, ]
	    normParamsZCopy<-normParamsZCopy[ norm.params.order, ]
	    stopifnot(all(row.names(normParamsZCopy)==metab.order))
	   
	    #Determine correction values  		  
	    corrValc <- zVarsc %*% t(normParamsZCopy)
	    corrVal <- zVars %*% t(normParamsZCopy)

	    #Subtract correction values from observed data
	    ctlNorm <- obsYc[,colnames(corrValc)] - corrValc 
	    obsNorm <- obsY[,colnames(corrVal)] - corrVal 
		  
	    #return list of results 
	    ans <- list(normParamsZ=normParamsZ, ctlNorm=ctlNorm, obsNorm=obsNorm, conv=model.conv)
	   return(ans)
				  
}


