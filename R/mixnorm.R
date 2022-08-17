#This function performs per-metabolite batch normalization using a mixture model with batch-specific thresholds and run order correction if desired

mixnorm <- function(ynames, batch="Batch", mxtrModel=NULL, cData, data, batchTvals=NULL, removeCorrection=NULL, nNA=5,minProp=0.2,method="BFGS", qc.sd.outliers=2){

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
		
		#confirm that all variables have enough observations to run models, in both the control and observational data  
		ynamesC<-unlist(lapply(ynames, function(metab){
		
			metab.vec<-cData[ , metab]
			if( sum(!is.na(metab.vec)) > minProp*length(metab.vec)){
			
				metab
			
			} else {
			
				warning(paste("There is too much missing control data to proceed with normalization for", metab))
				return(NULL)
			
			} 		
		
		}))
		
		ynamesO<-unlist(lapply(ynames, function(metab){
		
			metab.vec<-data[ , metab]
			if( sum(!is.na(metab.vec)) > minProp*length(metab.vec)){
			
				metab
			
			} else {
			
				warning(paste("There is too much missing experimental data to proceed with normalization for", metab))
				return(NULL)
			} 		
		
		}))
		
		if (is.null(ynamesC)){
		
			stop("All of the input metabolites are missing too much control data to proceed with normalization.")
		
		} else if (is.null(ynamesO)){
		
			stop("All of the input metabolites are missing too much experimental data to proceed with normalization.")
		
		}
		
		ynames<-intersect(ynamesC, ynamesO)
		
		#confirm that all variables included in the mixture model are in both data and cData, in compatible forms 
		model.vars<-all.vars(mxtrModelF)
		model.var.classes<-sapply(model.vars, function(x){ class(data[ ,x]) })
		cat.vars<-names(model.var.classes)[model.var.classes %in% c("factor", "character")]
		
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
		  cMin <- apply(obsYc,1, function(x){
				
					if (all(is.na(x))){
					
						NA
					
					} else {
					
						min(x, na.rm=TRUE)
					
					}
				})
		  oMin <- apply(obsY,1, function(x){
				
					if (all(is.na(x))){
					
						NA
					
					} else {
					
						min(x, na.rm=TRUE)
					
					}
				})
		  cMinB <- tapply(cMin,INDEX=cData[,batch], function(x){
					
					if (all(is.na(x))){
					
						NA
					
					} else {
					
						min(x, na.rm=TRUE)
					
					}		
		  
				})
		  
		  oMinB <- tapply(oMin,INDEX=data[,batch], function(x){
					
					if (all(is.na(x))){
					
						NA
					
					} else {
					
						min(x, na.rm=TRUE)
					
					}		
		  
				})
				
		  batchTvals <- pmin(cMinB,oMinB)
	    }

	    #Create Tvals to assign batch-specific threshold to individuals
	    Tvals <- batchTvals[cData[,batch]]
	    names(Tvals)<-paste("batch", batch, cData[ , batch])

	    #Estimate parameters from control data to use for normalization
	    normParams <- mxtrmod(ynames=ynames,mxtrModel=mxtrModel,Tvals=Tvals,
							nNA=nNA,minProp=minProp,method=method,
							data=cData, remove.outlier.sd=qc.sd.outliers)
														
		#determine the final metabolites for which normalization was able to run 
		norm.metabs<-as.character(normParams[ , ".id"])
		
		#set outlying qc data values to Inf
		cData[ , norm.metabs]<-apply(cData[ , norm.metabs, drop=FALSE], 2, function(qc.metab.vals){
		
			#determine whether the values fall outside the specified standard deviation threshold 
			#and set to infinite if yes 
			mean.metab.abundance<-mean(qc.metab.vals, na.rm=TRUE)
			sd.metab.abundance<-sd(qc.metab.vals, na.rm=TRUE)
			outlier<-qc.metab.vals > mean.metab.abundance + qc.sd.outliers*sd.metab.abundance | qc.metab.vals < mean.metab.abundance - qc.sd.outliers*sd.metab.abundance
			outlying.values<-!is.na(qc.metab.vals) & outlier
			qc.metab.vals[outlying.values]<-Inf
			qc.metab.vals
		
		})
		
		#grab parameter estimates for the continuous portion of the model 		
	    zIntI <- which(colnames(normParams)=="zInt")
	    if(length(zIntI)==0) zIntI <- NULL
	    normParamsZ <- normParams[,c(zIntI,grep("z_",colnames(normParams)))]
	    rownames(normParamsZ) <- normParams[,".id"]
	  
	    #Also, store model convergence information
	    model.conv<-normParams[ ,c(".id", "conv")]
	  
	    #add warning about missing data to model convergence data frame, if necessary 
		#loop through metabolites and set values equal to infinite where model predictors are missing too much data 
	    if ("predictors_missing_levels" %in% colnames(normParams)){
	  
			model.conv$predictors_missing_levels<-normParams$predictors_missing_levels
			missing.level.metabs<-as.character(normParams[!is.na(normParams$predictors_missing_levels), ".id"])
			
			#loop over metabolites with missing levels 
			for (metab in missing.level.metabs){
			
				missing.level.note<-normParams[normParams$".id"==metab, "predictors_missing_levels"]
				missing.level.note.list<-strsplit(missing.level.note, "; ", fixed=TRUE)
				
				#then loop over predictors missing levels 
				for (missing.level.note.element in missing.level.note.list){
				
					note.split<-strsplit(missing.level.note.element, ":", fixed=TRUE)[[1]]
					variable<-note.split[1]
					missing.levels<-strsplit(note.split[2], ", ", fixed=TRUE)[[1]]
					
					#then for each combination of metabolite and predictor level, set metabolite abundance to infinite 
					for (missing.level in missing.levels){
					
							cData[cData[ , variable]==missing.level & !is.na(cData[ , metab]), metab]<-Inf
							data[data[ , variable]==missing.level & !is.na(data[ , metab]), metab]<-Inf	
					
					}
				}
			}
	    }
		
		#check for instances where an particular combination of categorical predictors is entirely missing metabolite 
		#QC data (e.g., mother samples from batch 1), and remove data from the experimental data set if so
		if (length(cat.vars) > 0){
		
			cat.var.list<-lapply(cat.vars, function(var.name) cData[ , var.name])
			names(cat.var.list)<-cat.vars
			cat.var.combns<-expand.grid(cat.var.list)
			unique.cat.combns<-unique(cat.var.combns)
		
			for (metab in norm.metabs){
		
				for (row.num in 1:nrow(unique.cat.combns)){
				
					combn.test<-merge(cData[ , c(metab, cat.vars)], unique.cat.combns[row.num, , drop =F], by=cat.vars)
					combn.test.target<-combn.test[!is.na(combn.test[ , metab]) & !is.infinite(combn.test[ , metab]), ]
					if (nrow(combn.test.target)==0){
					
						data.copy<-data
						data.copy$row.number<-1:nrow(data.copy)
						remove.row.df<-merge(data.copy[ ,c(metab, cat.vars, "row.number")], unique.cat.combns[row.num, , drop = F], by=cat.vars)
						remove.row.df<-remove.row.df[!is.na(remove.row.df[ , metab]), ]
						problem.rows<-remove.row.df$row.number
						
						data[problem.rows, metab]<-Inf
					
					}
				
				}
			}	
		}		
		
		#now check for cases where a variable was excluded from a model 
		#and set metabolite values to 0 in those instances 
		if ("excluded_predictors" %in% colnames(normParams)){
		
			excluded.predictor.metabs<-normParams[!is.na(normParams$excluded_predictors), ".id"]
			for (metab in excluded.predictor.metabs){
			
				cData[!is.na(cData[ , metab]), metab]<-Inf
				data[!is.na(data[ , metab]), metab]<-Inf	
			
			}
		
			model.conv$excluded_predictors<-normParams$excluded_predictors
				
		}
		
		#recreate data frame containing response variable in control and observed data, using updated data  
	    updated.obsYc<- yvals(cData, norm.metabs)
	    updated.obsY <- yvals(data, norm.metabs)

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
		
		#first, deal with categorical variables with too much missing data for 
		#all input metabolites 
	    if (length(missing.params)>0){
	  
		  normParamsZCopy[ , missing.params]<-0
		  
		  #make sure columns are in the correct order 
		  normParamsZCopy<-normParamsZCopy[ , colnames(normParamsZCopy)[match(final.col.order, colnames(normParamsZCopy))]]

		#then identify and correct cases where some but not all metabolites 
		#are missing too much data 
	    } else if (na.vals){
	  
		  for (i in 1:nrow(normParamsZCopy)){
		
			if (any(is.na(normParamsZCopy[i, ]))){
			
				#convert NA parameter to 0 
				normParamsZCopy[i, is.na(normParamsZCopy[i, ])]<-0
				
				
				}
			}
			
			#make sure columns are in the correct order 
			normParamsZCopy<-normParamsZCopy[ , colnames(normParamsZCopy)[match(final.col.order, colnames(normParamsZCopy))]]
			
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
	    metab.order<-colnames(updated.obsY)
	    norm.params.order<-match(metab.order, row.names(normParamsZ))
	    normParamsZ<-normParamsZ[ norm.params.order, ]
	    normParamsZCopy<-normParamsZCopy[ norm.params.order, ]
	    stopifnot(all(row.names(normParamsZCopy)==metab.order))
		stopifnot(all(colnames(normParamsZCopy)==final.col.order))
	   
	    #Determine correction values  		  
	    corrValc <- zVarsc %*% t(normParamsZCopy)
	    corrVal <- zVars %*% t(normParamsZCopy)
		
	    #Subtract correction values from observed data
	    ctlNorm <- updated.obsYc[ ,colnames(corrValc)] - corrValc 
	    obsNorm <- updated.obsY[ ,colnames(corrVal)] - corrVal
	
	    #return list of results 
	    ans <- list(normParamsZ=normParamsZ, ctlNorm=ctlNorm, obsNorm=obsNorm, conv=model.conv)
	   return(ans)
				  
}


