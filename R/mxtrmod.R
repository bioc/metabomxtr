mxtrmod <-function(ynames,mxtrModel,Tvals=NULL,nNA=5,minProp=0.2,method="BFGS",data,fullModel=NULL, remove.outlier.sd=NULL){

  ###############################################################################################
  # Determines optimized estimates for all parameters in a mixture model.  
  #
  # Args: 
  #  ynames:   A character vector of the mixture model outcome names, e.g. metabolites.
  #            If the input data object is a matrix or data frame, these should be column
  #            names. If the input data object is an expression set, these should be row 
  #            names. Response variables should have a normal or lognormal distribution. 
  #            If lognormal,log transformed variables should be input. Missing values should 
  #            be denoted by NA.   
  # mxtrModel: A formula of the form ~x1+x2...|z1+z2..., where x's are the names of 
  #            covariates included in the discrete portion of the model and z's are names
  #            of covariates included in the continuous portion. For intercept only models, 
  #            enter 1 instead of covariate names on the appropriate side of the |.  
  #  Tvals:    A vector of thresholds below which continuous variables are not 
  #            observable. By default, this parameter will be set to the minimum
  #            of the response variable.                
  #  nNA:      The minimum number of unobserved response values for the 
  #            discrete portion of the model to be included in the likelihood. 
  #            Models for variables with fewer than nNA missing values will include only 
  #            the continuous portion of the likelihood. The default value is 5.
  #  minProp:  The minimum proportion of non-missing values in the response variable
  #            necessary to run the model. The default value is 0.2. Models
  #            will not be run if more than 80% of response variable values are missing. 
  #  method:   The method used to optimize the parameter estimates of the mixture model. 
  #            "BFGS" is the default method. Other options are documented in the 
  #            manual for the function 'optimx' in package optimx.    
  #  data:     The input data object. Matrices, data frames, and expression sets are 
  #            all acceptable classes. If a data frame or matrix, rows are subjects 
  #            and columns are metabolites or outcomes.   
  #  fullModel: A formula of the form ~x1+x2...|z1+z2..., where x's are the names of 
  #            covariates included in the discrete portion of the full model and z's are names
  #            of covariates included in the continuous portion. Input if the mxtrModel parameter
  #            represents a reduced model. 
  #  remove.outlier.sd: The maximum number of standard deviations from the mean that should be considered 
  #						non-outlying metabolite abundance values. Metabolite abundances greater than 
  #						this will be removed from modeling. This defaults to NULL, meaning
  #						all values will be included. 
  # Returns:
  #   A data frame with optimized estimates for all parameters in the mixture model, the negative 
  #   log likelihood, the optimization method used, whether the algorithm converged, and the 
  #  number of observations used.  
  #################################################################################################
  
  #first, determine the categorical variables in the model 
  #this should be based on the full model 
  if (!is.null(fullModel)){
  
	mxtrModelTarget<-fullModel
  
  } else{
  
	mxtrModelTarget<-mxtrModel
  
  }
  
  #in the input data are expression sets, convert to data frames 
  if (class(data)=="ExpressionSet"){

	  data<-data.frame(cbind(t(exprs(data)), pData(data)))

  }

  #identify categorical variables in model 
  model.vars<-all.vars(mxtrModelTarget)
  model.var.classes<-sapply(model.vars, function(x){ class(data[ ,x]) })
  cat.vars<-names(model.var.classes)[model.var.classes %in% c("factor", "character")]
  
  #convert character variables to factors 
  character.vars<-names(model.var.classes)[model.var.classes == "character"]
  for (character.var in character.vars){
  
	data[ , character.var]<-as.factor(data[ , character.var])
  
  }
  
  #determine the proper order of the categorical variables 
  cat.var.levels<-lapply(cat.vars, function(cat.var) levels(data[ , cat.var]))
  names(cat.var.levels)<-cat.vars
  
  #loop over input metabolites 
  results.list<-bplapply(ynames, function(var.name){
	
	#subset the data to just the metabolite and predictor variables 
	data<-data[ , c(var.name, model.vars)]

	#if requested, remove outlying values
	metab.vec<-data[ , var.name]
	if (!is.null(remove.outlier.sd)){
		
		metab.vec.no.inf<-metab.vec[!is.infinite(metab.vec)]
		mean.metab.abundance<-mean(metab.vec.no.inf, na.rm=TRUE)
		sd.metab.abundance<-sd(metab.vec.no.inf, na.rm=TRUE)
		remove.these.rows<- !is.na(metab.vec) & !is.infinite(metab.vec) & 
								(metab.vec > mean.metab.abundance + remove.outlier.sd*sd.metab.abundance |  metab.vec < mean.metab.abundance - remove.outlier.sd*sd.metab.abundance)
		metab.vec<-metab.vec[!remove.these.rows]
		data<-data[!remove.these.rows, ]
		
	}	
	
	#determine whether any of the levels of categorical variables are entirely missing
	missing.levels.check<-lapply(var.name, anyMissingLevels, cat.vars=cat.vars, dataset=data)
	names(missing.levels.check)<-var.name
	any.missing.levels<-sapply(missing.levels.check, any)
	missing.level.vars<-names(any.missing.levels)[any.missing.levels]
	no.missing.level.vars<-names(any.missing.levels)[!any.missing.levels]
   
	#if there are any variables entirely missing for at least one level of the categorical variables, 
	#remove those levels from the data and run a mixture model
	if (var.name %in% missing.level.vars){
		
		#confirm that we actually have enough data to model 
		if(sum(!is.na(metab.vec) & !is.infinite(metab.vec)) > minProp*length(metab.vec)){
		
			#identify actual levels that are missing 
			missing.level.list<-idMissingLevels(var.name, missing.levels.check, data)
			
			#make a copy for use later in determining if batch Tvals need to be removed 
			missing.level.list.full<-missing.level.list
			
			#determine if all or all but one level of a categorical variable is missing data 
			all.missing.level.vec<-allMissingLevels(missing.level.list, data)
			all.missing.varnames<-names(all.missing.level.vec)[all.missing.level.vec]
			
			#if there are categorical variables with completely, or all but one level missing data, 
			#remove that variable from the model. The effects of these variables cannot be estimated 
			#from the data, but that may not be clear from the optimization function
			if (any(all.missing.level.vec)){
			
				for (missing.cat.var in all.missing.varnames){
				
					#update the mxtrModel and fullModel objects 
					mxtrModel<-removeAllMissingCatVar(missing.cat.var, mxtrModel)
					
					if (!is.null(fullModel)){
					
						fullModel<-removeAllMissingCatVar(missing.cat.var, fullModel)
					
					}
					
					#print warning that the model was updated 
					warning.message<-paste0("There is not enough data to estimate the effect of ", missing.cat.var, " on ", var.name)
					warning(warning.message)
				
				}
				
				#after updating the models, remove the completely missing predictors from the list 
				#of predictors with at least one missing level 
				missing.level.list<-missing.level.list[! names(missing.level.list) %in% all.missing.varnames]
			
			}
			
			#for categorical variables having some but not all levels with completely missing data, remove missing levels  
			if (length(missing.level.list)>0){
			
				clean.data<-removeMissingLevels(missing.level.list, data)
			
			} else if (length(missing.level.list)==0) {
			
				clean.data<-data
			
			}
			
			#if the Tvals are provided, and they are batch Tvals, make sure all batches are present 
			if (any(grep("batch", names(Tvals)))){
			
				#get the name of the batch variable
				batch.varname<-strsplit(names(Tvals)[1], " ")[[1]][2]
				
				#first, check to determine if batch was previously removed from the model and set Tvals to NULL if so 
				if (batch.varname %in% all.missing.varnames){
				
					Tvals<-NULL
				
				#otherwise, if batch is missing a level, then remove the corresponding Tvals
				} else if (batch.varname %in% names(missing.level.list)){
				
					present.batches<-clean.data$batch
					present.batch.tval.names<-paste("batch", batch.varname, present.batches)
					Tvals<-Tvals[match(present.batch.tval.names,  names(Tvals))]
		
				}
			
			}
			
			#re-order data frame so we get correct reference levels for categorical variables 
			for (cat.var in cat.vars){
			
				current.cat.var.levels<-levels(clean.data[ , cat.var])
				ordered.current.levels<-cat.var.levels[[cat.var]][cat.var.levels[[cat.var]] %in% current.cat.var.levels]
				clean.data[ , cat.var]<-factor(clean.data[ , cat.var], levels=ordered.current.levels)
			
			}
						
			#run mixture model 
			mixmod.results<-runMxtrmod(ynames=var.name, mxtrModel=mxtrModel, Tvals=Tvals, nNA=nNA, minProp=minProp, method=method, data=clean.data, fullModel=fullModel)
					
			#add warning note about excluded levels
			warning.note<-c()
			for (element in names(missing.level.list.full)){
			
				if (! element %in% all.missing.varnames){
				
					missing.values<-missing.level.list[[element]]
					missing.value.note<-paste(element, paste(missing.values, collapse=", "), sep=":")
					warning.note<-c(warning.note, missing.value.note)
				
				} 
			
			}
			
			if (length(warning.note>0)){
			
				mixmod.results$predictors_missing_levels<-paste(warning.note, collapse="; ")
			
			}
					
			#add note about excluded variables 
			if (any(names(missing.level.list.full) %in% all.missing.varnames)){
			
				excluded.vars<-paste(names(missing.level.list.full)[names(missing.level.list.full) %in% all.missing.varnames], collapse=",")
				mixmod.results$excluded_predictors<-excluded.vars
			
			}		
		
		} else {
		
			warning(paste("There is not enough data to run the mixutre model for", ynames))
			return(NULL)
		
		}

	} else {
		
		clean.data<-data
	
		#make sure we still have the correct reference levels for categorical predictors 
		for (cat.var in cat.vars){
			
				current.cat.var.levels<-levels(clean.data[ , cat.var])
				ordered.current.levels<-cat.var.levels[[cat.var]][cat.var.levels[[cat.var]] %in% current.cat.var.levels]
				clean.data[ , cat.var]<-factor(clean.data[ , cat.var], levels=ordered.current.levels)
			
		}
		
		#make sure Tvals are still appropriate 
		if (any(grep("batch", names(Tvals)))){
		
			batch.varname<-strsplit(names(Tvals)[1], " ")[[1]][2]
			present.batches<-clean.data$batch
			present.batch.tval.names<-paste("batch", batch.varname, present.batches)
			Tvals<-Tvals[match(present.batch.tval.names,  names(Tvals))]

		}
		
		#then run mixture model 
		mixmod.results<-runMxtrmod(ynames=var.name, mxtrModel=mxtrModel, Tvals=Tvals, nNA=nNA, minProp=minProp, method=method, data=clean.data, fullModel=fullModel)
	
	}
  
	return(mixmod.results)
  
  })
  
  #put together final results table 
  final.result<-rbind.fill(results.list)  
  final.result[ , ".id"]<-as.character(final.result[ , ".id"])
  return(final.result)
  
 }

