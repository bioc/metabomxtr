mxtrmod <-function(ynames,mxtrModel,Tvals=NULL,nNA=5,minProp=0.2,method="BFGS",data,fullModel=NULL){
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
  #  checkMissing: Should factor/character variables in the input mixture model be checked to determine 
  #				   if sufficient data exists to estimate their effects? Defualts to TRUE. 
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
  
  #Make the model a Formula object   
  mxtrModelF<-Formula(mxtrModelTarget)
  
  #in the input data are expression sets, convert to data frames 
  if (class(data)=="ExpressionSet"){

	  data<-data.frame(cbind(t(exprs(data)), pData(data)))

  }

  #identify categorical variables in model 
  model.vars<-all.vars(mxtrModelF)
  model.var.classes<-sapply(model.vars,function(x){ class(data[ ,x])})
  cat.vars<-names(model.var.classes)[model.var.classes %in% c("factor", "character")]
  
  #determine whether any of the levels of categorical variables are entirely missing
  missing.levels.check<-lapply(ynames, anyMissingLevels, cat.vars=cat.vars, dataset=data)
  names(missing.levels.check)<-ynames
  any.missing.levels<-sapply(missing.levels.check, any)
  missing.level.vars<-names(any.missing.levels)[any.missing.levels]
  no.missing.level.vars<-names(any.missing.levels)[!any.missing.levels]

  #if there are any variables entirely missing for at least one level of the categorical variables, 
  #remove those levels from the data and run a mixture model
  if (length(missing.level.vars)>0){
  
	missing.level.res.list<-lapply(missing.level.vars, function(var.name){

		#identify actual levels that are missing 
		missing.level.list<-idMissingLevels(var.name, missing.levels.check, data)
		
		#remove those levels 
		clean.data<-removeMissingLevels(var.name, missing.level.list, data)
		
		#if the Tvals are provided, and they are batch Tvals, make sure all batches are present 
		if (any(grep("batch", names(Tvals)))){
		
			#get the name of the batch variable
			batch.varname<-strsplit(names(Tvals)[1], " ")[[1]][2]
			
			#if batch is missing a level, then remove the corresponding Tvals
			if (batch.varname %in% names(missing.level.list)){
			
				batch.missing.levels<-missing.level.list[[batch.varname]]
				batch.missing.tval.names<-paste("batch", batch.varname, batch.missing.levels)
				Tvals<-Tvals[! names(Tvals) %in% batch.missing.tval.names]
			
			}
		
		}
		
		#run mixture model 
		mixmod.results<-runMxtrmod(ynames=var.name, mxtrModel=mxtrModel, Tvals=Tvals, nNA=nNA, minProp=minProp, method=method, data=clean.data, fullModel=fullModel)
		
		#add warning note about excluded levels
		warning.note<-c()
		for (element in names(missing.level.list)){
		
			missing.values<-missing.level.list[[element]]
			missing.value.note<-paste(element, paste(missing.values, collapse=", "), sep=":")
			warning.note<-c(warning.note, missing.value.note)
		
		}
		mixmod.results$missing_levels_warning<-paste(warning.note, collapse=", ")
		return(mixmod.results)
		
	})
	
  }  
  
  if (length(no.missing.level.vars)>0){
  
  
	no.missing.results<-runMxtrmod(ynames=no.missing.level.vars, mxtrModel=mxtrModel, Tvals=Tvals, nNA=nNA, minProp=minProp, method=method, data=data, fullModel=fullModel)
  
  }
  
  #put together final results table 
  if (exists('no.missing.results') & exists('missing.level.res.list')){
  
	final.result<-rbind.fill(c(list(no.missing.results), missing.level.res.list))  
  
  } else if (exists('missing.level.res.list')){
    
	final.result<-rbind.fill(missing.level.res.list) 
  
  } else if (exists('no.missing.results')){
  
	final.result<-no.missing.results
  
  }
  
  return(final.result)
  
 }
