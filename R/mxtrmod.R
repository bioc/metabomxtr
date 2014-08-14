mxtrmod <-
function(ynames,mxtrModel,Tvals=NULL,nNA=5,minProp=0.2,
                                      method="BFGS",data,fullModel=NULL){
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
  # Returns:
  #   A data frame with optimized estimates for all parameters in the mixture model, the negative 
  #   log likelihood, the optimization method used, whether the algorithm converged, and the 
  #  number of observations used.  
  #################################################################################################
  
  #Define a data frame containing the response variables
  obsY<-yvals(data,ynames)  

  #Make the model a Formula object 
  mxtrModel1<-Formula(mxtrModel)

  #Define the design matrix for the discrete portion of the model
  xVars<-xdesign(data,mxtrModel1)
  
  #Define the design matrix for the continuous portion of the model 
  zVars<-zdesign(data,mxtrModel1)
  
  #If mxtrModel is a reduced model, make the fullModel parameter a formula object. 
  #Then define a design matrices for the full version of the mixture model
  #and merge the two design matrices together. 
  if(!is.null(fullModel)){
    fullMod<-Formula(fullModel)
    xVarsFull<-xdesign(data,fullMod)
    zVarsFull<-zdesign(data,fullMod)
    allCovariates<-cbind(xVarsFull,zVarsFull)
   } else {
  #Otherwise, merge the original desgin matrices
  	allCovariates<-cbind(xVars,zVars)
   }
  
  #Determine if any of the observations have missing covariate values
  keep.these<-apply(allCovariates,1,function(x){all(!is.na(x))})
  
  #Remove observations from the data frame of response variables that 
  #do not have fully observed covariates
  obsY<-obsY[keep.these, ,drop=FALSE]
 
  #Remove observations from the discrete design matrix that 
  #do not have fully observed covariates
  xVars<-xVars[keep.these, ,drop=FALSE]
  
  #Remove observations from the continuous design matrix that 
  #do not have fully observed covariates
  zVars<-zVars[keep.these, ,drop=FALSE]

  #Run a mixture model on each response variable, to return optimized
  #parameter estimates and the negative log likelihood 
  result<-apply(obsY,2,function(obsY) {

  #Determine if there are enough missing values to include the discrete 
  #component of the mixture model.

  includeDiscrete <- sum(is.na(obsY))>nNA

  #Determine if there are enough non-missing values present to run the model 
  if(sum(!is.na(obsY))>minProp*length(obsY)){

  #If enough non-missing values are present, specify starting values for the model
  myStart<-mxtrmodstart(obsY=obsY,xVars=xVars,zVars=zVars,includeDiscrete=includeDiscrete)

  #Specify Tvals. By default, the function will use the minimum of the observed
  #responses.
  if(is.null(Tvals)){
    Tvals <- rep(min(obsY,na.rm=TRUE),length(obsY))
  }
  #Optimize the function.
  optimResult <- optimx(myStart,
                fn=mxtrmodLL,
                obsY=obsY,xVars=xVars,zVars=zVars,Tvals=Tvals,includeDiscrete=includeDiscrete,
		    method=method)

  #Get the optimized parameter estimates. 
  modelPars <- as.vector(coef(optimResult))

  #Name the parameters 
  if (dim(xVars)[2]==1){
    xVarNames<-NULL
  } else {
  xVarNames<-paste("x_",colnames(xVars)[-1],sep="")
  }
  if (dim(zVars)[2]==1){
    zVarNames<-NULL
  } else {
  zVarNames<-paste("z_",colnames(zVars)[-1],sep="")
  }
  #If the discrete component is included, name both continuous and discrete variables.
  namedPars<- modelPars
  if(includeDiscrete){
    names(namedPars) <- c("xInt",xVarNames,"zInt",zVarNames,"sigma")

  #Otherwise, only name the continuous variables. 
  } else if(!includeDiscrete){
    names(namedPars) <- c("zInt",zVarNames,"sigma")
  }   
  namedPars<-t(namedPars)
   
  #Calculate the negative log-likelihood of the model. 
  negLL <- mxtrmodLL(modelPars,obsY=obsY,xVars=xVars,zVars=zVars,Tvals=Tvals,includeDiscrete=includeDiscrete)
  
  #Combine the final results into a dataframe. 
  ans <- data.frame(namedPars, method=rownames(optimResult),
                    conv=optimResult$convcode,negLL=negLL,nObs=length(obsY))
  } else ans <- NULL
  ans
})
ldply(result,rbind)
}
