mxtrmodLRT<-function(fullmod,redmod,adj=NULL){
   #########################################################################################
  # Runs likelihood ratio tests on full vs. reduced mixture models.  
  #
  # Args: 
  #  fullmod:  The output data frame produced by the mxtrmod function for the full model. 
  #  redmod:   The output data frame produced by the mxtrmod function for the reduced model.               
  #     adj:   The adjustment method for multiple comparisons. The default is set to NULL. 
  #            Options for adjustment methods are described in the documentation for the function 
  #            mt.rawp2adjp in the multtest package. 
  # Returns:
  #   A data frame containing the response variables (i.e. metabolites), negative log 
  #   likelihoods of full and reduced models, chi square statistics, degrees of freedom, 
  #   p-values, and, if requested, adjusted p-values.    
  ###########################################################################################
  
  #Error Checking 

  #check that fullmod and redmod have the same response variables
  if (!all(fullmod$.id==redmod$.id)){
  stop('fullmod and redmod do not have the same response variables')   
  }

  #check that redmod has fewer columns (is a reduced version) than fullmod
  if(dim(fullmod)[2]<=dim(redmod)[2]){
  stop('redmod is not a reduced version of fullmod')
  }

  #check that optimization algorithms are the same
  if (!all(fullmod$method==redmod$method)){
  stop('models use different optimization methods')
  }
  
  #check that the full and reduced models have the same number of observations
  if (!all(fullmod$nObs==redmod$nObs)){
  stop('at least one of the full models does not have the same number of observations as the reduced model')
  }
  
  #Run likelihood ratio tests

  #get the response variable names 
  ynames<-fullmod$.id	

  #get full model negative log likelihood 
  fullmodnegLL<-fullmod$negLL

  #get reduced model negative log likelihood
  redmodnegLL<-redmod$negLL

  #get the chi square value for the likelihood ratio test
  chisqval<-2*(redmodnegLL-fullmodnegLL)

  #get the number of parameters in the reduced model 
  redparms<-apply(redmod,1,function(x){length(x)-sum(is.na(x))-4})

  #get the number of parameters in the full model 
  fullparms<-apply(fullmod,1,function(x){length(x)-sum(is.na(x))-4})

  #get the degrees of freedom 
  degfr<-fullparms-redparms

  #find p-vals
  chisqP<-pchisq(chisqval,df=degfr,lower.tail=FALSE)

  #find adjusted p-values, if desired
  if (!is.null(adj)){
  adjp<-mt.rawp2adjp(chisqP,proc=adj)
  adjp<-adjp$adjp[order(adjp$index),2]
  } else {
  adjp<-NULL
  }

  #combine results into data frame

  if (!is.null(adjp)){
  results<-data.frame(.id=ynames,negLLFull=fullmodnegLL,negLLRed=redmodnegLL,
			     chisq=chisqval,df=degfr,p=chisqP,adjP=adjp)
  } else {
  results<-data.frame(.id=ynames,negLLFull=fullmodnegLL,negLLRed=redmodnegLL,
			     chisq=chisqval,df=degfr,p=chisqP)
	}
  #output result
  results
}
	