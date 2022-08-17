anyMissingLevels<-function(yname, cat.vars, dataset){

  ##############################################################
  #This function determines whether any level of a categorical 
  #variable in the mixture model has entirely missing metabolite
  #values.  
  #
  #Args:
  #	yname: A character string corresponding to the metabolite column 
  #		   name. 
  #	cat.vars: a character vector of the names of categorical 
  #			  variables in the mixture model 
  # dataset: a data frame containing metabolite levels and 
  #			 categorical predictors 
  #
  #Returns:
  #	A logical vector indicating whether each categorical variable 
  # has at least one level with entirely missing metabolite data 
  ####################################################################

  #first, start by removing rows that have all missing outcome values
  yname.vec<-dataset[ , yname]
  subset.data<-dataset[!is.na(yname.vec), cat.vars, drop=FALSE]
  
  #now, for each level of cat vars, determine whether there were any levels with entirely missing outcomes 
  missing.levels<-sapply(cat.vars, function(cat.varname){
  
	cat.var.levels<-unique(subset.data[ ,cat.varname])
	original.cat.var.levels<-unique(dataset[ ,cat.varname])
	if (length(cat.var.levels)==length(original.cat.var.levels)){
	
		FALSE
	
	} else {
	
		TRUE
	
	}
  })
  
  return(missing.levels)
  
}