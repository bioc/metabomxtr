#This function determines whether any level of a categorical variable has entirely missing values  

anyMissingLevels<-function(yname, cat.vars, dataset){

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