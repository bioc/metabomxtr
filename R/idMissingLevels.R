#This function determines the specific level of a categorical variable with completely 
#missing outcome data, given that the categorical variable is known to be missing outcome data for at 
#least one level    

idMissingLevels<-function(yname, missing.levels.list, dataset){

	#get the element of the missing levels list corresponding to yname
	target.list<-missing.levels.list[[yname]]
	
	#identify the variables with any missing levels 
	missing.level.vars<-names(target.list)[target.list]
	
	#loop over those variables and pull out the actual values of the missing levels 
	missing.level.vals<-lapply(missing.level.vars, function(missing.var){
	
		yname.vec<-dataset[ , yname]
		subset.data<-dataset[!is.na(yname.vec), missing.var, drop=F]
		susbet.levels<-unique(subset.data[ , missing.var])
		original.levels<-unique(dataset[ , missing.var])
		missing.levels<-original.levels[! original.levels %in% susbet.levels]
		return(missing.levels)
	
	
	})
	names(missing.level.vals)<-missing.level.vars
	return(missing.level.vals)

}