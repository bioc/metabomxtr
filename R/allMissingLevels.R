allMissingLevels<-function(missing.levels.list, dataset){

	###################################################################
	#This function determines whether a categorical variable is missing 
	#metabolite data for all or all but one level, and thus whether 
	#that variable needs to be removed from the model. 
	#
	#Args:
	#	missing.levels.list: a list output by function idMissingLevels 
	#						 indicating which categorical variable
	#						 level predictors have no corresponding 
	#						 metabolite data, for a single metabolite. 
	#	dataset: a data frame containing the categorical predictors and
	#			 metabolite values. 
	#
	#Returns:
	#	A list indicating whether categorical model variables 
	#	have entirely missing metabolite data for all, or all 
	#	but one level 
	########################################################################

	all.missing.list<-unlist(lapply(names(missing.levels.list), function(cat.var){
	
		#get missing levels
		missing.levels<-missing.levels.list[[cat.var]]
		
		#Determine how many levels of the categorical variable have data present 
		cat.var.vec<-dataset[ ,cat.var]
		present.levels<-unique(cat.var.vec[! cat.var.vec %in% missing.levels])
		
		#if all levels are missing data, or all but one, return true, otherwise false 
		if (length(present.levels)<=1){
		
			TRUE
		
		} else if (length(present.levels)>2){
		
			FALSE
		
		}
		
	}))
	names(all.missing.list)<-names(missing.levels.list)
	return(all.missing.list)

}