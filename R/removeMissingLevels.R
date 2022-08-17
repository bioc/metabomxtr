removeMissingLevels<-function(missing.levels.list, dataset){

	######################################################################################
	#This function removes rows of the input dataset corresponding to 
	#levels of factor/categorical variables with entirely missing outcome data, 
	#and then re-levels the factor.
	#
	#Args:
	#	missing.levels.list: a list output by function idMissingLevels
	#						indicating which levels of categorical are completely 
	#						missing metabolite values for a specific metabolite 
	#	dataset: a dataframe containing the categorical variables in missing.levels.list
	#
	#Returns:
	#	A dataset with missing levels of categorical variables removed, and the factor 
	#	re-leveled based on the existing data 
	########################################################################################

	for (cat.var in names(missing.levels.list)){
	
		#get missing levels
		missing.levels<-missing.levels.list[[cat.var]]
		
		#remove those levels from the input data 
		cat.var.vec<-dataset[ ,cat.var]
		dataset<-dataset[! cat.var.vec %in% missing.levels, ]
		
		#re-level the factor 
		present.levels<-unique(cat.var.vec[! cat.var.vec %in% missing.levels])
		dataset[ ,cat.var]<-factor(dataset[ ,cat.var], levels=present.levels)
		
	}

	#return the subset data 
	return(dataset)

}