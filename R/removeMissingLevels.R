#This function removes rows of the input dataset corresponding to 
#levels of factor/categorical variables with entirely missing outcome data, 
#and then re-levels the factor  

removeMissingLevels<-function(yname, missing.levels.list, dataset){

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