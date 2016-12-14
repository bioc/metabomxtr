addBatchMeans<-function(metab.col, batch.col, in.data, sample.type.name, analytical.type.name, color.col=NULL){

	###############################################################################################
	#This function takes a data frame of metabolite abundances and batch information, 
	#calculates mean metabolite abundance per batch, and appends the means to the 
	#bottom of the input data. 
	#
	#Arguments: 
	#
	# metab.col: A character value indicating the name of the variable in the input data set
	#			 that corresponds to metabolite abundance. 
	#
	# batch.col: A character value indicating the name of the variable in the input data set 
	#			 corresponding to batch. 
	#
	# in.data: A data frame containing the variables specified in metab.col, batch.col, and color.col
	#
	# sample.type.name: A character value indicating the type of sample: "Quality Control" or "Experimental"
	#
	# analytical.type.name: A character value indicating whether the data have been normalized. Either 
	#						"Before Normalization" or "After Normalization". 
	#
	# color.col: A character value indicating the name of the variable in the input data set 
	#			 corresponding to a grouping variable (for instance mom vs. baby samples). 
	#
	#Returns:
	#	A data frame of metabolite abundances, with batch specific means appended to the end. 
	######################################################################################################
	
	#grab target columns 
	in.data<-as.data.frame(in.data)
	in.data$sample.type<-sample.type.name
	in.data$analytical.type<-analytical.type.name
	target.cols<-c(metab.col, batch.col, 'sample.type', 'analytical.type', color.col, 'outlier')
	data.target<-in.data[ , target.cols]
	data.target$point.type<-"Observed Data"
	colnames(data.target)[1]<-"metabolite.abundance"
	
	#get rid of infinite rows
	any.non.inf.non.outlying.data<-data.target[!is.infinite(data.target$metabolite.abundance) & !data.target$outlier, ]
	
	if (nrow(any.non.inf.non.outlying.data)>0){
	
		#to get means more easily, set infinite values to NA
		non.inf.data<-data.target
		non.inf.data[is.infinite(non.inf.data$metabolite.abundance), 'metabolite.abundance']<-NA
		
		#also set outliers to NA
		non.inf.data[non.inf.data$outlier, 'metabolite.abundance']<-NA
	
		#get batch specific means
		mean.annotation<-ddply(non.inf.data, as.quoted(c(batch.col, color.col)), summarize, metabolite.abundance=mean(substitute(metabolite.abundance),na.rm=T))
		mean.annotation$sample.type<-sample.type.name
		mean.annotation$analytical.type<-analytical.type.name
		mean.annotation$point.type<-"Batch Specific Mean"
		mean.annotation$outlier<-FALSE
	
		#merge into a single df 
		target.data.cols<-colnames(data.target)
		final.data<-rbind(data.target, mean.annotation[ , target.data.cols])
		return(final.data)
	
	}

}
