addOutlierInfo<-function(metab.col, in.data, outlier.sd.thresh=2){
	
	##################################################################################
	#This function adds a column to a dataframe of metabolite 
	#abundance values indicating whether a observations
	#for a specific metabolite are outliers 
	#
	#Arguments:
	#
	#metab.col: A character value indicating the name of the target variable 
	#				in the input data set that corresponds to metabolite abundance.
	#
	#in.data: a data frame containing the variable specified in argument metab.col
	#
	#outlier.sd.thresh: The number of standard deviations from the mean metabolite 
	#					abundance an observation must be to be considered an outlier. 
	###########################################################################################
	
	#make sure input is a data frame 
	in.data<-as.data.frame(in.data)
	
	#get metabolite mean and standard deviation 
	metab.vec<-in.data[ , metab.col]
	non.inf.metab.vec<-metab.vec[!is.infinite(metab.vec)]
	mean.val<-mean(non.inf.metab.vec, na.rm=TRUE)
	sd.val<-sd(non.inf.metab.vec, na.rm=TRUE)
	
	#identify outliers
	outlier.rows<- !is.na(metab.vec) & !is.infinite(metab.vec) & 
								(metab.vec > mean.val + outlier.sd.thresh*sd.val | metab.vec < mean.val - outlier.sd.thresh*sd.val)
	
	#add in outlier data 
	in.data$outlier<-outlier.rows	
	return(in.data)
	
}
