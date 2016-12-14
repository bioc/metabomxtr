metabplot<-function(metab.name, batch="Batch", raw.obs.data=NULL, raw.cont.data=NULL, norm.obs.data=NULL, norm.cont.data=NULL, color.var=NULL, cont.outlier.sd.thresh=2, norm.outlier.sd.thresh=4){

	#################################################################################################################
	#This function plots metabolite abudance data, both before and after normalization
	#with mixnorm. 
	#
	#Arguments:
	#
	#metab.name: A character value indicating the name of the variable in the input data sets
	#			 that corresponds to metabolite abundance. If more than one data set is input, 
	#			 the metab.name variable must be the same across all inputs. 
	#
	#batch: A character value indicating the name of the variable in the input data sets that indicates batch.
	#       If not specified, this argument defaults to "Batch". 
	#
	#raw.obs.data: A data frame of raw (non.normalized) experimental sample
	#				metabolite abundances and predictor variables. If input, 
	#				this must include the variables specified in metab.name, 
	#				batch, and color.var. 
	#
	#raw.cont.data: A data frame of raw (non.normalized) quality control sample
	#				metabolite abundances and predictor variables.If input, 
	#				this must include the variables specified in metab.name, 
	#				batch, and color.var.  
	#
	#norm.obs.data: A data frame or matrix of normalized experimental sample
	#				metabolite abundances output by function mixnorm. If raw.obs.data is specified, 
	#				the variables specified in batch and color.var do not need to be included in 
	#				this object, and it is assumed that the row order in norm.obs.data is identical 
	#				to raw.obs.data. If raw.obs.data is not specified, then this object must contain  
	#				the variables specified in batch and color.var. 
	#
	#norm.cont.data: A data frame or matrix of normalized quality control sample
	#				metabolite abundances output by function mixnorm. If raw.cont.data is specified, 
	#				the variables specified in batch and color.var do not need to be included in 
	#				this object, and it is assumed that the row order in norm.cont.data is identical 
	#				to raw.cont.data. If raw.cont.data is not specified, then this object must contain  
	#				the variables specified in batch and color.var. 
	#
	#color.var: A character value indicating the name of the variable that will be used to 
	#			define the color and grouping of data points in output plots. 
	#
	#cont.outlier.sd.thresh: The number of standard deviations from the mean metabolite abundance a point must be to be 
	#					considered an outlier for the raw (non-normalized) control data. This should match argument 
	#					qc.sd.outliers from function mixnorm. This defaults to 2 standard deviations, which is the same 
	#					as the default for qc.sd.outliers.
	#
	#norm.outlier.sd.thresh: The number of standard deviations from the mean metabolite abundance a point must be to be 
	#					considered an outlier for the normalized experimental data. This defaults to 4 standard deviations. 
	#Returns:
	# A plot of metabolite abundance by batch. 
	###############################################################################################################################################
	
	#make sure there is at least some data input 
	if (is.null(raw.obs.data) & is.null(raw.cont.data) & is.null(norm.obs.data) & is.null(norm.cont.data)){
	
		stop("No data were input to be plotted.")
	
	} 
	
	#get batch specific means and bind together input data 
	if (!is.null(raw.obs.data)){
	
		raw.obs.data$outlier<-FALSE
		raw.obs.batch.means.data<-addBatchMeans(metab.name, batch, raw.obs.data, sample.type.name="Experimental", analytical.type.name="Before Normalization", color.var)
		
	} else {
	
		raw.obs.batch.means.data<-NULL
	
	}
	
	if (!is.null(raw.cont.data)){
	
		raw.cont.data.with.outliers<-addOutlierInfo(metab.name, raw.cont.data, cont.outlier.sd.thresh)
		raw.cont.batch.means.data<-addBatchMeans(metab.name, batch, raw.cont.data.with.outliers, sample.type.name="Quality Control", analytical.type.name="Before Normalization", color.var)
		raw.cont.batch.means.data$point.type<-ifelse(raw.cont.batch.means.data$outlier, "Excluded Outlier", 
															raw.cont.batch.means.data$point.type)
	
	} else{
	
		raw.cont.batch.means.data<-NULL
	
	}
	
	if (!is.null(norm.obs.data)){
		
		norm.obs.data<-as.data.frame(norm.obs.data)
				
		#since the normalized data is from mixnorm, 
		#need to add on predictor data if not present 
		if (any(! c(batch, color.var) %in% colnames(norm.obs.data))){
		
			norm.obs.data<-cbind(norm.obs.data, raw.obs.data[ ,  c(batch, color.var), drop=F])		
		}
		
		norm.obs.data.with.outliers<-addOutlierInfo(metab.name, norm.obs.data, norm.outlier.sd.thresh)
		norm.obs.batch.means.data<-addBatchMeans(metab.name, batch, norm.obs.data.with.outliers, sample.type.name="Experimental", analytical.type.name="After Normalization", color.var)
		norm.obs.batch.means.data$point.type<-ifelse(norm.obs.batch.means.data$outlier, "Potential Outlier", 
															norm.obs.batch.means.data$point.type)
	
	} else {
	
		norm.obs.batch.means.data<-NULL
	
	}
	
	if (!is.null(norm.cont.data)){
	
		norm.cont.data<-as.data.frame(norm.cont.data)
				
		#since the normalized data is from mixnorm, 
		#need to add on predictor data if not present 
		if (any(! c(batch, color.var) %in% colnames(norm.cont.data))){
		
			norm.cont.data<-cbind(norm.cont.data, raw.cont.data[ ,  c(batch, color.var), drop=F])		
		}
		
		norm.cont.data$outlier<-FALSE
		norm.cont.batch.means.data<-addBatchMeans(metab.name, batch, norm.cont.data, sample.type.name="Quality Control", analytical.type.name="After Normalization", color.var)
	
	} else {
	
		norm.cont.batch.means.data<-NULL
	
	}
		
	#put into a list 
	data.list<-list(raw.obs.batch.means.data, raw.cont.batch.means.data, norm.obs.batch.means.data, norm.cont.batch.means.data)
		
	#put into a single data frame 
	all.data<-do.call("rbind", data.list)
	
	#re-level so before normalization is the reference
	all.data$analytical.type<-factor(all.data$analytical.type,levels=c("Before Normalization","After Normalization"))
	
	#set all infinite values to NA for plotting 
	all.data$metabolite.abundance[is.infinite(all.data$metabolite.abundance)]<-NA
	
	#plot data 
	
	#first, specify shapes for plot points 
	if (length(unique(all.data$point.type))==3 & "Excluded Outlier" %in% all.data$point.type){
	
		point.numbers<-c(8, 13, 16)
	
	} else if ((length(unique(all.data$point.type))==3 & "Potential Outlier" %in% all.data$point.type)){
	
		point.numbers<-c(8, 16, 4)
	
	} else if (length(unique(all.data$point.type))==4){ 
	
		point.numbers<-c(8, 13, 16, 4)
	
	} else {
	
		point.numbers<-c(8, 16)
	
	}
	if (!is.null(color.var)){
	
		ggplot(all.data, aes_string(x=batch, y='metabolite.abundance', color=color.var, shape='point.type', group=color.var))+
			geom_point() + geom_line(data=all.data[which(all.data$point.type=="Batch Specific Mean") , ]) +
			scale_color_discrete(name="Sample Type") + scale_shape_manual(name="Point Type", values=point.numbers) +
			xlab("Analytic Batch")+ylab(metab.name)+facet_grid(sample.type~analytical.type) + theme(axis.text.x=element_text(angle=90, hjust = 1))
	
	} else {
	
		ggplot(all.data, aes_string(x=batch, y='metabolite.abundance', shape='point.type', group=1))+
			geom_point() + geom_line(data=all.data[which(all.data$point.type=="Batch Specific Mean") , ]) + 
			scale_shape_manual(name="Point Type", values=point.numbers) +
			xlab("Analytic Batch")+ylab(metab.name)+facet_grid(sample.type~analytical.type)	+ theme(axis.text.x=element_text(angle=90, hjust = 1))
	
	}
	
}