\name{metabplot}
\alias{metabplot}
\title{
A function to plot metabolite abundance before and after normalization.
}
\description{
This function plots metabolite abundance before and after normalization. 
}
\usage{
metabplot(metab.name, batch="Batch", raw.obs.data=NULL, raw.cont.data=NULL, norm.obs.data=NULL, norm.cont.data=NULL, color.var=NULL, cont.outlier.sd.thresh=2, norm.outlier.sd.thresh=4)
}
\arguments{
  \item{metab.name}{
A character value indicating the name of the variable in the input data sets that corresponds to metabolite abundance. If more than one data set is input, the metab.name variable must be the same across all inputs. 
}
  \item{batch}{
A character value indicating the name of the variable in the input data sets that corresponds to batch. If not specified, this argument defaults to "Batch". 
}
  \item{raw.obs.data}{
A data frame of raw (non-normalized) experimental sample metabolite abundances and predictor variables. If input, this must include the variables specified in metab.name, batch, and color.var. 
}
  \item{raw.cont.data}{
A data frame of raw (non-normalized) quality control sample metabolite abundances and predictor variables. If input, this must include the variables specified in metab.name, batch, and color.var. 
}
  \item{norm.obs.data}{
A data frame or matrix of normalized experimental sample metabolite abundances output by function mixnorm. If raw.obs.data is specified, the variables specified in batch and color.var do not need to be included in this object, and it is assumed that the row order in norm.obs.data is identical to raw.obs.data. If raw.obs.data is not specified, then this object must contain the variables specified in batch and color.var. 
}
  \item{norm.cont.data}{
A data frame or matrix of normalized quality control sample metabolite abundances output by function mixnorm. If raw.cont.data is specified, the variables specified in batch and color.var do not need to be included in this object, and it is assumed that the row order in norm.cont.data is identical to raw.obs.data. If raw.cont.data is not specified, then this object must contain the variables specified in batch and color.var. 
}
  \item{color.var}{
A character value indicating the name of the variable that will be used to define the color and grouping of data points in output plots. 
}
  \item{cont.outlier.sd.thresh}{
The number of standard deviations from the mean a point must be to be considered an outlier for the raw quality control data. If plotting data output by function mixnorm, this should match the value for argument qc.sd.outliers in mixnorm. This argument defaults to 2, which is the same as the default for qc.sd.outliers. These outlying points will be indicated in the plots as "Excluded Outlier", and represent data points that were not included in estimating batch and other technical effects during normalization. 
}
  \item{norm.outlier.sd.thresh}{
The number of standard deviations from the mean a point must be to be considered an outlier for the normalized experimental data. This argument defaults to 4. These outlying points will be indicated in the plots as "Potential Outlier", and represent data points that remain outliers after normalization that the user may wish to remove before downstream analysis. 
}
}
\details{
This function aims to show a 4 panelled plot of quality control and experimental sample metabolite abundance data both before and after normalization wit h function mixnorm. However, it will produce a plot as long as at least one of raw.obs.data, raw.cont.data, norm.obs.data, or norm.cont.data is specified. The user may include as few as one or as many as all 4 of these arguments. 
}
\value{
Returns a graph of metabolite abudance.
}
\references{
Nodzenski M, Muehlbauer MJ, Bain JR, Reisetter AC, Lowe WL Jr, Scholtens DM. Metabomxtr: an R package for mixture-model analysis of non-targeted metabolomics data. Bioinformatics. 2014 Nov 15;30(22):3287-8.
}
\author{
Michael Nodzenski
}
\examples{
data(euMetabCData)
data(euMetabData)

ynames <- c("betahydroxybutyrate","pyruvic_acid","malonic_acid","aspartic_acid")

#in this example, batch minima specified in batchTvals were calculated from the full data set for this experiment that is not available here
euMetabNorm <- mixnorm(ynames,
                        batch="batch",
                        mxtrModel=~pheno+batch|pheno+batch,
                        batchTvals=c(10.76,11.51,11.36,10.31,11.90),
                        cData=euMetabCData,
                        data=euMetabData, 
                        qc.sd.outliers=2)
                        
 #plot results 
 metab.name.list<-names(euMetabNorm$obsNorm)	                
 plot.list<-lapply(metab.name.list, metabplot, batch="batch", raw.obs.data=euMetabData, raw.cont.data=euMetabCData, norm.obs.data=euMetabNorm$obsNorm, 	
 			norm.cont.data=euMetabNorm$ctlNorm, color.var="pheno" )
plot.list[[1]]
               
}

