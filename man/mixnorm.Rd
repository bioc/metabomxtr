\name{mixnorm}
\alias{mixnorm}
\title{
A function to perform per-metabolite batch normalization using a mixture model with batch-specific thresholds and run order correction if desired.
}
\description{
This function performs per-metabolite batch normalization using a mixture model with batch-specific thresholds and run order correction if desired.
}
\usage{
mixnorm(ynames, batch = "Batch", mxtrModel = NULL, batchTvals = NULL, correctSampleType = TRUE, sampleType = NULL, nNA = 5, minProp = 0.2, method = "BFGS", cData, data)
}
\arguments{
  \item{ynames}{
A character vector of the mixture model outcome names, e.g. metabolites. If the input data object is a matrix or data frame, these should be column names. If the input data object is an expression set, these should be row names. Response variables should have a normal or lognormal distribution. If lognormal,log transformed variables should be input. Missing values should be denoted by NA.
}
  \item{batch}{
A character value indicating the name of the variable in cData and data that indicates batch.
}
  \item{mxtrModel}{
A formula of the form ~x1+x2...|z1+z2..., where x's are the names of covariates included in the discrete portion of the model and z's are names of covariates included in the continuous portion. For intercept only models, enter 1 instead of covariate names on the appropriate side of the |. The covariate names must be the same for cData and data.  The default model includes a variable 'batch' for both discrete and continuous model components.
}
  \item{batchTvals}{
A vector, the length of batch, of thresholds below which continuous variables are not observable. The default is the minimum across all response variables (metabolites) in a given batch.
}
  \item{correctSampleType}{
By default, when there are control samples of different types, normalized values are returned with a location shift correction for control sample type. If correctSampleType=FALSE, then normalized values will be returned without this location shift correction.  In this case, sampleType must be specified and must be included in at least one portion of the mxtrModel formula. sampleType must be coded as a factor with the same levels in both cData and data.
}
  \item{sampleType}{
The name of the factor variable indicating the type of control sample. The variable must be the same in both cData and data.
}
  \item{nNA}{
The minimum number of unobserved values needed to be present for the discrete portion of the model likelihood to be calculated. Models for variables with fewer than nNA missing values will include only the continuous portion. The default value is 5.
}
  \item{minProp}{
The minimum proportion of non-missing data in the response variable necessary to run the model. The default value is 0.2. Models will not be run if more than 80\% of response variable values are missing.
}
  \item{method}{
The method used to optimize the parameter estimates of the mixture model. "BFGS" is the default method. Other options are documented in the manual for the function 'optimx' in package optimx.
}
  \item{cData}{
The input data object of control data to estimate normalization parameters. Matrices, data frames, and expression sets are all acceptable classes. If a data frame or matrix, rows are subjects and columns are metabolites or outcomes.
}
  \item{data}{
The input data object for observed values to be normalized (i.e. not controls). Matrices, data frames, and expression sets are all acceptable classes. If a data frame or matrix, rows are subjects and columns are metabolites or outcomes.
}
}
\details{
This function adapts the mxtrmod function in a normalization context in which aliquots from one or more control samples are run with each batch in a series of non-targeted metabolomics assays. The function accepts a data frame of log2 peak areas from control samples and a separate data frame of log2 peak areas from samples of analytical interest.  
}
\value{
Returns a list with the following components:
\item{normParamsZ}{A data frame of the per-metabolite parameter estimates from the mixture model that are subtracted from the observed values to created the normalized data set.}
\item{ctlNorm}{A data frame of normalized values for the control samples.}
\item{obsNorm}{A data frame of normalized values for the samples of analytical interest.}
}
\references{
Nodzenski M, Muehlbauer MJ, Bain JR, Reisetter AC, Lowe WL Jr, Scholtens DM. Metabomxtr: an R package for mixture-model analysis of non-targeted metabolomics data. Bioinformatics. 2014 Nov 15;30(22):3287-8.
}
\author{
Denise Scholtens, Michael Nodzenski, Anna Reisetter
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
                        data=euMetabData)
}

