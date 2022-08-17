\name{mixnorm}
\alias{mixnorm}
\title{
A function to perform per-metabolite batch normalization using a mixture model with batch-specific thresholds and run order correction if desired.
}
\description{
This function performs per-metabolite batch normalization using a mixture model with batch-specific thresholds and run order correction if desired.
}
\usage{
mixnorm(ynames, batch = "Batch", mxtrModel=NULL, cData, data, batchTvals = NULL, removeCorrection=NULL, nNA = 5, minProp = 0.2, method = "BFGS", qc.sd.outliers=2 )
}
\arguments{
  \item{ynames}{
A character vector of the mixture model outcome names, e.g. metabolites. If the input data object is a matrix or data frame, these should be column names. If the input data object is an expression set, these should be row names. Response variables should have a normal or lognormal distribution. If lognormal,log transformed variables should be input. Missing values should be denoted by NA.
}
  \item{batch}{
A character value indicating the name of the variable in cData and data that indicates batch. If not specified, this argument defaults to "Batch"'.
}
  \item{mxtrModel}{
A formula of the form ~x1+x2...|z1+z2..., where x's are the names of covariates included in the discrete portion of the model and z's are names of covariates included in the continuous portion. The covariate names must be the same for cData and data.  The default model includes a variable specified in argument 'batch' for both discrete and continuous model components.If manually specified, mixture models must include at minimum batch in both the continuous and discrete portions. Models with covariates containing missing values will not run. See documentation for mxtrmod for additional details. 
}
  \item{cData}{
The input data object of control data to estimate normalization parameters. Matrices, data frames, and expression sets are acceptable classes. If a data frame or matrix, rows are subjects and columns are metabolites or outcomes.
}
  \item{data}{
The input data object for observed values to be normalized (i.e. not controls). Matrices, data frames, and expression sets are acceptable classes. If a data frame or matrix, rows are subjects and columns are metabolites or outcomes.
}
  \item{batchTvals}{
A vector of thresholds below which continuous variables are not observable. If specifying this argument, it must be the length of the unique levels of batch, with the thresholds in the same order that batch levels appear in cData. For instance, if there are 5 batches in cData, and they appear in order from 1 to 5, then argument batchTvals would need to be a numeric vector of length 5 with batch specific thresholds specified in that order. The default is the minimum across all response variables (metabolites) for each batch. ***Note that this default behavior may yield different normalization results for the same metabolite depending on the other metabolites entered as part of argument ynames. This should not be a concern if using all metabolite names as part of the ynames argument***. If running normalization on a subset of the full list of metabolites, to get the same results as the full set, the user should manually enter the minimum observed abundance across metabolites for each batch. 
}
  \item{removeCorrection}{
A character vector of variable names from mxtrModel whose effects should be estimated, but not subtracted from the non-normalized data. This parameter may be useful when data sets contain control samples of different types, for instance mothers and babies. In those instances, sample type may be an important covariate with respect to accurately estimating batch effects, necessitating inclusion in the mixture model, but it may not be of interest to actually subtract the estimated sample effect from the non-normalized data. If not specified, all estimated effects from the mixture model will be subtracted from the non-normalized data.
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
  \item{qc.sd.outliers}{
The maximum number of standard deviations from the mean that should be considered non-outlying metabolite abundance values for the control data. Metabolite abundances greater than this will be removed from modeling batch effects to avoid having outliers unduly influence normalization parameters. This defaults to 2 standard deviations from the mean. To not exclude outliers, enter Inf for this argument. 
}
}
\details{
This function adapts the mxtrmod function in a normalization context in which aliquots from one or more control samples are run with each batch in a series of non-targeted metabolomics assays. The function accepts a data frame of log2 peak areas from control samples and a separate data frame of log2 peak areas from samples of analytical interest.  
}
\value{
Returns a list with the following components:
\item{normParamsZ}{A data frame of the per-metabolite parameter estimates from the mixture model that are subtracted from the observed values to created the normalized data set. If  a parameter estimated is NA, that parameter level was used as the reference. This occurs when the outcome (metabolite) values are completely missing for a particular level of a categorical variable, and the specific missing level can be found in the conv element of the function outcome.}
\item{ctlNorm}{A data frame of normalized values for the control samples. Observations that were NA in the input data will remain NA, but observations considered outliers based on argument qc.sd.outliers  will not be normalized and instead coded as Inf.}
\item{obsNorm}{A data frame of normalized values for the samples of analytical interest. Observations that were NA in the input data will remain NA, but observations that could not be normalized due to too much missing quality control data, or too many outlying values in the quality control data, will not be normalized and instead coded Inf. This may occur if there is not enough data to estimated a specific batch effect (i.e., the effect of one out of 5 batches could not be estimated) or if there is not enough data to estimate the effect of a categorical predictor (i.e., if there is only QC data present for maternal samples, we can't estimate the effect relative to baby samples). }
\item{conv}{A data frame indicating whether models converged (indicated by a 0). This also indicates whether any categorical predictor levels were not modeled, and whether any categorical predictors were omitted from the model entirely. Specific levels of categorical predictors will be excluded from models when quality control metabolite data are entirely missing for that level. These variables and corresponding levels are reported in column "predictors_missing_levels". Categorical variables will be completely removed from models when quality control metabolite values are missing entirely or present for only one level of a categorical variable. These omitted variables are reported in column "excluded_predictors". Both columns indicate variables, or specific levels of variables, whose effects could not be accounted for in normalization. Normalized values will therefore not be output for metabolites with these warnings. For example, if metabolites were assayed in 20 total batches, and quality control data were completely missing for batch 2,  the effect of batch 2 cannot be estimated but the effects of the other 19 batches can and will be estimated. In the normalized data, non-missing metabolite values in batch 2 will be re-coded as Inf, but the correction will be made for the remaining batches. If columns predictors_missing_levels and excluded_predictors are not present, then no exclusions were made.}
}
\references{
Moulton LH, Halsey NA. A mixture model with detection limits for regression analyses of antibody response to vaccine. Biometrics. 1995 Dec;51(4):1570-8.
Nodzenski M, Muehlbauer MJ, Bain JR, Reisetter AC, Lowe WL Jr, Scholtens DM. Metabomxtr: an R package for mixture-model analysis of non-targeted metabolomics data. Bioinformatics. 2014 Nov 15;30(22):3287-8.
Reisetter AC, Muehlbauer MJ, Bain JR, Nodzenski M, Stevens RD, Ilkayeva O, Metzger BE, Newgard CB, Lowe WL Jr, Scholtens DM. Mixture model normalization for non-targeted gas chromatography/mass spectrometry metabolomics data. BMC Bioinformatics. 2017 Feb 2;18(1):84.
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

