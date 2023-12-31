\name{mxtrmod}
\alias{mxtrmod}
\title{
A function to return optimized parameter estimates and the negative log-likelihood of mixture models for truncated normal or lognormal data 
}
\description{
This function returns optimized parameter estimates and the negative log-likelihood of mixture models for truncated normal or lognormal data, via call to function runMxtrmod, after properly accounting for factor variable predictors with entirely missing outcome data. 
}
\usage{
mxtrmod(ynames, mxtrModel, Tvals=NULL, nNA=5, minProp=0.2, method="BFGS", data, fullModel=NULL, remove.outlier.sd=NULL)
}
\arguments{
  \item{ynames}{
A character vector of the mixture model outcome names, e.g. metabolites. If the input data object is a matrix or data frame, these should be column names. If the input data object is an expression set, these should be row names. Response variables should have normal or lognormal distributions. If lognormal, log transformed variables should be input. Missing values should be denoted by NA.  
}
  \item{mxtrModel}{
A formula of the form ~x1+x2...|z1+z2..., where x's are the names of covariates included in the discrete portion of the model and z's are names of covariates included in the continuous portion. For intercept only models, enter 1 instead of covariate names on the appropriate side of the |. 
}
  \item{Tvals}{
A vector of thresholds below which continuous variables are not observable. By default, this parameter will be set to the minimum of the response variable. 
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
  \item{data}{
The input data object. Matrices, data frames, and expression sets are all acceptable classes. If a data frame or matrix, rows are subjects and columns are metabolites or outcomes.  
}
\item{fullModel}{
A formula of the form ~x1+x2...|z1+z2..., where x's are the names of covariates included in the discrete portion of the full model and z's are names of covariates included in the continuous portion. Input if the mxtrModel parameter represents a reduced model.  
}
\item{remove.outlier.sd}{
The maximum number of standard deviations from the mean that should be considered non-outlying metabolite abundance values. Metabolite abundances greater than this will be removed from modeling to avoid having outliers unduly influence model parameters. This defaults to NULL, meaning all data will be used in modeling. 
}
}
\value{
 Returns a data frame containing optimized estimates for all parameters in the mixture model, the negative log likelihood of the model, the optimization method used, the total number of observations used, whether the algorithm converged, whether any categorical predictor levels were not modeled, and whether any categorical predictors were omitted from the model entirely. Specific levels of categorical predictors will be excluded from models when metabolite data are entirely missing for that level. These variables and corresponding levels are reported in column "predictors_missing_levels". Categorical variables will be completely removed from models when metabolite values are missing entirely or present for only one level of the categorical variable. These omitted variables are reported in column "excluded_predictors". If those columns are not present in the output, no predictor exclusions occured. 

}
\references{
Moulton LH, Halsey NA. A mixture model with detection limits for regression analyses of antibody response to vaccine. Biometrics. 1995 Dec;51(4):1570-8.
Nodzenski M, Muehlbauer MJ, Bain JR, Reisetter AC, Lowe WL Jr, Scholtens DM. Metabomxtr: an R package for mixture-model analysis of non-targeted metabolomics data. Bioinformatics. 2014 Nov 15;30(22):3287-8.
}
\author{
Michael Nodzenski, Anna Reisetter, Denise Scholtens
}
\note{
This function may generate warning messages about production of NaNs, but the function is still operating normally. 
}
\examples{
#Create sample data frame 
set.seed(123)
yvar<-rlnorm(200)
these<-sample(1:100,20)
yvar[these]<-NA
logyvar<-log(yvar)
y2var<-rlnorm(200)
those<-sample(1:200,25)
y2var[those]<-NA
logy2var<-log(y2var)
pred1<-sample(0:1,200,replace=TRUE)
pred2<-sample(1:10,200,replace=TRUE)
pred3<-sample(0:1,200,replace=TRUE)
pred3miss<-sample(1:200,50)
pred3[pred3miss]<-NA
testdata<-data.frame(cbind(yvar,y2var,logyvar,logy2var,pred1,pred2,pred3))

#Get the names of the response variables 
ynames<-names(testdata)[3:4]

#Run a mixture model on each response variable 
mod<-~pred1+pred2+pred3|pred1+pred2+pred3
mxtrmod(ynames=ynames,mxtrModel=mod,data=testdata)

#Create example expression set
#Specify the response variables
exprsobs<-t(testdata[,3:4])

#Specify the phenotype data 
exprspheno<-testdata[,5:7]

#make phenotype data an annotated data frame 
phenoData <- new("AnnotatedDataFrame",data=exprspheno)

#combine into example expression set 
testexpr<-ExpressionSet(assayData=exprsobs,phenoData=phenoData)

#Get the names of the response variables 
ynames<-rownames(exprs(testexpr))

#Run the mixture model on each response variable 
mxtrmod(ynames=ynames,mxtrModel=mod,data=testexpr)

#Load the data set from the package
data(metabdata)

#Select the response variables
ynames<-names(metabdata)[11:17]

#Run the mixture models 
mod2<-~PHENO|PHENO+age_ogtt_mc+parity12+ga_ogtt_wks_mc
mxtrmod(ynames,mxtrModel=mod2,data=metabdata)
}

