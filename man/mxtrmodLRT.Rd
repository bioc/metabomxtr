\name{mxtrmodLRT}
\alias{mxtrmodLRT}
\title{
A function to run likelihood ratio tests on full vs. reduced mixture models.
}
\description{
This function runs likelihood ratio tests on full vs. reduced mixture models. Input arguments are data frame outputs from the mxtrmod function.  
}
\usage{
mxtrmodLRT(fullmod, redmod, adj = NULL)
}
\arguments{
  \item{fullmod}{
The output data frame from the mxtrmod function on the full mixture model. 
}
  \item{redmod}{
The output data frame from the mxtrmod function on the reduced mixture model. 
}
  \item{adj}{
The adjustment method for multiple comparisons. The default is set to NULL. Options for adjustment methods are described in the documentation for the function mt.rawp2adjp in the multtest package. 
}
}
\value{
A data frame containing the response variables (i.e. metabolites), negative log likelihoods of full and reduced models, chi square statistics, degrees of freedom, p-values, and, if requested, adjusted p-values.    
}
\references{
Moulton LH, Halsey NA. A mixture model with detection limits for regression analyses of antibody response to vaccine. Biometrics. 1995 Dec;51(4):1570-8.
}
\author{
Michael Nodzenski, Anna Reisetter, Denise Scholtens
}
\examples{
#Create sample data
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

#Run a full mixture model on each response variable 
fullMod<-~pred1+pred2+pred3|pred1+pred2+pred3
fullModRes<-mxtrmod(ynames=ynames,mxtrModel=fullMod,data=testdata)
fullModRes

#Run a reduced mixture model on each response variable
redMod<-~pred2|pred2
redModRes<-mxtrmod(ynames=ynames,mxtrModel=redMod,data=testdata,fullModel=fullMod)
redModRes

#Compare models using likelihood ratio test
mxtrmodLRT(fullModRes,redModRes)

}

