\name{mxtrmodstart}
\alias{mxtrmodstart}
\title{
A function to generate starting parameter estimates for the optimization of mixture model parameters.
}
\description{
This function returns starting parameter estimates for the optimization of the mixture model parameters. The intercept of the continuous portion is set to the mean of the observed responses and the intercept of the discrete portion is set to the log odds of having observed a response. All other parameter starting values are set to zero. 
}
\usage{
mxtrmodstart(obsY, xVars, zVars, includeDiscrete)
}
\arguments{
  \item{obsY}{
A vector containing the response variable, which must be normally or log normally distributed. If lognormal, log transformed Y's should be input as the response variable. Missing Y values should be indicated by NA. 
}
  \item{xVars}{
The design matrix for the covariates included in the discrete portion of the model. 
}
  \item{zVars}{
The design matrix for the covariates inlcuded in the continuous portion of the model. 
}
  \item{includeDiscrete}{
A logical indicator for whether or not to include the discrete portion of the model.
}
}
\value{
A vector containing the starting values for each parameter in the mixture model function, to be used as starting points when optimizing the parameter estimates. 
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
testdata<-data.frame(cbind(yvar,y2var,logyvar,logy2var,pred1,pred2))

#Create a vector of response variables
obsY<-testdata$logyvar

#Create the design matrix for the discrete portion of the model
xVars<-model.matrix(~pred1+pred2,data=testdata)

#Create the design matrix for the continuous portion of the model 
zVars<-model.matrix(~pred1+pred2,data=testdata)

#Determine if the discrete portion should be inlcuded in the model 
includeDiscrete<-sum(is.na(obsY))>5

#Calculate starting values
mxtrmodstart(obsY=obsY,xVars=xVars,zVars=zVars,includeDiscrete=includeDiscrete)
}


