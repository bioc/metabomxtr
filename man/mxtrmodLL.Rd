\name{mxtrmodLL}
\alias{mxtrmodLL}
\title{
A function to return the negative log-likelihood of mixture models of truncated, normal or lognormal data 
}
\description{
This function returns the negative log-likelihood of the specified mixture model. 
}
\usage{
mxtrmodLL(params, obsY, xVars, zVars, Tvals, includeDiscrete)
}
\arguments{
  \item{params}{
A vector of parameter estimates for all paramters in the mixture model. 
}
  \item{obsY}{
A vector containing the response variable, which must be normally or log normally distributed. If lognormal, log transformed Y's should be input as the response variable. Missing Y values should be indicated by NA. 
}
  \item{xVars}{
The design matrix for the covariates included in the discrete portion of the model.  
}
  \item{zVars}{
The design matrix for the covariates inlcuded in the continuous portion of the model. 
}
  \item{Tvals}{
A vector of thresholds below which continuous variables are not observable. 
}
  \item{includeDiscrete}{
A logical indicator for whether or not to include the discrete portion of the model.
}
}
\value{
Returns the negative log-likelihood of the specified mixture model. 
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

#Create a vector of starting values for the function 
startvals<-c(2,0,0,1.5,0,0,2)

#Create a vector of response variables
obsY<-testdata$logyvar

#Create the design matrix for the discrete portion of the model
xVars<-model.matrix(~pred1+pred2,data=testdata)

#Create the design matrix for the continuous portion of the model 
zVars<-model.matrix(~pred1+pred2,data=testdata)

#Create the Tvals vector 
Tvals<-rep(min(obsY,na.rm=TRUE),length(obsY))

#Determine if the discrete portion should be inlcuded in the model 
includeDiscrete<-sum(is.na(obsY))>5

#Calculate the negative log-likelihood
mxtrmodLL(params=startvals,obsY=obsY,xVars=xVars,zVars=zVars,Tvals=Tvals,
          includeDiscrete=includeDiscrete)
}

