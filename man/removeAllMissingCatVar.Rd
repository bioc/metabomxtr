\name{removeAllMissingCatVar}
\alias{removeAllMissingCatVar}
\title{
A function to remove categorical variables with insufficient data to estimate effects from mixture models
}
\description{
This function removes categorical variables from mixture models. This is needed when metabolite values are entirely missing or only present for one level of the categorical variable, making it impossible to estimate effects for the categorical variable. 
}
\usage{
removeAllMissingCatVar(cat.varname,  mxtrModel)
}
\arguments{
  \item{cat.varname}{
A character string indicating a categorical predictor variable in the mixture model.
}
  \item{mxtrModel}{
A mixture model formula of class formula.
}
}
\value{
 Returns a mixture model formula of class formula with the variable indicated in cat.varname omitted. 

}

\author{
Michael Nodzenski
}

\examples{

#specifiy mixture model 
mix.model<- ~ var1 + var2 | var1 + var2

#remove var1 
removeAllMissingCatVar("var1", mix.model)

}

