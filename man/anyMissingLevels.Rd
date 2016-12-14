\name{anyMissingLevels}
\alias{anyMissingLevels}
\title{
A function to determine whether any level of a categorical variable has completely missing outcome (metabolite) data 
}
\description{
This function determines whether any level of a categorical predictor variable has completely missing outcome data.
}
\usage{
anyMissingLevels( yname,  cat.vars,  dataset)
}
\arguments{
  \item{yname}{
A character string indicating the outcome variable, i.e., a metabolite.
}
  \item{cat.vars}{
A character vector of categorical variable names. 
}
  \item{dataset}{
A data frame containing the yname and cat.vars variables. 
}
}
\value{
 Returns a list indicating whether each of the variables specified in argument cat.vars has at least one level with completely missing values for the outcome variable specified in yname.

}

\author{
Michael Nodzenski
}

\examples{

#create example analysis data 
data(euMetabCData)
example.data<-euMetabCData
example.data[example.data$batch==1, "aspartic_acid"]<-NA

#check to determine if aspartic acid values are entirely missing for 
#any level of batch or pheno 
anyMissingLevels( "aspartic_acid", c( "batch",  "pheno"), example.data)


}

