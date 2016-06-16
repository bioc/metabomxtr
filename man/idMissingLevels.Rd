\name{idMissingLevels}
\alias{idMissingLevels}
\title{
A function to determine the levels of a categorical variable with completely missing outcome data. 
}
\description{
This function determines the levels of a categorical variable with completely missing outcome data. 
}
\usage{
idMissingLevels( yname,  missing.levels.list,  dataset)
}
\arguments{
  \item{yname}{
A character string indicating the outcome variable, i.e., a metabolite.
}
  \item{missing.levels.list}{
A named list indicating whether categorical predictors in argument dataset have levels where argument yname is completely missing. The names of the list objects should be metabolite names, and should include argument yname. 
}
  \item{dataset}{
A data frame containing the variables specified in yname and missing.levels.list. 
}
}
\value{
 Returns a list of the specific levels of the categorical variables included in missing.levels.list that are missing all values for the metabolite specified in yname.

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
missing.levels.check<-lapply( "aspartic_acid", anyMissingLevels, cat.vars=c("pheno", "batch"), dataset=example.data)
names(missing.levels.check)<-"aspartic_acid"
missing.levels.check

#find the specific missing level
idMissingLevels( "aspartic_acid", missing.levels.check, example.data)


}

