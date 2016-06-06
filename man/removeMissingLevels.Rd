\name{removeMissingLevels}
\alias{removeMissingLevels}
\title{
A function to remove levels of categorical variables with completely missing outcome data from data frames or matrices. 
}
\description{
This function removes levels of categorical variables with completely missing outcome data from data frames or matrices. 
}
\usage{
removeMissingLevels( yname,  missing.levels.list,  dataset)
}
\arguments{
  \item{yname}{
A character string indicating the outcome variable, i.e., a metabolite.
}
  \item{missing.levels.list}{
A list, output by function idMissingLevels, containing the specific levels of categorical variables that are missing all outcome data. 
}
  \item{dataset}{
A data frame containing the variables specified in yname and missing.levels.list. 
}
}
\value{
 Returns a data frame with the levels indicated in missing.levels.list removed. 

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
missing.levels<-idMissingLevels( "aspartic_acid", missing.levels.check, example.data)

#remove the missing level 
cleaned.data<-removeMissingLevels("aspartic_acid", missing.levels, example.data)


}

