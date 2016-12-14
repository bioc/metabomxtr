\name{allMissingLevels}
\alias{allMissingLevels}
\title{
A function to determine whether metabolite levels are present for at most one level of a categorical variable. 
}
\description{
This function determines metabolite data are present for at most one level of a categorical predictor variable, and thus whether that predictor needs to be removed from the mixture model. 
}
\usage{
allMissingLevels(missing.levels.list, dataset)
}
\arguments{
  \item{missing.levels.list}{
A list output by function idMissingLevels indicating which categorical variables have no corresponding metabolite data for at least one level. 
}
  \item{dataset}{
A data frame containing the variables specified in missing.levels.list. 
}
}
\value{
 Returns a list indicating whether categorical variables in missing.levels.list have metabolite data present for at most one level. 

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
missing.level.ids<-idMissingLevels( "aspartic_acid", missing.levels.check, example.data)

#check to see if those variables have completely missing data 
allMissingLevels(missing.level.ids, example.data)

}

