\name{addOutlierInfo}
\alias{addOutlierInfo}
\title{
A function to determine whether specific metabolite observations are outliers. 
}
\description{
This function adds a column to a dataframe of metabolite abundance values indicating whether a observations for a specific metabolite are outliers. 
}
\usage{
addOutlierInfo(metab.col, in.data, outlier.sd.thresh=2)
}
\arguments{
  \item{metab.col}{
A character value indicating the name of the variable in the input data set that corresponds to metabolite abundance.
}
  \item{in.data}{
A data frame containing the variables specified in metab.col.
}
  \item{outlier.sd.thresh}{
The number of standard deviations from the mean a point must be to be considered an outlier. This argument defaults to 2. 
}
}
\details{
This function adds a column to a dataframe of metabolite abundance values indicating whether a observations for a specific metabolite are outliers. 
}
\value{
Returns a dataframe. 
}
\references{
Nodzenski M, Muehlbauer MJ, Bain JR, Reisetter AC, Lowe WL Jr, Scholtens DM. Metabomxtr: an R package for mixture-model analysis of non-targeted metabolomics data. Bioinformatics. 2014 Nov 15;30(22):3287-8.
}
\author{
Michael Nodzenski
}
\examples{
data(euMetabCData)
outlier.data<-addOutlierInfo("pyruvic_acid",  euMetabCData, outlier.sd.thresh=2)
              
}

