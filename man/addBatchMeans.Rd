\name{addBatchMeans}
\alias{addBatchMeans}
\title{
A function to append batch specific mean metabolite abundances to data frames containing raw values. 
}
\description{
This function appends batch specific mean metabolite abundances to data frames containing raw values. 
}
\usage{
addBatchMeans(metab.col, batch.col, in.data, sample.type.name, analytical.type.name, color.col=NULL)
}
\arguments{
  \item{metab.col}{
A character value indicating the name of the variable in the input data set that corresponds to metabolite abundance.
}
  \item{batch.col}{
A character value indicating the name of the variable in the input data sets that corresponds to batch. If not specified, this argument defaults to "Batch". 
}
  \item{in.data}{
A data frame containing the variables specified in metab.col, batch.col, and color.col.
}
  \item{sample.type.name}{
A character value indicating the type of sample: "Quality Control" or "Experimental".
}
  \item{analytical.type.name}{
A character value indicating whether the data have been normalized. Either "Before Normalization" or "After Normalization".
}
  \item{color.col}{
A character value indicating the name of the variable in the input data set corresponding to a grouping variable (for instance mom vs. baby samples).
}
}
\details{
This function adds rows reporting batch specific mean metabolite abundances to data frames of raw values. This is not meant to be a stand alone function but rather work as a part of function metabplot. 
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
euMetabCData$outlier<-FALSE
cdata.with.batch.means<-addBatchMeans("pyruvic_acid", "batch", euMetabCData, sample.type.name="Quality Control", analytical.type.name="Before Normalization", color.col="pheno")
              
}

