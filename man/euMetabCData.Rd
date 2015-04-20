\name{euMetabCData}
\alias{euMetabCData}
\docType{data}
\title{
A sample data set of truncated metabolomics data
}
\description{
This data set contains a subset of non-targeted GC/MS data for mother / baby pairs of Northern European ancestry who participated in the Hyperglycemia and Adverse Pregnancy Outcome (HAPO) Metabolomics study. 
}
\usage{data("euMetabCData")}
\format{
  A data frame with 30 observations on the following 6 variables.
  \describe{
    \item{\code{batch}}{A factor variable for the batch in which the sample was processed with levels \code{1} \code{2} \code{3} \code{4} \code{5}.}
    \item{\code{pheno}}{A factor variable indicating whether the serum sample was from a mother or baby with levels \code{BABY} \code{MOM}.}
    \item{\code{betahydroxybutyrate}}{A numeric vector of non-normalized GC/MS log2 transformed peak areas for this metabolite.}
    \item{\code{pyruvic_acid}}{A numeric vector of non-normalized GC/MS log2 transformed peak areas for this metabolite.}
    \item{\code{malonic_acid}}{A numeric vector of non-normalized GC/MS log2 transformed peak areas for this metabolite.}
    \item{\code{aspartic_acid}}{A numeric vector of non-normalized GC/MS log2 transformed peak areas for this metabolite.}
  }
}
\details{
The 30 rows correspond to 3 mom and 3 baby control samples run at the beginning, middle and end of each batch.  All control samples were drawn from an identical pool. 
}

\keyword{datasets}
