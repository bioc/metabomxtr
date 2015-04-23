\name{euMetabData}
\alias{euMetabData}
\docType{data}
\title{
A sample data set of truncated metabolomics data.
}
\description{
This data set contains a subset of non-targeted GC/MS data for mother / baby pairs of Northern European ancestry who participated in the Hyperglycemia and Adverse Pregnancy Outcome (HAPO) Metabolomics study. 
}
\usage{data("euMetabData")}
\format{
  A data frame with 120 observations on the following 6 variables.
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
The 120 rows correspond to mother / baby sample triples of analytical interest.  24 of these samples were run in each batch.    Row names ending in "_mf" are for fasting maternal samples.  Row names ending in "_m1" are for maternal samples at 1-hour into an oral glucose tolerance test.  Row names ending in "_bc" are for samples of baby cord blood collected at birth.
}

\keyword{datasets}
