\name{metabdata}
\alias{metabdata}
\docType{data}
\title{
A sample data set of truncated metabolomics data
}
\description{
This data set contains log transformed metabolite levels and phenotype data of sample 115 of pregnant women. Metabolite levels are contained in columns 11-59, and missing values are indicated by NA. All metabolites included contain at least 5 missing values. Columns 1:10 represent phenotypic data.  
}
\usage{data(metabdata)}
\format{
  A data frame with 115 observations on the following 59 variables.
  \describe{
    \item{\code{sample}}{a numeric vector representing sample number.}
    \item{\code{PHENO}}{a factor representing high and low fasting blood glucose, with levels \code{MomHighFPG} \code{MomLowFPG}.}
    \item{\code{age_ogtt}}{a numeric vector representing the woman's age in years when the oral glucose tolerance test was performed.}
    \item{\code{age_ogtt_mc}}{a numeric vector representing mean centered age_ogtt.}
    \item{\code{ga_ogtt_wks}}{a numeric vector representing gestational age in weeks when the oral glucose tolerance test was performed.}
    \item{\code{ga_ogtt_wks_mc}}{a numeric vector representing mean centered ga_ogtt_wks.}
    \item{\code{FCg}}{a factor representing field center where the data were collected with levels \code{BCD} \code{E} \code{F} \code{P}.}
    \item{\code{parity12}}{a numeric vector with 1 indicative of previous pregnancy and 0 otherwise.}
    \item{\code{storageTimesYears}}{a numeric vector representing the number of years the metabolite sample had been stored prior to assay.}
    \item{\code{storageTimesYears_mc}}{a numeric vector representing mean centered storageTimesYears.}
    \item{\code{ketovaline}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{alpha.ketoglutaric.acid}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{ketoleucine}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{acetoacetate}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{aldohexose}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{beta.alanine}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{methylmalonic.acid}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{creatinine}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{hexuronic.acid}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{ethanolamine}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{glutamine}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{glycolic.acid}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{isoleucine}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{malonic.acid}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{ribose}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{phenylalanine}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{pyruvic.acid}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{hexitol.1}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{lysine}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{disaccharide.1}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{tyrosine}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{leucine}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{hexitol.2}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{disaccharide.2}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{ornithine}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{disaccharide.3}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{beta.tocopherol}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{hexitol.3}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{benzene.1.2.4.triol}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{heptadecane}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{nonadecane}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{tetradecanedioic.acid}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{pentadecanoic.acid}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{undecane}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{methyl.heptadecanoate}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{hydrocarbon}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{deoxyhexose}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{glucose}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{pentose.sugar}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{hexitol.4}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{beta.sitosterol}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{X1.5.anhydroglucitol}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{threonine}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{proline}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{campesterol}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{X6.deoxy.glucose}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{erythronic.acid}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{methyl.myristate}}{a numeric vector representing log2 transformed metabolite abundance.}
    \item{\code{methyl.eicosanoate}}{a numeric vector representing log2 transformed metabolite abundance.}
  }
}

\source{
Scholtens DM, Muehlbauer MJ, Daya NR, Stevens RD, Dyer AR, Lowe LP, Metzger BE, Newgard CB, Bain JR, Lowe WL Jr; HAPO Study Cooperative Research Group. Metabolomics reveals broad-scale metabolic perturbations in hyperglycemic mothers during pregnancy. Diabetes Care. 2014 Jan; 37(1):158-66.
}

\keyword{datasets}
