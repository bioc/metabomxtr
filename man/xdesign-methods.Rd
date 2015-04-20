\name{xdesign-methods}
\docType{methods}
\alias{xdesign}
\alias{xdesign-methods}
\alias{xdesign,data.frame-method}
\alias{xdesign,ExpressionSet-method}
\alias{xdesign,matrix-method}
\title{ ~~ Methods for Function \code{xdesign}  ~~}
\description{
 ~~ Methods for function \code{xdesign}  ~~
}
\section{Methods}{
\describe{

\item{\code{signature(x="data.frame",m="ANY")}}{
The columns of data frame 'x' specified in the input Formula object 'm' are converted to the design matrix for the discrete portion of the mixture model.  
}

\item{\code{signature(x="ExpressionSet",m="ANY")}}{
The columns of the phenoData section of expression set 'x' specified in the input Formula object 'm' are converted to the design matrix for the discrete portion of the mixture model.  
}

\item{\code{signature(x="matrix",m="ANY")}}{
The columns of matrix 'x' specified in the input Formula object 'm' are converted to the design matrix for the discrete portion of the mixture model.  
}
}}
\keyword{methods}

