\name{zdesign-methods}
\docType{methods}
\alias{zdesign}
\alias{zdesign-methods}
\alias{zdesign,data.frame-method}
\alias{zdesign,ExpressionSet-method}
\alias{zdesign,matrix-method}
\title{ ~~ Methods for Function \code{zdesign}  ~~}
\description{
 ~~ Methods for function \code{zdesign}  ~~
}
\section{Methods}{
\describe{

\item{\code{signature(x="data.frame",m="ANY")}}{
The columns of data frame 'x' specified in the input Formula object 'm' are converted to the design matrix for the continuous portion of the mixture model.  
}

\item{\code{signature(x="ExpressionSet",m="ANY")}}{
The columns of the phenoData section of expression set 'x' specified in the input Formula object 'm' are converted to the design matrix for the continuous portion of the mixture model.  
}

\item{\code{signature(x="matrix",m="ANY")}}{
The columns of matrix 'x' specified in the input Formula object 'm' are converted to the design matrix for the continuous portion of the mixture model.  
}
}}
\keyword{methods}

