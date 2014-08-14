\name{yvals-methods}
\docType{methods}
\alias{yvals}
\alias{yvals-methods}
\alias{yvals,data.frame,character-method}
\alias{yvals,ExpressionSet,character-method}
\alias{yvals,matrix,character-method}
\title{ ~~ Methods for Function \code{yvals}  ~~}
\description{
 ~~ Methods for function \code{yvals}  ~~
}
\section{Methods}{
\describe{

\item{\code{signature(y="data.frame",n="character")}}{
The columns of data frame 'y' with names 'n' are converted to a data frame, whose columns are to be used as the response variables in the specified mixture model. 
}

\item{\code{signature(y="ExpressionSet",n="character")}}{
The rows of the assayData section of expression set 'y' with names 'n' are converted to a data frame, whose columns are to be used as the response variables in the specified mixture model.
}

\item{\code{signature(y="matrix",n="character")}}{
The columns of matrix 'y' with names 'n' are converted to a data frame, whose columns are to be used as the response variables in the specified mixture model. 
}
}}
\keyword{methods}

