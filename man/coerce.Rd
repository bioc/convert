\name{coerce}
\alias{coerce,RGList,marrayRaw-method}
\alias{coerce,marrayRaw,RGList-method}
\alias{coerce,MAList,marrayNorm-method}
\alias{coerce,marrayNorm,MAList-method}
\alias{coerce,RGList,exprSet-method}
\alias{coerce,MAList,exprSet-method}
\alias{coerce,marrayRaw,exprSet-method}
\alias{coerce,marrayNorm,exprSet-method}
\title{Convert Data Objects}

\description{
Convert between limma, marrayClasses and Biobase data objects.
}

\details{
Objects can be converted (coerced) from one class to another using \code{as(object, Class)} where \code{object} is an object to convert and \code{Class} is the name of the class to convert to.
The following conversions are provided:

\tabular{ll}{
From:\tab To:\cr
\code{RGList} \tab \code{marrayRaw}\cr
\code{marrayRaw} \tab \code{RGList}\cr
\code{MAList} \tab \code{marrayNorm}\cr
\code{marrayNorm} \tab \code{MAList}\cr
\code{RGList} \tab \code{exprSet}\cr
\code{MAList} \tab \code{exprSet}\cr
\code{marrayRaw} \tab \code{exprSet}\cr
\code{marrayNorm} \tab \code{exprSet}\cr
}
}

\author{Gordon Smyth}

\seealso{
\code{\link[methods]{as}} in the \code{methods} package.
}

\examples{
x <- new("RGList")
y <- as(x,"marrayRaw")

x <- new("marrayRaw")
x$maGf <- matrix(rnorm(5*2),5,2)
x$maRf <- matrix(rnorm(5*2),5,2)
y <- as(x,"RGList")
y <- as(x,"exprSet")

x <- new("MAList")
y <- as(x,"marrayNorm")

x <- new("marrayNorm")
x$maM <- matrix(rnorm(5*2),5,2)
x$maA <- matrix(rnorm(5*2),5,2)
y <- as(x,"MAList")
y <- as(x,"exprSet")

x <- new("MAList")
x$M <- matrix(rnorm(5*2),5,2)
y <- as(x,"exprSet")
}

\keyword{classes}
\keyword{data}
