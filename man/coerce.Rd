\name{coerce}
\alias{coerce,RGList,marrayRaw-method}
\alias{coerce,marrayRaw,RGList-method}
\alias{coerce,MAList,marrayNorm-method}
\alias{coerce,marrayNorm,MAList-method}
\alias{coerce,RGList,exprSet-method}
\alias{coerce,MAList,exprSet-method}
\title{Convert Data Objects}

\description{
Convert between limma, marrayClasses and Biobase data objects.
}

\details{
Objects can be converted (coerced) from one class to another using \code{as(object, Class)} where \code{object} is an object to convert and \code{Class} is the name of the class to convert to.
}

\author{Gordon Smyth}

\seealso{
\code{\link[methods]{as}} in the \code{methods} package.
}

\examples{
x <- new("RGList")
y <- as(x,"marrayRaw")

x <- new("marrayRaw")
y <- as(x,"RGList")

x <- new("MAList")
y <- as(x,"marrayNorm")

x <- new("marrayNorm")
y <- as(x,"MAList")

x <- new("MAList")
x$M <- matrix(rnorm(5*2),5,2)
y <- as(x,"exprSet")
}

\keyword{classes}
\keyword{data}
