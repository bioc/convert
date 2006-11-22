\name{coerce}
\alias{convert}
\alias{coerce,RGList,marrayRaw-method}
\alias{coerce,marrayRaw,RGList-method}
\alias{coerce,MAList,marrayNorm-method}
\alias{coerce,marrayNorm,MAList-method}
\alias{coerce,RGList,exprSet-method}
\alias{coerce,RGList,ExpressionSet-method}
\alias{coerce,MAList,exprSet-method}
\alias{coerce,MAList,ExpressionSet-method}
\alias{coerce,marrayRaw,exprSet-method}
\alias{coerce,marrayRaw,ExpressionSet-method}
\alias{coerce,marrayNorm,ExpressionSet-method}
\title{Convert Data Objects}

\description{
Convert between limma, marray and Biobase data objects.
}

\details{
Objects can be converted (coerced) from one class to another using
\code{as(object, Class)} where \code{object} is an object to convert
and \code{Class} is the name of the class to convert to.  The following
conversions are provided:
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
Note that \code{RGList} and \code{marrayRaw} objects are coerced to
\code{exprSet} in such a way that green and red channel intensities are stored in alternate
columns of the \code{exprs} slot expression matrix.  This means that the resulting
\code{exprSet} object has twice as many columns as the original object.
On the other hand, \code{MAList} and \code{marrayNorm} are coerced so
that the resulting \code{exprs} slot contains log-ratios (M-values)
and the \code{exprSet} object has the same number of columns as the
original object.
In this case, information on the A-values is lost.
}

\author{Gordon Smyth}

\seealso{
\code{\link[methods]{as}} in the \code{methods} package.
}

\examples{
x <- new("RGList")
y <- as(x,"marrayRaw")

x <- new("marrayRaw")
x@maGf <- matrix(rnorm(5*2),5,2)
x@maRf <- matrix(rnorm(5*2),5,2)
y <- as(x,"RGList")
y <- as(x,"exprSet")

x <- new("MAList")
y <- as(x,"marrayNorm")

##we construct a reasonably complete fake, small 
##instance of the marrayNorm class
x <- new("marrayNorm")
x@maM <- matrix(rnorm(5*2),5,2)
dimnames(x@maM) = list(paste("gene",1:5, sep=""), c("S1", "S2"))
x@maA <- matrix(rnorm(5*2),5,2)
dimnames(x@maA) = list(paste("gene",1:5, sep=""), c("S1", "S2"))
maT = new("marrayInfo", maLabels=c("S1", "S2"),
  maInfo= data.frame(slide=c("S1", "S2"), Cy3=c("T", "C"),
  Cy5=c("C", "T"), row.names=c("S1", "S2")))
maTargets(x) = maT
y <- as(x,"MAList")
y <- as(x,"exprSet")

x <- new("MAList")
x$M <- matrix(rnorm(5*2),5,2)
dimnames(x$M) = list(paste("gene",1:5, sep=""), c("S1", "S2"))
x$targets = data.frame(slide=c("S1", "S2"), Cy3=c("T", "C"),
  Cy5=c("C", "T"), row.names=c("S1", "S2"))
y <- as(x,"ExpressionSet")
}

\keyword{classes}
\keyword{data}
