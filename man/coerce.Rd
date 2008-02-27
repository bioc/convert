\name{coerce}
\alias{convert}
\alias{coerce,RGList,marrayRaw-method}
\alias{coerce,marrayRaw,RGList-method}
\alias{coerce,MAList,marrayNorm-method}
\alias{coerce,marrayNorm,MAList-method}
\alias{coerce,marrayRaw,NChannelSet-method}
\alias{coerce,RGList,NChannelSet-method}
\alias{coerce,MAList,ExpressionSet-method}
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
\code{RGList} \tab \code{NChannelSet}\cr
\code{marrayRaw} \tab \code{NChannelSet}\cr
\code{MAList} \tab \code{ExpressionSet}\cr
\code{marrayNorm} \tab \code{ExpressionSet}\cr
}

\code{RGList} and \code{marrayRaw} are coerced to
\code{NChannelSet}. Channel values are not transformed.

\code{MAList} and \code{marrayNorm} are coerced so that the
\code{ExpressionSet} slot contains log-ratios (M-values) and the
\code{ExpressionSet} object has the same number of columns as the
original object.  In this case, information on the A-values is lost.

There is intentionally no conversion from \code{RGList} or
\code{marrayRaw} to \code{ExpressionSet}, as \code{ExpressionSet} is
intended for expression values, not intensities.
}

\author{Gordon Smyth and others}

\seealso{
\code{\link[methods]{as}} in the \code{methods} package.
}

\examples{
##first set up some fake intensity matrices
testRed <- matrix(rnorm(5*2),5,2, 
  dimnames=list(paste("gene",1:5, sep=""), c("S1", "S2")))
testGreen <- matrix(rnorm(5*2),5,2, 
  dimnames=list(paste("gene",1:5, sep=""), c("S1", "S2")))

##some sample/target info
testTarget <- data.frame(slide=c("S1", "S2"), Cy3=c("T", "C"),
  Cy5=c("C", "T"), row.names=c("S1", "S2"))

maT <- new("marrayInfo", maLabels=c("S1", "S2"),
  maInfo= testTarget)

##now create instances and convert
x <- new("RGList")
x$R <- testRed
x$G <- testGreen
y <- as(x,"marrayRaw")
z <- as(x, "NChannelSet")

x <- new("marrayRaw")
x@maGf <- testGreen
x@maRf <- testRed
x@maTargets = maT
y <- as(x,"RGList")
z <- as(x, "NChannelSet")

x <- new("MAList")
y <- as(x,"marrayNorm")

##we construct a reasonably complete fake, small 
##instance of the marrayNorm class
x <- new("marrayNorm")
x@maM <- testRed
x@maA <- testGreen
maTargets(x) = maT
y <- as(x,"MAList")
y <- as(x, "ExpressionSet")

x <- new("MAList")
x$M <- testRed
x$A <- testGreen
x$targets <- testTarget
y <- as(x,"ExpressionSet")
}
\keyword{classes}
\keyword{data}
