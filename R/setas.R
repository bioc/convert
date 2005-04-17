#	SETAS.R

setAs("marrayRaw", "RGList", function(from, to) {
#	Gordon Smyth
#	20 Dec 2003

	y <- new(to)
	ifposlen <- function(x) if(length(x)) return(x) else return(NULL)
	y$R <- ifposlen(from@maRf)
	y$G <- ifposlen(from@maGf)
	y$Rb <- ifposlen(from@maRb)
	y$Gb <- ifposlen(from@maGb)
	y$weights <- ifposlen(from@maW)
	y$printer$ngrid.r <- ifposlen(from@maLayout@maNgr)
	y$printer$ngrid.c <- ifposlen(from@maLayout@maNgc)
	y$printer$nspot.r <- ifposlen(from@maLayout@maNsr)
	y$printer$nspot.c <- ifposlen(from@maLayout@maNsc)
	y$printer$notes <- ifposlen(from@maLayout@maNotes)
	y$genes <- ifposlen(from@maGnames@maInfo)
	y$genes$Labels <- ifposlen(from@maGnames@maLabels)
	attr(y$genes,"notes") <- ifposlen(from@maGnames@maNotes)
	y$genes$Sub <- ifposlen(from@maLayout@maSub)
	y$genes$Plate <- ifposlen(from@maLayout@maPlate)
	y$genes$Controls <- ifposlen(from@maLayout@maControls)
	y$targets <- ifposlen(from@maTargets@maInfo)
	y$targets$Labels <- ifposlen(from@maTargets@maLabels)
	y$notes <- ifposlen(from@maNotes)
	y
})

setAs("RGList", "marrayRaw", function(from, to) {
#	Gordon Smyth
#	20 Dec 2003
#      Last modified by Jean Yang, June 30, 2004. Add conversation of control status
	y <- new(to)
	if(!is.null(from$G)) y@maGf <- from$G
	if(!is.null(from$R)) y@maRf <- from$R
	if(!is.null(from$Gb)) y@maGb <- from$Gb
	if(!is.null(from$Rb)) y@maRb <- from$Rb
	if(!is.null(from$weights)) y@maW <- from$weights
	if(!is.null(from$printer$ngrid.r)) y@maLayout@maNgr <- from$printer$ngrid.r
	if(!is.null(from$printer$ngrid.c)) y@maLayout@maNgc <- from$printer$ngrid.c
	if(!is.null(from$printer$nspot.r)) y@maLayout@maNsr <- from$printer$nspot.r
	if(!is.null(from$printer$nspot.c)) y@maLayout@maNsc <- from$printer$nspot.c
	if(!is.null(from$genes)) y@maGnames@maInfo <- from$genes
        if(!is.null(from$genes$Status)) y@maLayout@maControls <- as.factor(from$genes$Status)
	if(!is.null(from$targets)) y@maTargets@maInfo <- from$targets
	if(!is.null(from$notes)) y@maNotes <- "Converted from RGList object"
	y
})

setAs("marrayNorm", "MAList", function(from, to)
#	Gordon Smyth
#	20 Dec 2003. Last modified 18 April 2005.
{
	y <- new(to)
	ifposlen <- function(x) if(length(x)) return(x) else return(NULL)
	y$A <- ifposlen(from@maA)
	y$M <- ifposlen(from@maM)
	y$weights <- ifposlen(from@maW)
	y$printer$ngrid.r <- ifposlen(from@maLayout@maNgr)
	y$printer$ngrid.c <- ifposlen(from@maLayout@maNgc)
	y$printer$nspot.r <- ifposlen(from@maLayout@maNsr)
	y$printer$nspot.c <- ifposlen(from@maLayout@maNsc)
	y$printer$notes <- ifposlen(from@maLayout@maNotes)
	y$genes <- ifposlen(from@maGnames@maInfo)
	y$genes$Labels <- ifposlen(from@maGnames@maLabels)
	attr(y$genes,"notes") <- ifposlen(from@maGnames@maNotes)
	y$genes$Sub <- ifposlen(from@maLayout@maSub)
	y$genes$Plate <- ifposlen(from@maLayout@maPlate)
	y$genes$Controls <- ifposlen(from@maLayout@maControls)
	y$targets <- ifposlen(from@maTargets@maInfo)
	y$targets$Labels <- ifposlen(from@maTargets@maLabels)
	y$notes <- ifposlen(from@maNotes)
        ##modified with Gordon's suggestions
        ##y$maNormCall <- from@maNormCall
	##if(as.character(y$maNormCall)=="<undef>") y$maNormCall <- NULL
        if(from@maNormCalll != getClass("call")@prototype) y$maNormCall <- from@maNormCall
	y
})

setAs("MAList", "marrayNorm", function(from, to)
#	Gordon Smyth
#	20 Dec 2003. Last modified 7 March 2004.
{
	y <- new(to)
	if(!is.null(from$A)) y@maA <- from$A
	if(!is.null(from$M)) y@maM <- from$M
	if(!is.null(from$weights)) y@maW <- from$weights
	if(!is.null(from$printer$ngrid.r)) y@maLayout@maNgr <- from$printer$ngrid.r
	if(!is.null(from$printer$ngrid.c)) y@maLayout@maNgc <- from$printer$ngrid.c
	if(!is.null(from$printer$nspot.r)) y@maLayout@maNsr <- from$printer$nspot.r
	if(!is.null(from$printer$nspot.c)) y@maLayout@maNsc <- from$printer$nspot.c
	if(!is.null(from$genes)) y@maGnames@maInfo <- from$genes
	if(!is.null(from$targets)) y@maTargets@maInfo <- from$targets
	y@maNotes <- "Converted from MAList object"
	y
})

setAs("RGList", "exprSet", function(from, to)
#	Gordon Smyth
#	7 March 2004.  Last modified 28 June 2004.
{
	y <- new(to)
#	Assemble green and red intensities into alternate columns
	d <- dim(from)
	exprs <- array(0,c(d,2))
	exprs[,,1] <- from$G
	exprs[,,2] <- from$R
	exprs <- aperm(exprs,c(1,3,2))
	dim(exprs) <- c(d[1],2*d[2])
	y@exprs <- exprs
	if(!is.null(from$targets)) y@phenoData@pData <- targetsA2C(from$targets)
	y@notes <- "Converted from RGList object, exprs are green/red intensites in odd/even columns"
    y
})

setAs("MAList", "exprSet", function(from, to)
#	Gordon Smyth
#	7 March 2004. Last modified 16 March 2004.
{
	y <- new(to)
	if(!is.null(from$M)) y@exprs <- as.matrix(from$M)
	if(!is.null(from$targets)) y@phenoData@pData <- from$targets
	y@notes <- "Converted from MAList object, exprs are M-values"
    y
})

setAs("marrayRaw", "exprSet", function(from)
## Assemble green and red intensities into alternate columns
## Jean Yang
## 15 March 2004.
#	Modified by Gordon Smyth 28 June 2004.
{
  eset<-new("exprSet")
  d <- dim(from@maGf)
  exprs <- array(0,c(d,2))
  exprs[,,1] <- maLG(from)
  exprs[,,2] <- maLR(from)
  exprs <- aperm(exprs,c(1,3,2))
  dim(exprs) <- c(d[1],2*d[2])
  eset@exprs <- exprs
  if(length(from@maTargets@maInfo)) {
    targets <- targetsA2C(maInfo(maTargets(from)),grep=TRUE)
    pdata <- new("phenoData", pData=targets, varLabels=as.list(names(targets)))
    eset@phenoData <- pdata
  }
  eset@notes <- paste(from@maNotes, ":: Converted from marrayRaw object, exprs are green/red log-intensites in odd/even columns")
  eset
})

setAs("marrayNorm", "exprSet", function(from)
## Jean Yang
## 15 March 2004
{
  eset <- new("exprSet", exprs=maM(from))
  targets <- maInfo(maTargets(from))
  eset@phenoData  <- new("phenoData", pData=targets, varLabels=as.list(names(targets)))
  eset@notes <- paste(from@maNotes, ":: Converted from marrayNorm object, exprs are log-ratios")
  eset
})

