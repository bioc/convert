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

##There is intentionally no setAs function with this signature
##Robert Gentleman, Nov 22, 2006
#setAs("marrayRaw", "ExpressionSet", function(from)


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
        if(from@maNormCall != getClass("call")@prototype) y$maNormCall <- from@maNormCall
	y
})

setAs("marrayNorm", "ExpressionSet", function(from){
    ##Robert Gentleman
    ##22 November 2006
    nM <- new("MIAME")
    notes(nM) <- paste(from@maNotes,
                       ":: Converted from marrayNorm object, exprs are log-ratios")
    exprs <- maM(from)
    colnames(exprs) <- NULL
    if (any(duplicated(maLabels(maGnames(from))))) {
        rownames(exprs) <- seq_len(nrow(exprs))
        fData <- data.frame(maGnames=maLabels(maGnames(from)))
    } else {
        rownames(exprs) <- maLabels(maGnames(from))
        fData <- data.frame(row.names=seq_len(nrow(exprs)))
    }

    new("ExpressionSet", exprs=exprs,
        phenoData=new("AnnotatedDataFrame", data=maInfo(maTargets(from))),
        featureData=new("AnnotatedDataFrame", data=fData),
        experimentData=nM)
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

setAs("MAList", "ExpressionSet", function(from)
#  Robert Gentleman
#  22 November 2006
{
    theExperimentData <- new("MIAME")
    notes(theExperimentData) <- list("Converted from MAList object, exprs are M-values")
    theExprs <- as.matrix(from$M)
    thePhenoData <- from$targets
    rownames(thePhenoData) <- colnames(theExprs)
    new("ExpressionSet", exprs = theExprs,
        phenoData = new("AnnotatedDataFrame", data=thePhenoData),
        experimentData = theExperimentData)
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

##there is intentionally no setAs function of this type
##setAs("RGList", "ExpressionSet", function(from)
##Robert Gentleman
##Nov 22, 2006
