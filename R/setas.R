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
        if(from@maNormCall != getClass("call")@prototype) y$maNormCall <- from@maNormCall
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
        pD <- targetsA2C(from$targets)
        varL = as.list(names(pD))
        names(varL) = names(pD)
	if(!is.null(from$targets)) 
              phenoData(y) <- new("phenoData", pData = pD,
                  varLabels = varL)
	y@notes <- "Converted from RGList object, exprs are green/red intensites in odd/even columns"
    y
})

##there is intentionally no setAs function of this type
##setAs("RGList", "ExpressionSet", function(from)
##Robert Gentleman
##Nov 22, 2006


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

setAs("MAList", "ExpressionSet", function(from)
#  Robert Gentleman
#  22 November 2006
{
    nM = new("MIAME")
    notes(nM) = list("Converted from MAList object, exprs are M-values")
    new("ExpressionSet", exprs = as.matrix(from$M),
        phenoData = new("AnnotatedDataFrame", data=from$targets),
        experimentData = nM)
})

setAs("marrayRaw", "NChannelSet", function(from)
## Martin Morgan
## 27 August, 2007
## Modified from marrayRaw -> exprSet
{
    ## assayData
    elts <- list(R=maRf(from), G=maGf(from))
    if (length(maRb(from))>0)
        elts[["Rb"]] <- maRb(from)
    if (length(maGb(from))>0)
        elts[["Gb"]] <- maGb(from)
    assayData <-
        do.call("assayDataNew",
                c(storage.mode="lockedEnvironment", elts))
    ## phenoData, featureData
    pData <- 
        if (length(maInfo(maTargets(from)))>0) {
            data=maInfo(maTargets(from))
        } else {
            data=data.frame(rep(0, ncol(from)))[,FALSE]
        }
    phenoData <- new("AnnotatedDataFrame", data=pData)
    fData <- 
        if (length(maInfo(maGnames(from)))>0) {
            maInfo(maGnames(from))
        } else {
            data.frame(rep(0, nrow(from)))[,FALSE]
        }
    
    if (!is.null(rownames(assayData[["R"]])))
        row.names(fData) <- rownames(assayData[["R"]])
    featureData <- new("AnnotatedDataFrame", data=fData)
    ## experimentData
    experimentData <- new("MIAME",
                          other=list("::Converted from marrayRaw object"))
    ## NChannelSet
    obj <- new("NChannelSet",
               assayData=assayData,
               phenoData=phenoData,
               featureData=featureData,
               experimentData=experimentData)
    ## adjustments
    if (!is.null(obj[["Names"]])) {
        phenoData(obj)[["FileName",
                        labelDescription="Source file name"]] <-
                            sampleNames(obj)[[1]]
        sampleNames(obj) <- obj[["Names"]]
    }
    lbls <- maLabels(maGnames(from))
    if (length(lbls)==nrow(from)) {
        if (any(duplicated(lbls)))
            featureData(obj)[["maLabels",
                              labelDescription="marrayRaw gene names"]] <-
                                  lbls
        else
            featureNames(obj) <- lbls
    }
    obj
})

setAs("RGList", "NChannelSet", function(from)
## Martin Morgan
## 27 August, 2007
{
    ## assayData
    assayData <- with(from, {
        if (!exists("other", inherits=FALSE)) {
            elts <- list(R=R, G=G)
        } else {
            if (is.null(names(other)) ||
                !all(sapply(names(other), nzchar)))
                stop(paste("RGList 'other' elements must be named, found '",
                           paste(names(other), collapse="', '"),
                           "'", sep=""))
            bad <- names(other) %in% c("R", "G", "Rb", "Gb")
            if (any(bad))
                stop(paste("RGList 'other' elements contain reserved names '",
                           paste(names(other)[bad],
                                 collapse="', '"),
                           "'", sep=""))
                
            elts <- c(R=R, G=G, other)
        }
        if (exists("Rb", inherits=FALSE))
            elts[["Rb"]] <- Rb
        if (exists("Gb", inherits=FALSE))
            elts[["Gb"]] <- Gb
        do.call("assayDataNew",
                c(storage.mode="lockedEnvironment", elts))
    })
    ## phenoData
    phenoData <- 
        if (!is.null(from$target)) {
            new("AnnotatedDataFrame",
                data=data.frame(FileName=from$target),
                varMetadata=data.frame(
                  labelDescription="Source file name"))
        } else {
            new("AnnotatedDataFrame",
                data=data.frame(rep(0, ncol(from)))[,FALSE])
    }
    ## featureData
    fData <-
        if (!is.null(from$genes))
            from$genes
        else
            data.frame(x=rep(0,nrow(from)))[,FALSE]
    if (!is.null(rownames(assayData[["R"]])))
        row.names(fData) <- rownames(assayData[["R"]])
    featureData <- new("AnnotatedDataFrame", data = fData)
    if (!is.null(from$weights))
        if ("weights" %in% names(df))
            warning("RGList 'genes' contains column 'weights'; 'wt.fun' weights discarded")
        else
            featureData[["weights",
                         labelDescription="calculated, from RGList"]] <-
                             from$weights
    ## experimentData
    other <- 
        if (!is.null(from$source))
            list(source=from$source)
        else
            list()
    experimentData <- new("MIAME",
                          other=c("converted from marrayRaw", other))
    new("NChannelSet",
        assayData = assayData,
        featureData = featureData,
        phenoData = phenoData,
        experimentData=experimentData)
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

##There is intentionally no setAs function with this signature
##Robert Gentleman, Nov 22, 2006
#setAs("marrayRaw", "ExpressionSet", function(from)

setAs("marrayNorm", "exprSet", function(from)
## Jean Yang
## 15 March 2004
{
  eset <- new("exprSet")
  eset@exprs <- maM(from)
  targets <- maInfo(maTargets(from))
  eset@phenoData  <- new("phenoData", pData=targets, varLabels=as.list(names(targets)))
  eset@notes <- paste(from@maNotes, ":: Converted from marrayNorm object, exprs are log-ratios")
  eset
})

setAs("marrayNorm", "ExpressionSet", function(from){
	##Robert Gentleman
	##22 November 2006
	nM = new("MIAME")
	notes(nM) = paste(from@maNotes, ":: Converted from marrayNorm object, exprs are log-ratios")
    exprs <- maM(from)
    colnames(exprs) <- NULL
    rownames(exprs) <- maLabels(maGnames(from))
	new("ExpressionSet", exprs=exprs,
		phenoData=new("AnnotatedDataFrame", data=maInfo(maTargets(from))),
		experimentData=nM)
})
