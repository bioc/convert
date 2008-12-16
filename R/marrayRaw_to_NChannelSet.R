## Martin Morgan
## 27 August, 2007
## Modified from marrayRaw -> NChannelSet
setAs("marrayRaw", "NChannelSet", function(from) {
  
    ## assayData
    elts <- list(R=maRf(from), G=maGf(from))
    if (length(maRb(from))>0)
        elts[["Rb"]] <- maRb(from)
    if (length(maGb(from))>0)
        elts[["Gb"]] <- maGb(from)
    assayData <-
        do.call(assayDataNew,
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
