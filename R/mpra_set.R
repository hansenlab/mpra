setClass("MPRASet", contains = "SummarizedExperiment")

setValidity("MPRASet", function(object) {
    ## Required information: DNA, RNA, and element ID
    msg <- validMsg(NULL, .check_assay_names(object, c("DNA", "RNA")))
    msg <- validMsg(msg, !is.null(eid(object)))
    ## If barcode is supplied, then all elements must be unique
    bc <- barcode(object)
    if (!is.null(bc)) {
        msg <- validMsg(msg, !any(duplicated(bc)))
    }
    if (is.null(msg)) TRUE else msg
})

MPRASet <- function(DNA = new("matrix"), RNA = new("matrix"),
                    barcode = new("DNAStringSet"), eid = new("character"),
                    eseq = new("DNAStringSet"), ...) {
    assays <- SimpleList(DNA = DNA, RNA = RNA)
	rowData <- DataFrame(barcode = barcode, eid = eid, eseq = eseq)
    new("MPRASet",
        SummarizedExperiment(assays = assays, rowData = rowData, ...)
    )
}

setMethod("show", signature(object = "MPRASet"),
          function(object) {
    callNextMethod()
})

getDNA <- function(object, aggregate = FALSE) {
    .is_mpra_or_stop(object)
    raw <- assay(object, "DNA")
    if (aggregate) {
        eid <- eid(object)
        by_out <- by(raw, eid, colSums)
        agg <- do.call("rbind", by_out)
        rownames(agg) <- names(by_out)
        return(agg)
    } else {
        return(raw)
    }
}

getRNA <- function(object, aggregate = FALSE) {
    .is_mpra_or_stop(object)
    raw <- assay(object, "RNA")
    if (aggregate) {
        eid <- eid(object)
        by_out <- by(raw, eid, colSums)
        agg <- do.call("rbind", by_out)
        rownames(agg) <- names(by_out)
        return(agg)
    } else {
        return(raw)
    }
}

barcode <- function(object) {
    .is_mpra_or_stop(object)
    rowData(object)$barcode
}

eid <- function(object) {
    .is_mpra_or_stop(object)
    rowData(object)$eid
}

eseq <- function(object) {
    .is_mpra_or_stop(object)
    rowData(object)$eseq
}
