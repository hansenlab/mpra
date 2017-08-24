setClass("MPRAset", contains = "SummarizedExperiment")

setValidity("MPRAset", function(object) {
    ## Required assays
    msg <- validMsg(NULL, .check_assay_names(object, c("dna", "rna")))
    ## If barcode is supplied, then all elements must be unique
    bc <- barcode(object)
    if (!is.null(bc)) {
        msg <- validMsg(msg, !any(duplicated(bc)))
    }
    if (is.null(msg)) TRUE else msg
})

MPRAset <- function(dna = new("matrix"), rna = new("matrix"),
                    barcode = new("DNAStringSet"), eid = new("character"),
                    eseq = new("DNAStringSet"), ...) {
    assays <- SimpleList(dna = dna, rna = rna)
	rowData <- DataFrame(barcode = barcode, eid = eid, eseq = eseq)
    new("MPRAset",
        SummarizedExperiment(assays = assays, rowData = rowData, ...)
    )
}

setMethod("show", signature(object = "MPRAset"),
          function(object) {
    callNextMethod()
})

dna <- function(object, aggregate = FALSE) {
    .is_mpra_or_stop(object)
    raw <- assay(object, "dna")
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

rna <- function(object) {
    .is_mpra_or_stop(object)
    raw <- assay(object, "rna")
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
