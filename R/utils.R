.is_mpra_or_stop <- function(object) {
    if (!is(object, "MPRASet"))
        stop("object is of class '", class(object), "', but needs to be of class 'MPRASet'")
}

.check_assay_names <- function(object, names) {
    nms <- names(assays(object, withDimnames = FALSE))
    if(!all(names %in% nms))
        return(sprintf("object of class '%s' needs to have assay slots with names '%s'",
                       class(object), paste0(names, collapse = ", ")))
    else
        NULL
}

.show.barcodePresence <- function(object) {
    if ("barcode" %in% names(rowData(object)))
        cat("Barcodes present")
    else
        cat("No barcodes present")
}
