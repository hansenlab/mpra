mpralm <- function(object, design, aggregate = c("mean", "sum", "none"),
                   normalize = TRUE, block = NULL,
                   model_type = c("indep_groups", "corr_groups"),
                   plot = TRUE, ...) {
    .is_mpra_or_stop(object)
    if (nrow(design) != ncol(object)) {
        stop("Rows of design must correspond to the columns of object")
    }
    
    model_type <- match.arg(model_type)
    aggregate <- match.arg(aggregate)
    
    if (model_type=="indep_groups") {
        fit <- .fit_standard(object = object, design = design,
                             aggregate = aggregate, normalize = normalize,
                             plot = plot, ...)
    } else if (model_type=="corr_groups") {
        if (is.null(block)) {
            stop("'block' must be supplied for the corr_groups model type")
        }
        fit <- .fit_corr(object = object, design = design, aggregate = aggregate,
                         normalize = normalize, block = block, plot = plot, ...)
    }
    return(fit)
}

get_precision_weights <- function(logr, design, log_dna, span = 0.4,
                                  plot = TRUE, ...) {
    if (nrow(design) != ncol(logr)) {
        stop("Rows of design must correspond to the columns of logr")
    }

    ## Obtain element-specific residual SDs
    fit <- lmFit(logr, design = design, ...)
    s <- fit$sigma
    x <- rowMeans(log_dna, na.rm = TRUE)
    y <- sqrt(s)
    ## Lowess fitting
    lo <- lowess(x, y, f = span)
    if (plot) {
        plot(x, y, pch = 16, col = alpha("black", 0.25),
             xlab = "Mean(log2(dna+1))", ylab = "sqrt(sd(log-ratio))")
        lines(lo, lwd = 3, col = "red")
    }
    loFun <- approxfun(lo, rule = 2)
    ## Use mean log DNA to get estimated sqrt(SD) to
    ## convert to precision weights
    fittedvals <- log_dna
    w <- 1/loFun(fittedvals)^4
    dim(w) <- dim(fittedvals)
    rownames(w) <- rownames(logr)
    colnames(w) <- colnames(logr)

    return(w)
}

compute_logratio <- function(object, aggregate = c("mean", "sum", "none")) {
    .is_mpra_or_stop(object)

    aggregate <- match.arg(aggregate)

    if (aggregate %in% c("sum", "none")) {
        ## Do aggregation even with option "none" to ensure 
        ## matching ordering of eids in logr and log_dna
        dna <- getDNA(object, aggregate = TRUE)
        rna <- getRNA(object, aggregate = TRUE)
        logr <- log2(rna + 1) - log2(dna + 1)
    } else if (aggregate=="mean") {
        dna <- getDNA(object, aggregate = FALSE)
        rna <- getRNA(object, aggregate = FALSE)
        eid <- getEid(object)
        logr <- log2(rna + 1) - log2(dna + 1)
        
        by_out <- by(logr, eid, colMeans, na.rm = TRUE)
        logr <- do.call("rbind", by_out)
        rownames(logr) <- names(by_out)
    }
    return(logr)
}

normalize_counts <- function(object, block = NULL) {
    .is_mpra_or_stop(object)

    ## Perform total count normalization
    dna <- getDNA(object, aggregate = FALSE)
    rna <- getRNA(object, aggregate = FALSE)

    if (is.null(block)) {
        libsizes_dna <- colSums(dna, na.rm = TRUE)
        libsizes_rna <- colSums(rna, na.rm = TRUE)
    } else {
        libsizes_dna <- tapply(colSums(dna, na.rm = TRUE), block,
                               sum, na.rm = TRUE)
        libsizes_dna <- libsizes_dna[block]
        libsizes_rna <- tapply(colSums(rna, na.rm = TRUE), block,
                               sum, na.rm = TRUE)
        libsizes_rna <- libsizes_rna[block]
    }
    dna_norm <- round(sweep(dna, 2, libsizes_dna, FUN = "/")*10e6)
    rna_norm <- round(sweep(rna, 2, libsizes_rna, FUN = "/")*10e6)
    
    assay(object, "DNA") <- dna_norm
    assay(object, "RNA") <- rna_norm

    return(object)
}

.fit_standard <- function(object, design, aggregate = c("mean", "sum", "none"),
                          normalize = TRUE, return_elist = FALSE,
                          return_weights = FALSE, plot = TRUE, span = 0.4, ...) {
    .is_mpra_or_stop(object)
    if (nrow(design) != ncol(object)) {
        stop("Rows of design must correspond to the columns of object")
    }

    aggregate <- match.arg(aggregate)

    if (normalize) {
        object <- normalize_counts(object)
    }
    logr <- compute_logratio(object, aggregate = aggregate)
    log_dna <- log2(getDNA(object, aggregate = TRUE) + 1)
    
    ## Estimate mean-variance relationship to get precision weights
    w <- get_precision_weights(logr = logr, design = design, log_dna = log_dna,
                               span = span, plot = plot, ...)
    
    elist <- new("EList", list(E = logr, weights = w, design = design))
    
    if (return_weights) {
        return(w)
    }
    if (return_elist) {
        return(elist)
    } 
    fit <- lmFit(elist, design)
    fit <- eBayes(fit)
    fit
}

.fit_corr <- function(object, design, aggregate = c("mean", "sum", "none"),
                     normalize = TRUE, block = NULL, return_elist = FALSE,
                     return_weights = FALSE, plot = TRUE, span = 0.4, ...) {
    .is_mpra_or_stop(object)
    if (nrow(design) != ncol(object)) {
        stop("Rows of design must correspond to the columns of object")
    }

    aggregate <- match.arg(aggregate)

    if (normalize) {
        object <- normalize_counts(object, block)
    }
    logr <- compute_logratio(object, aggregate = aggregate)
    log_dna <- log2(getDNA(object, aggregate = TRUE) + 1)

    ## Estimate mean-variance relationship to get precision weights
    w <- get_precision_weights(logr = logr, design = design, log_dna = log_dna,
                               span = span, plot = plot, ...)

    ## Estimate correlation between element versions that are paired
    corfit <- duplicateCorrelation(logr, design = design,
                                   ndups = 1, block = block)

    elist <- new("EList", list(E = logr, weights = w, design = design))

    if (return_weights) {
        return(w)
    }
    if (return_elist) {
        return(elist)
    } 

    fit <- lmFit(elist, design, block = block, correlation = corfit$consensus)
    fit <- eBayes(fit)
    fit
}
