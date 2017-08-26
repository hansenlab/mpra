mpralm <- function(object, design, block = NULL, model_type = c("indep_groups", "corr_groups"), ...) {
    .is_mpra_or_stop(object)

    if (model_type=="indep_groups") {
        fit <- fit_standard(object = object, design = design, ...)
    } else if (model_type=="corr_groups") {
        if (is.null(block)) {
            stop("'block' must be supplied for the corr_groups model type")
        }
        fit <- fit_corr(object = object, design = design, block = block, ...)
    }
    return(fit)
}

get_precision_weights <- function(logr, log_dna, span = 0.4, plot = TRUE, ...) {
    ## Obtain element-specific residual SDs
    fit <- lmFit(logr, design = design, ...)
    s <- fit$sigma
    x <- rowMeans(log_dna, na.rm = TRUE)
    y <- sqrt(s)
    ## Lowess fitting
    lo <- lowess(x, y, f = span)
    if (plot) {
        plot(x, y, pch = 16, col = alpha("black", 0.25), xlab = "Mean(log2(dna+1))", ylab = "sqrt(sd(log-ratio))")
        lines(lo, lwd = 3, col = "red")
    }
    loFun <- approxfun(lo, rule = 2)
    ## Use mean log DNA to get estimated sqrt(SD) to convert to precision weights
    fittedvals <- log_dna
    w <- 1/loFun(fittedvals)^4
    dim(w) <- dim(fittedvals)
    rownames(w) <- rownames(logr)
    colnames(w) <- colnames(logr)

    return(w)
}

fit_standard <- function(object, design, return_elist = FALSE, return_weights = FALSE, plot = TRUE, span = 0.4, ...) {
    log_dna <- log2(dna(object) + 1)
    logr <- log2(rna(object) + 1) - log_dna

    ## Estimate mean-variance relationship to get precision weights
    w <- get_precision_weights(logr = logr, log_dna = log_dna,
                               span = span, plot = plot, ...)

    elist <- new("EList", list(E = logr, weights = w, design = design))

    if (return_weights) {
        return(w)
    } else if (return_elist) {
        return(elist)
    } else {
        fit <- lmFit(elist, design)
        fit <- eBayes(fit)
        return(fit)
    }
}

fit_corr <- function(object, design, block = NULL, plot = TRUE, span = 0.4, ...) {
    log_dna <- log2(dna(object) + 1)
    logr <- log2(rna(object) + 1) - log_dna

    ## Estimate mean-variance relationship to get precision weights
    w <- get_precision_weights(logr = logr, log_dna = log_dna,
                               span = span, plot = plot, ...)

    ## Estimate correlation between element versions that are paired
    corfit <- duplicateCorrelation(logr, design = design, ndups = 1, block = block)

    elist <- new("EList", list(E = logr, weights = w, design = design))
    fit <- lmFit(elist, design, block = block, correlation = corfit$consensus)
    fit <- eBayes(fit)
    
    return(fit)
}
