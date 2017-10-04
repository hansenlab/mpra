test_agg <- function() {
    data(mpraSetExample)
    data(mpraSetAggExample)
    agg_manual <- getRNA(mpraSetExample, aggregate = TRUE)
    agg_pre <- getRNA(mpraSetAggExample)
    stopifnot(identical(dim(agg_manual), dim(agg_pre)))
    agg_pre <- agg_pre[rownames(agg_manual),]
    colnames(agg_pre) <- gsub("cond_sample_", "", colnames(agg_pre))
    checkEquals(agg_manual, agg_pre)
}

