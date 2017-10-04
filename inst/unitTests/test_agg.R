test_agg <- function() {
    data(mpraSetExample)
    data(mpraSetAggExample)
    checkEquals(getRNA(mpraSetExample, aggregate = TRUE), getDNA(mpraSetAggExample))
}

