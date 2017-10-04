test_agg <- function() {
    data(mpraSetExample)
    data(mpraSetAggExample)
    checkEquals(getRNA(mpraSetExample, aggregate = "sum"), getDNA(mpraSetAggExample))
}

