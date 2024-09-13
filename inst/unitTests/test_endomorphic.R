test_endo <- function() {
    data(mpraSetExample)
    m <- mpraSetExample[1:10000,]
    colData(m)$condition <- factor(c("MT","MT","MT","WT","WT","WT"),
                                   levels=c("WT","MT"))
    rowData(m)$eid <- paste("e",as.numeric(factor(rowData(m)$eid)))
    design <- model.matrix(~condition, colData(m))
    fit <- mpralm(object = m, design = design, aggregate = "sum",
                  normalize = TRUE, model_type = "indep_groups",
                  plot = FALSE, endomorphic = TRUE, coef = 2)
    checkTrue(is(fit, "MPRASet"))
    tab <- topTable(attr(fit, "MArrayLM"), coef = 2, number = Inf)
    checkEquals(tab[rowData(fit)$eid,"logFC"], rowData(fit)$logFC)
}
