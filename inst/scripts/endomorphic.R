library(mpra)
data(mpraSetExample)
m <- mpraSetExample

colData(m)$condition <- factor(c("MT","MT","MT","WT","WT","WT"),
                                            levels=c("WT","MT"))
random_ids <- apply(
  matrix(sample(letters, 10 * length(unique(rowData(m)$eid)),
                replace = TRUE), ncol=10),
  1, paste0, collapse="")

rowData(m)$eid <- as.numeric(factor(rowData(m)$eid))
rowData(m)$eid <- random_ids[ rowData(m)$eid ]
rowData(m)$score <- rnorm(nrow(m))
rowData(m)

design <- model.matrix(~condition, colData(m))

fit <- mpralm(
  object = m,
  design = design,
  aggregate = "sum",
  normalize = TRUE,
  model_type = "indep_groups",
  plot = FALSE
)

class(fit)

tab <- topTable(fit, coef = 2, number = Inf)
head(tab)
