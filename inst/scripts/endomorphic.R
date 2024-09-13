# library(mpra)
devtools::load_all()
data(mpraSetExample)
m <- mpraSetExample[1:10000,] # for speed
colData(m)$condition <- factor(c("MT","MT","MT","WT","WT","WT"),
                               levels=c("WT","MT"))

# no aggregation
rowData(m)$eid <- paste0("e",seq_len(nrow(m)))
rowData(m)$score <- rnorm(nrow(m))
rowData(m)

design <- model.matrix(~condition, colData(m))

fit <- mpralm(object = m, design = design, aggregate = "none",
              normalize = TRUE, model_type = "indep_groups",
              plot = FALSE, endomorphic = TRUE, coef = 2)

# previously, an MArrayLM object
class(fit) # MPRASet object

tab <- topTable(attr(fit, "MArrayLM"), coef = 2, number = Inf)
head(tab)

all.equal(tab[rowData(fit)$eid,"logFC"], rowData(fit)$logFC)

# now try it with aggregation
data(mpraSetExample)
m <- mpraSetExample[1:10000,] # for speed
colData(m)$condition <- factor(c("MT","MT","MT","WT","WT","WT"),
                               levels=c("WT","MT"))
rowData(m)$eid <- paste("e",as.numeric(factor(rowData(m)$eid)))

fit <- mpralm(object = m, design = design, aggregate = "sum",
              normalize = TRUE, model_type = "indep_groups",
              plot = FALSE, endomorphic = TRUE, coef = 2)

class(fit) # MPRASet object

tab <- topTable(attr(fit, "MArrayLM"), coef = 2, number = Inf)
head(tab)

all.equal(tab[rowData(fit)$eid,"logFC"], rowData(fit)$logFC)
