library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(mpra)

files <- list.files("GSE83894", pattern = "-[DR]NA", full.names = TRUE)
samples <- sapply(files %>% str_split("_"), function(x) {
    str_sub(str_split(x[2], "-")[[1]][1], 3, 3)
}) %>% as.integer
cond <- sapply(files %>% str_split("_"), function(x) {
    str_sub(str_split(x[2], "-")[[1]][1], 1, 2)
})
count_type <- sapply(files %>% str_split("_"), function(x) {
    str_sub(str_split(x[2], "-")[[1]][2], 1, 3)
})
new_colnames <- paste0(cond, "_", samples, "_", count_type)

counts <- lapply(seq_along(files), function(i) {
    read_tsv(files[i], col_names = c("barcode", new_colnames[i], "eid"))
})
counts <- Reduce(function(data1, data2) { full_join(data1, data2) }, counts)
counts <- counts %>%
    mutate(bcid = barcode) %>%
    mutate(eid = str_replace(eid, ":[:digit:]*$", "")) %>%
    gather(key = type, value = count, -barcode, -eid, -bcid) %>%
    separate(type, into = c("condition", "sample", "type"), sep = "_") %>%
    spread(key = type, value = count) %>%
    dplyr::rename(rna = RNA, dna = DNA)

## Barcode level counts
dna_bc <- counts %>%
    select(eid, bcid, condition, sample, dna) %>%
    unite(col = cond_sample, condition, sample) %>%
    spread(key = cond_sample, value = dna)
dna_mat_bc <- as.matrix(dna_bc[,3:ncol(dna_bc)])
rownames(dna_mat_bc) <- dna_bc$bcid
rna_bc <- counts %>%
    select(eid, bcid, condition, sample, rna) %>%
    unite(col = cond_sample, condition, sample) %>%
    spread(key = cond_sample, value = rna)
rna_mat_bc <- as.matrix(rna_bc[,3:ncol(rna_bc)])
rownames(rna_mat_bc) <- rna_bc$bcid

## Check that rows are identical in DNA and RNA
identical(rownames(dna_mat_bc), rownames(rna_mat_bc))
## Check that columns are identical in DNA and RNA
identical(colnames(dna_mat_bc), colnames(rna_mat_bc))

mpraSetExample <- MPRASet(DNA = dna_mat_bc, RNA = rna_mat_bc,
                          eid = dna_bc$eid, barcode = dna_bc$bcid,
                          eseq = NULL
                  )

save(mpraSetExample, file = "../../data/mpraSetExample.rda", compress = "xz")

## Aggregated counts
counts_summ <- counts %>%
    group_by(eid, sample, condition) %>%
    summarize(agg_rna = sum(rna, na.rm = TRUE),
              agg_dna = sum(dna, na.rm = TRUE))
dna <- counts_summ %>%
    select(eid, sample, condition, agg_dna) %>%
    unite(col = cond_sample, condition, sample) %>%
    spread(key = cond_sample, value = agg_dna, sep = "_")
rna <- counts_summ %>%
    select(eid, sample, condition, agg_rna) %>%
    unite(col = cond_sample, condition, sample) %>%
    spread(key = cond_sample, value = agg_rna, sep = "_")

dna_mat <- as.matrix(dna[,2:ncol(dna)])
rownames(dna_mat) <- dna$eid
rna_mat <- as.matrix(rna[,2:ncol(rna)])
rownames(rna_mat) <- rna$eid

## Check that rows are identical in DNA and RNA
identical(rownames(dna_mat), rownames(rna_mat))
## Check that columns are identical in DNA and RNA
identical(colnames(dna_mat), colnames(rna_mat))


mpraSetAggExample <- MPRASet(DNA = dna_mat, RNA = rna_mat,
                          eid = rownames(dna_mat), barcode = NULL, eseq = NULL
                  )

save(mpraSetAggExample, file = "../../data/mpraSetAggExample.rda", compress = "xz")

get_counts_GSE75661 <- function(file) {
    counts <- read_tsv(file)
    counts <- counts %>% 
        gather(key = type, value = count, -Oligo) %>%
        separate(type, into = c("type", "sample"), sep = "_") %>%
        mutate(sample = str_replace(sample, "r", "") %>% as.numeric, bcid = 1) %>%
        spread(key = type, value = count) %>%
        dplyr::rename(eid = Oligo, dna = Plasmid)

    counts_na12878 <- counts %>%
        select(eid, sample, bcid, dna, NA12878) %>%
        dplyr::rename(rna = NA12878)
    counts_na19239 <- counts %>%
        select(eid, sample, bcid, dna, NA19239) %>%
        dplyr::rename(rna = NA19239) %>%
        filter(!is.na(rna))

    ## 7.5k pool wasn't analyzed in HepG2
    has_hepg2 <- any(str_detect(colnames(counts), "Hep"))
    if (any(has_hepg2)) {
        counts_hepg2 <- counts %>%
            select(eid, sample, bcid, dna, HepG2) %>%
            rename(rna = HepG2)
        list(hepg2 = counts_hepg2, na12878 = counts_na12878, na19239 = counts_na19239)
    } else {
        list(na12878 = counts_na12878, na19239 = counts_na19239)
    }
}

file <- "~/Desktop/mpra/GSE75661/GSE75661_7.5k_collapsed_counts.txt"
counts <- get_counts_GSE75661(file)
counts <- counts$na12878

counts_summ <- counts %>%
    mutate(snp_id = str_replace(eid, "_[AB]$", ""),
           allele = str_extract(eid, "[AB]$")) %>%
    select(snp_id, allele, sample, bcid, dna, rna) %>%
    group_by(snp_id, allele, sample) %>%
    summarize(agg_rna = sum(rna, na.rm = TRUE),
              agg_dna = sum(dna, na.rm = TRUE)) %>%
    filter(!is.na(allele))
dna <- counts_summ %>%
    select(snp_id, allele, sample, agg_dna) %>%
    unite(col = allele, allele, sample) %>%
    spread(key = allele, value = agg_dna, sep = "_")
rna <- counts_summ %>%
    select(snp_id, allele, sample, agg_rna) %>%
    unite(col = allele, allele, sample) %>%
    spread(key = allele, value = agg_rna, sep = "_")

cat("Row orders are identical in DNA and RNA:", identical(dna$snp_id, rna$snp_id), "\n")
cat("Columns are identical in DNA and RNA:", identical(colnames(dna), colnames(rna)), "\n")

dna_mat <- as.matrix(dna[,2:ncol(dna)])
rownames(dna_mat) <- dna$snp_id
rna_mat <- as.matrix(rna[,2:ncol(rna)])
rownames(rna_mat) <- rna$snp_id

mpraSetAllelicExample <- MPRASet(DNA = dna_mat, RNA = rna_mat,
                          eid = rownames(dna_mat), barcode = NULL, eseq = NULL
                  )

save(mpraSetAllelicExample, file = "../../data/mpraSetAllelicExample.rda", compress = "xz")
