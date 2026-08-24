#annotate CpGs
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(minfi)
library(dplyr)

anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

my_cpgs <- c("cg23956565", "cg00119811", "cg27585641", "cg21221899")

cpg_anno <- anno %>%
  as.data.frame() %>%
  filter(Name %in% my_cpgs) %>%
  dplyr::select(
    Name,
    chr,
    pos,
    UCSC_RefGene_Name,
    UCSC_RefGene_Group,
    Relation_to_Island
  ) %>%
  mutate(
    promoter_based = grepl("TSS200|TSS1500|5'UTR|1stExon", UCSC_RefGene_Group)
  )

cpg_anno
