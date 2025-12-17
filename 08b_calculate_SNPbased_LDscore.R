#!/usr/bin/env Rscript
# R script: convert_bed_to_bim.R

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Ensure two arguments are provided
if (length(args) < 2) {
  stop("Usage: Rscript convert_bed_to_bim.R <subdir> <filename>")
}

# Define arguments
subdir <- args[1]
filename <- args[2]

base_dir <- "/home/projects/hearing_loss/clsaARHL_SA/str"

lds_seg = read.table(paste0(base_dir, "/heritability/", subdir, "/input/", filename),header=T,colClasses=c("character",rep("numeric",8)))
quartiles=summary(lds_seg$ldscore_SNP)

lb1 = which(lds_seg$ldscore_SNP <= quartiles[2])
lb2 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP <= quartiles[3])
lb3 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP <= quartiles[5])
lb4 = which(lds_seg$ldscore_SNP > quartiles[5])

lb1_snp = lds_seg$SNP[lb1]
lb2_snp = lds_seg$SNP[lb2]
lb3_snp = lds_seg$SNP[lb3]
lb4_snp = lds_seg$SNP[lb4]

write.table(lb1_snp, paste0(base_dir, "/heritability/", subdir, "/input/snp_group1.txt"),
            row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(lb2_snp, paste0(base_dir, "/heritability/", subdir, "/input/snp_group2.txt"),
            row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(lb3_snp, paste0(base_dir, "/heritability/", subdir, "/input/snp_group3.txt"),
            row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(lb4_snp, paste0(base_dir, "/heritability/", subdir, "/input/snp_group4.txt"),
            row.names = FALSE, quote = FALSE, col.names = FALSE)
