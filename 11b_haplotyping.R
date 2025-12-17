#hapltype analysis of CHR 5

setwd("/home/projects/hearing_loss/clsaARHL_SA/str/haplotype/")
#install.packages("GHap")
library(GHap)

# Path to VCF
vcf_file <- "/home/projects/hearing_loss/clsaARHL_SA/str/haplotype/filtered_imputed_STR_SNPs_chr5_subset_biallelic"

#prepare sample file (this should have population column and ID)
samples <- readLines("/home/projects/hearing_loss/clsaARHL_SA/str/haplotype/temp_samples.txt")
ghap_samples <- data.frame(POP = rep("all", length(samples)), ID = samples)
write.table(ghap_samples, 
            "/home/projects/hearing_loss/clsaARHL_SA/str/haplotype/filtered_imputed_STR_SNPs_chr5_subset_biallelic.sample",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")


# Convert to GHap format
ghap.vcf2phase(
  input.files = vcf_file,
  out.file = "chr5_window"  # prefix for output files
)

#remove duplicates
markers <- read.table("chr5_window.markers", header = FALSE, stringsAsFactors = FALSE)
sum(duplicated(markers$V2))  # should match 4348
duplicated_ids <- markers$V2[duplicated(markers$V2)]
head(duplicated_ids)

#Make marker IDs unique
markers$V2 <- make.unique(markers$V2)
write.table(markers, "chr5_window_unique.markers", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

#compress the data to get the phaseb binary format
ghap.compress(samples.file = "chr5_window.samples", markers.file = "chr5_window_unique.markers", phase.file = "chr5_window.phase", out.file = "chr5_window")


#load the phased data
phase <- ghap.loadphase(samples.file = "chr5_window.samples", markers.file = "chr5_window_unique.markers", phaseb.file = "chr5_window.phaseb")

#subset markers with MAF > 0.01??
freq <- ghap.freq(phase, type = 'maf')
mkr <- names(freq)[which(freq > 0.01)]
# phase <- ghap.subset(object = phase, ids = unique(phase$id), variants = mkr) # ids are usually duplicated by GHap to store halotypes
# 
# #define haplotype blocks
# blocks.mkr <- ghap.blockgen(phase, windowsize = 15, slide = 2, unit = "marker")
# 
# #generate the haplotypes (keep haplotypes with frequency >0.01)
# ghap.haplotyping(object = phase, blocks = blocks.mkr, outfile = "chr5_window", binary = TRUE, freq = c(0.01, 1), ncores = 5)
# 
# # Load haplotype genotypes using prefix 
# haplo <- ghap.loadhaplo("chr5_window")
# 
# # Convert to plink 
# ghap.hap2plink(haplo, outfile = "chr5_window")
# 
# #update the .fam file to match plink format
# haplo_fam <- read.table("chr5_window.fam", header = FALSE)
# haplo_fam <- haplo_fam[, 2:6]
# haplo_fam.1 <- separate(haplo_fam, col = V2, into = c("V1", "V2"), sep = "_")
# haplo_fam.1$V1 <- as.numeric(haplo_fam.1$V1)
# haplo_fam.1$V2 <- as.numeric(haplo_fam.1$V1)
# 
# #save the file
# 
# write.table(haplo_fam.1, 
#             "/home/projects/hearing_loss/clsaARHL_SA/str/haplotype/chr5_window.fam",
#             quote = FALSE, row.names = FALSE, col.names = FALSE)
# 
# #read the output file
# assoc_file1 <- read.table("chr5_window_met_gwas_sorted_noNA_add", header = TRUE)
# 
# #find the haploblock that contains my top STR Ens:chr5:73782133 and the top GWAS variant chr5:73780686
# str_pos <- 73782133 #lowest p B3171 B3172
# str_pos2 <-73780643 # B3171 B3172
# str_pos3 <-73778077 #highest pip B3171 B3172
# snp_pos <- 73780686 #snp with lowest p B3171 B3172
# snp_pos2 <- 73776529 #second sig missense
# 
# str_block <- blocks.mkr[blocks.mkr$BP1 <= str_pos4 & blocks.mkr$BP2 >= str_pos4, ]
# snp_block <- blocks.mkr[blocks.mkr$BP1 <= snp_pos & blocks.mkr$BP2 >= snp_pos, ]
# 
# snp_block2 <- blocks.mkr[blocks.mkr$BP1 <= snp_pos2 & blocks.mkr$BP2 >= snp_pos2, ]

#########################
##try a haploblock with the two missense and the str with the highest pip
#load the phased data
phase <- ghap.loadphase(samples.file = "chr5_window.samples", markers.file = "chr5_window_unique.markers", phaseb.file = "chr5_window.phaseb")

#EnsTR:chr5:73778077 STR can have multiple copies depending on the number of alt alleles after being converted to bi-allelic format
# Suppose the STR of interest
target_str <- "EnsTR:chr5:73778077"

# Find all markers that start with this STR (ignores suffixes)
idx <- grep(paste0("^", target_str), phase$marker)

# Show all matching markers
str_mkr <-  phase$marker[idx]

#subset markers
sub_mkr <- c("5:73780686:C:A", "5:73776529:T:C", str_mkr)

phase_subset <- ghap.subset(phase, ids = unique(phase$id), variants = sub_mkr)

#define haplotype blocks
blocks.mkr_subset <- ghap.blockgen(phase_subset, windowsize = 8, slide = 2, unit = "marker") #we have only 8 markers
 
#generate the haplotypes 
ghap.haplotyping(object = phase_subset, blocks = blocks.mkr_subset, outfile = "top_var", binary = TRUE, freq = c(0.01, 1))

# Load haplotype genotypes using prefix 
haplo <- ghap.loadhaplo("top_var")

# Convert to plink 
ghap.hap2plink(haplo, outfile = "top_var")

#read the output file
assoc_file <- read.table("top_var_met_gwas_sorted_noNA_add", header = TRUE)

#scale the effect size
assoc_file$scaled_beta <- scale(assoc_file$BETA)

#############
###
#Targeted haplotype analysis
#read phased data
# phase <- ghap.loadphase(samples.file = "chr5_window.samples", markers.file = "chr5_window_unique.markers", phaseb.file = "chr5_window.phaseb")
# 
# #try removing maf <0.01 and low frequency haplotypes
# freq <- ghap.freq(phase, type = 'maf')
# mkr <- names(freq)[which(freq > 0.01)]
# phase <- ghap.subset(object = phase, ids = unique(phase$id), variants = mkr)
# 
# 
# # Suppose the STR of interest
# target_str <- "EnsTR:chr5:73778077"
# 
# # Find all markers that start with this STR (ignores suffixes)
# idx <- grep(paste0("^", target_str), phase$marker)
# 
# # Show all matching markers
# str_mkr <-  phase$marker[idx]
# 
# # Loop from 0 to 5
# for (i in 0:5) {
# 
#   # Define marker set (i+1)
#   this_mkr <- c("5:73780686:C:A", "5:73776529:T:C", str_mkr[i+1])
# 
#   # Subset
#   phase_subset <- ghap.subset(phase, ids = unique(phase$id), variants = this_mkr)
# 
#   #reove low frequency markers
# 
#   # Define haplotype blocks
#   blocks_subset <- ghap.blockgen(phase_subset, windowsize = 3, slide = 2, unit = "marker")
# 
#   # Generate haplotypes and save to file
#   ghap.haplotyping(object = phase_subset,
#                    blocks = blocks_subset,
#                    outfile = paste0("top_var_maf01.", i),
#                    binary = TRUE,
#                    freq = c(0.01, 1))
# 
#   # Load haplotype genotypes
#   haplo <- ghap.loadhaplo(paste0("top_var_maf01.", i))
# 
#   # Convert to PLINK
#   ghap.hap2plink(haplo, outfile = paste0("top_var_maf01.", i))
# }
# 
# 
# assoc_file0maf <- read.table("top_var_maf01.0_met_gwas_sorted_noNA_add", header = TRUE)
# assoc_file1maf <- read.table("top_var_maf01.1_met_gwas_sorted_noNA_add", header = TRUE)
# assoc_file2maf <- read.table("top_var_maf01.2_met_gwas_sorted_noNA_add", header = TRUE)
# assoc_file3maf <- read.table("top_var_maf01.3_met_gwas_sorted_noNA_add", header = TRUE) # no alt
# assoc_file4maf <- read.table("top_var_maf01.4_met_gwas_sorted_noNA_add", header = TRUE) # no alt
# assoc_file5maf <- read.table("top_var_maf01.5_met_gwas_sorted_noNA_add", header = TRUE) # no alt
# 
# 
# 
