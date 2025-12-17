#install.packages("Matrix", type = "source")
#install.packages("susieR")
library(susieR)
library(caret)
library(data.table)
#Load sum stats
# chr5_sumStats <- read_tsv("/home/projects/hearing_loss/clsaARHL_SA/str/chr5/results/chr5_met_association_results.tsv")
# chr5_sumStats <- chr5_sumStats %>% 
#   select(chrom, pos, p_met_phenotype, coeff_met_phenotype,se_met_phenotype)
# 
# #keep STRs within 1mbp of top SNP chr5_73782133
# chr5_sumStats_1mbp <- filter(chr5_sumStats, pos >= 72782133 & pos <= 74782133)
# #remove missing values
# chr5_sumStats_1mbp <- na.omit(chr5_sumStats_1mbp)
# #filtered positions
# chr5_position_to_keep <- as.vector(chr5_sumStats_1mbp$pos)

#########################################
#now get the genotypes (from the dosage data) 

# Step 1: Read raw dosage file
dosage_raw <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/str_finemapping_susieR/chr5_window_1mb_dosage.tsv", header=FALSE)

chr5_position_to_keep <- 72782133:74782133
#keep only values > 0 
dosage_filt <- dosage_raw %>% 
  filter(V2 %in% chr5_position_to_keep)

# Step 2: Function to convert one cell from comma-separated string to numeric sum
convert_cell <- function(cell_str) {
  # Handle missing or empty cells gracefully
  if(is.na(cell_str) || cell_str == "") return(NA_real_)
  # Split string by comma
  nums <- as.numeric(unlist(strsplit(cell_str, ",")))
  # Replace NA with 0 if you want, or leave as NA to propagate missingness
  nums[is.na(nums)] <- 0
  # Sum numeric values to get single dosage per cell
  sum(nums)
}

# Step 3: Extract matrix of dosage strings (all columns except first 3)
dosage_strings <- dosage_filt[, -(1:3)]

# Step 4: Apply conversion over all cells
dosage_numeric <- as.data.frame(lapply(dosage_strings, function(col) sapply(col, convert_cell)))

# Step 5: Assign row names (STR IDs) and column names (sample IDs)
rownames(dosage_numeric) <- dosage_filt[[3]]  # Assuming 3rd column is STR_ID
colnames(dosage_numeric) <- colnames(dosage_filt)[-(1:3)]
#transpose for susie
dosage_numeric_t <- t(dosage_numeric)

# Step 6: (Optional) Convert to matrix for susie
dosage_matrix <- as.matrix(dosage_numeric_t)
n= nrow(dosage_matrix)
dim (dosage_matrix)

# Now dosage_matrix is ready: rows = samples, columns = STRs, values = numeric dosages

#get the phenotype
met_phenotype <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/phenotype/met_merge_filt.txt", header = FALSE)
met_phenotype <- met_phenotype[,1]

#fit susie
fit <- susie(dosage_matrix, met_phenotype, L = 10, verbose = F)
summary(fit)
pip <- as.data.frame(fit$pip)
susie_plot(fit, y = "PIP")


#use summary stats 
#get the LD (correlation) 
#chr5_R <- cor(dosage_matrix)
#high_cor <- findCorrelation(chr5_R, cutoff = 0.9)
 

#try to add SNPs in this region
snp_genotype <- fread("/home/projects/hearing_loss/clsaARHL_SA/str/str_finemapping_susieR/chr5_snps_1mb.raw")
snp_genotype <- snp_genotype[, -(1:6)]
#populate the missing values with the mean 
snp_genotype_imp <- apply(snp_genotype, 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
})
snp_genotype_imp <- as.matrix(snp_genotype_imp)  # ensure it's a matrix
sum(is.na(snp_genotype_imp))
#merge with strs

snp_str_merge <- cbind(scale(snp_genotype_imp), scale(dosage_matrix))
fit2 <- susie(snp_str_merge, met_phenotype, L = 10)
summary(fit2)
pip2 <- as.data.frame(fit2$pip)
susie_plot(fit2, y = "PIP")
write.table(pip2, "/home/projects/hearing_loss/clsaARHL_SA/str/str_finemapping_susieR/pip_str_snp_not_pruned.txt", quote = F)

#get the results
pip2["EnsTR:chr5:73782133", ]
pip2["EnsTR:chr5:73778077", ]
pip2["chr5:73780686",]
pip2["chr5:73776529",] 


#prune variants

# Set correlation threshold
threshold <- 0.90

# Find upper triangle of correlation matrix
cor_mat <- cor(dosage_matrix)
#square to get r2
cor_mat2 <- cor_mat^2
cor_upper <- cor_mat2
cor_upper[lower.tri(cor_upper, diag = TRUE)] <- NA

# Identify highly correlated pairs
high_ld_pairs <- which(abs(cor_upper) > threshold, arr.ind = TRUE)

# Keep only one variant from each highly correlated pair
vars_to_remove <- unique(rownames(cor_upper)[high_ld_pairs[,1]])
# Prune the matrix
dosage_matrix_pruned <- dosage_matrix[, !(colnames(dosage_matrix) %in% vars_to_remove)]


fit3 <- susie(dosage_matrix_pruned, met_phenotype, L = 10, coverage = 0.3)
summary(fit3)
pip3 <- as.data.frame(fit3$pip)
susie_plot(fit3, y = "PIP")

#save results
write.table(pip3, "/home/projects/hearing_loss/clsaARHL_SA/str/str_finemapping_susieR/pip_str.txt", quote = F)


###prune STRs and SNPs
threshold <- 0.90

cor_mat_comb <- cor(snp_str_merge)
#square to get r2
cor_mat_comb2 <- cor_mat_comb^2
cor_upper_comb <- cor_mat_comb2
cor_upper_comb[lower.tri(cor_upper_comb, diag = TRUE)] <- NA
cor_mat_comb3 <- as.data.frame(cor_mat_comb2)

# Identify highly correlated pairs
high_ld_pairs_comb <- which(abs(cor_upper_comb) > threshold, arr.ind = TRUE)

# Keep only one variant from each highly correlated pair
vars_to_remove_comb <- unique(rownames(cor_upper_comb)[high_ld_pairs_comb[,1]])

# Prune the matrix
snp_str_merge_pruned <- snp_str_merge[, !(colnames(snp_str_merge) %in% vars_to_remove_comb)]


fit4 <- susie(snp_str_merge_pruned, met_phenotype, L = 10, coverage = 0.2)
summary(fit4)
pip4 <- as.data.frame(fit4$pip)
susie_plot(fit4, y = "PIP")

write.table(pip4, "/home/projects/hearing_loss/clsaARHL_SA/str/str_finemapping_susieR/pip_str_snp.txt", quote = F)

#get the results
pip3["EnsTR:chr5:73782133", ]
pip4["EnsTR:chr5:73782133", ]

pip3["EnsTR:chr5:73778077", ]
pip4["EnsTR:chr5:73778077", ]


pip4["chr5:73780686",]
pip4["chr5:73776529",] #this SNP was pruned out
pip2["chr5:73776529",] #0.006


#LD between the two top variants
cor_mat_comb3[ "chr5:73780686:C:A_A", "EnsTR:chr5:73782133"] #0.7

cor_mat_comb3["EnsTR:chr5:73778077", "EnsTR:chr5:73754888"]
cor_mat_comb3["EnsTR:chr5:73782133", "EnsTR:chr5:73754888"]


