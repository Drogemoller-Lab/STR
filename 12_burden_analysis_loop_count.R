#!/usr/bin/env Rscript
.libPaths("/home/ahmeds26@med.umanitoba.ca/R/x86_64-pc-linux-gnu-library/4.4")

#!/usr/bin/env Rscript

# Complete STR expansion counting analysis
# Part 1: Process each chromosome and save individual burden files
# Part 2: Read and combine all chromosome files

library(vcfR)
library(data.table)
library(stringr)
library(tidyr)

# Function to get allele length given genotype index
get_allele_length <- function(idx, ref_len, alt_lens, motif_len) {
  if (is.na(idx)) return(NA_real_)
  if (idx == 0) return(ref_len / motif_len)
  if (idx > length(alt_lens) || is.na(alt_lens[idx])) return(NA_real_)
  return(alt_lens[idx] / motif_len)
}

# Define analysis parameters
length_thresholds <- c(1, 5, 10, 20)  # repeat units longer than reference
frequency_cutoffs <- c(1, 5, 10, 100, Inf)  # expansion frequency cutoffs (Inf = no cutoff)
frequency_names <- c("once", "le5", "le10", "le100", "no_cutoff")

# UPDATED FUNCTION: Count number of STR expansions (not sum lengths)
calculate_str_expansion_count <- function(longer_allele_matrix, ref_lengths, period, length_thresh, freq_cutoff) {
  if (is.data.frame(longer_allele_matrix)) {
    longer_allele_matrix <- as.matrix(longer_allele_matrix)
  }
  
  # Filter valid variants
  valid_variants <- !is.na(period) & period > 0 & !is.na(ref_lengths) & ref_lengths > 0
  if (sum(valid_variants) == 0) {
    return(rep(0, ncol(longer_allele_matrix)))
  }
  
  longer_allele_subset <- longer_allele_matrix[valid_variants, , drop = FALSE]
  period_subset <- period[valid_variants]
  ref_lengths_subset <- ref_lengths[valid_variants]
  ref_lengths_ru <- ref_lengths_subset / period_subset
  
  # Calculate expansion sizes
  expansion_matrix <- sweep(longer_allele_subset, 1, ref_lengths_ru, FUN = "-")
  
  # Apply length threshold filter (≥threshold RU longer than reference)
  expanded_sites <- expansion_matrix >= length_thresh
  expanded_sites[is.na(expanded_sites)] <- FALSE
  
  # Apply frequency filter
  if (!is.infinite(freq_cutoff)) {
    site_frequencies <- rowSums(expanded_sites, na.rm = TRUE)
    sites_to_keep <- site_frequencies <= freq_cutoff
    expanded_sites[!sites_to_keep, ] <- FALSE
  }
  
  # COUNT number of expansion sites, not sum lengths
  expansion_count_per_sample <- colSums(expanded_sites, na.rm = TRUE)
  
  return(expansion_count_per_sample)
}

# Define chromosomes and output directory
chromosomes <- 1:22
out_dir <- "/home/projects/hearing_loss/clsaARHL_SA/str/burden_analysis/"

cat("=== PART 1: PROCESSING INDIVIDUAL CHROMOSOMES ===\n")
cat("Counting STR expansion sites (not summing repeat unit lengths)\n\n")

# Loop through each chromosome and save individual files
for (chr in chromosomes) {
  cat("Processing chromosome", chr, "...\n")
  
  # Build file paths
  vcf_file <- paste0("/home/projects/hearing_loss/clsaARHL_SA/str/chr", chr,
                     "/input/annotated_filtered_imputed_STR_SNPs_chr", chr, ".vcf.gz")
  burden_file <- paste0(out_dir, "chr", chr, "_expansion_count.tsv")
  
  # Skip if burden file already exists
  if (file.exists(burden_file)) {
    cat("Count file already exists:", burden_file, ". Skipping...\n")
    next
  }
  
  # Check if VCF file exists
  if (!file.exists(vcf_file)) {
    cat("Warning: VCF file does not exist:", vcf_file, "\n")
    next
  }
  
  tryCatch({
    # Read and process VCF file
    cat("  Reading VCF file...\n")
    vcf <- read.vcfR(vcf_file)
    
    ref <- getREF(vcf)
    alt <- getALT(vcf)
    period <- as.numeric(extract.info(vcf, "PERIOD"))
    
    gt <- as.data.frame(extract.gt(vcf))
    gt_matrix <- as.matrix(gt)
    
    cat("  Dataset:", nrow(gt_matrix), "variants x", ncol(gt_matrix), "samples\n")
    
    # Fix sample IDs if they're in "1", "2" format
    sample_ids <- colnames(gt)
    if (all(grepl("^[0-9]+$", sample_ids))) {
      colnames(gt_matrix) <- paste0(sample_ids, "_", sample_ids)
      colnames(gt) <- paste0(sample_ids, "_", sample_ids)
      cat("  Fixed sample IDs to 1_1 format\n")
    }
    
    # Parse genotypes
    cat("  Parsing genotypes...\n")
    gt1_matrix <- matrix(NA_integer_, nrow = nrow(gt_matrix), ncol = ncol(gt_matrix))
    gt2_matrix <- matrix(NA_integer_, nrow = nrow(gt_matrix), ncol = ncol(gt_matrix))
    
    for (j in seq_len(ncol(gt_matrix))) {
      if (j %% 1000 == 0) cat("    Processing sample", j, "of", ncol(gt_matrix), "\n")
      
      splits <- strsplit(gt_matrix[, j], "\\|", fixed = FALSE)
      splits_matrix <- do.call(rbind, splits)
      
      gt1_matrix[, j] <- as.integer(splits_matrix[, 1])
      gt2_matrix[, j] <- as.integer(splits_matrix[, 2])
    }
    
    colnames(gt1_matrix) <- colnames(gt)
    colnames(gt2_matrix) <- colnames(gt)
    
    # Calculate allele lengths
    cat("  Calculating allele lengths...\n")
    ref_lengths <- nchar(as.character(ref))
    alt_list <- strsplit(as.character(alt), ",")
    max_alts <- max(lengths(alt_list))
    
    alt_lengths <- matrix(NA_real_, nrow = length(alt_list), ncol = max_alts)
    for (i in seq_along(alt_list)) {
      if (length(alt_list[[i]]) > 0) {
        alt_lengths[i, seq_along(alt_list[[i]])] <- nchar(alt_list[[i]])
      }
    }
    
    allele1_len <- matrix(NA_real_, nrow = nrow(gt1_matrix), ncol = ncol(gt1_matrix))
    allele2_len <- matrix(NA_real_, nrow = nrow(gt2_matrix), ncol = ncol(gt2_matrix))
    
    for (i in seq_len(nrow(gt1_matrix))) {
      if (i %% 1000 == 0) cat("    Processing variant", i, "of", nrow(gt1_matrix), "\n")
      
      motif_len <- period[i]
      if (is.na(motif_len) || motif_len == 0) next
      
      ref_len <- ref_lengths[i]
      alt_lens <- alt_lengths[i, ]
      
      allele1_len[i, ] <- vapply(gt1_matrix[i, ], get_allele_length,
                                 numeric(1), ref_len, alt_lens, motif_len)
      allele2_len[i, ] <- vapply(gt2_matrix[i, ], get_allele_length,
                                 numeric(1), ref_len, alt_lens, motif_len)
    }
    
    # Take the longer allele
    cat("  Computing longer alleles...\n")
    longer_allele <- pmax(allele1_len, allele2_len, na.rm = TRUE)
    colnames(longer_allele) <- colnames(gt)
    
    # Calculate expansion counts for all combinations
    cat("  Counting STR expansions for all threshold/frequency combinations...\n")
    chr_results <- data.frame(sample_id = colnames(gt), stringsAsFactors = FALSE)
    
    for (i in seq_along(length_thresholds)) {
      len_thresh <- length_thresholds[i]
      cat("    Length threshold >=", len_thresh, "RU\n")
      
      for (j in seq_along(frequency_cutoffs)) {
        freq_cutoff <- frequency_cutoffs[j]
        freq_name <- frequency_names[j]
        
        # Use counting function instead of length summing
        expansion_count <- calculate_str_expansion_count(longer_allele, ref_lengths, period, len_thresh, freq_cutoff)
        
        col_name <- paste0("chr", chr, "_ge", len_thresh, "_", freq_name)
        chr_results[col_name] <- expansion_count  # No rounding needed for counts
        
        cat("      ", freq_name, ": mean =", round(mean(expansion_count), 1), "expansions\n")
      }
    }
    
    # Set row names and save
    rownames(chr_results) <- chr_results$sample_id
    chr_results$sample_id <- NULL
    
    write.table(chr_results, burden_file, sep="\t", quote=FALSE, row.names=TRUE)
    
    cat("Chromosome", chr, "expansion count file saved:", burden_file, "\n")
    cat("Dimensions:", nrow(chr_results), "samples x", ncol(chr_results), "expansion count measures\n\n")
    
    # Clean up memory
    rm(vcf, gt, gt_matrix, gt1_matrix, gt2_matrix, longer_allele, chr_results)
    gc()
    
  }, error = function(e) {
    cat("ERROR processing chromosome", chr, ":", e$message, "\n")
    cat("Continuing with next chromosome...\n\n")
  })
}

cat("=== PART 2: COMBINING ALL CHROMOSOME COUNT FILES ===\n")

# Read and combine all count files
count_files <- paste0(out_dir, "chr", chromosomes, "_expansion_count.tsv")
existing_files <- count_files[file.exists(count_files)]

if (length(existing_files) == 0) {
  stop("No expansion count files found to combine!")
}

cat("Found", length(existing_files), "expansion count files to combine:\n")
for (f in existing_files) {
  cat("  ", basename(f), "\n")
}

# Read first file to get sample structure
cat("\nReading first file to initialize combined results...\n")
first_file <- existing_files[1]
combined_results <- read.table(first_file, sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)

cat("Sample structure from first file:", nrow(combined_results), "samples x", ncol(combined_results), "columns\n")

# Read and merge remaining files
for (i in 2:length(existing_files)) {
  file_path <- existing_files[i]
  cat("Reading and merging", basename(file_path), "...\n")

  chr_data <- read.table(file_path, sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)

  # Check sample consistency
  if (!identical(rownames(combined_results), rownames(chr_data))) {
    cat("  Warning: Sample order differs. Reordering to match...\n")
    chr_data <- chr_data[rownames(combined_results), , drop=FALSE]
  }

  # Add columns from this chromosome
  combined_results <- cbind(combined_results, chr_data)
}

# Calculate genome-wide totals for each threshold/frequency combination
cat("\nCalculating genome-wide expansion count totals...\n")

for (i in seq_along(length_thresholds)) {
  len_thresh <- length_thresholds[i]

  for (j in seq_along(frequency_cutoffs)) {
    freq_name <- frequency_names[j]

    # Find all chromosome columns for this combination
    pattern <- paste0("chr[0-9]+_ge", len_thresh, "_", freq_name, "$")
    chr_cols <- grep(pattern, names(combined_results), value=TRUE)

    if (length(chr_cols) > 0) {
      total_col_name <- paste0("total_ge", len_thresh, "_", freq_name)
      combined_results[total_col_name] <- rowSums(combined_results[, chr_cols, drop=FALSE], na.rm=TRUE)
      cat("  ", total_col_name, ": summed", length(chr_cols), "chromosomes\n")
    }
  }
}

# Reorder columns: totals first, then chromosomes
total_cols <- grep("^total_", names(combined_results), value=TRUE)
chr_cols <- grep("^chr[0-9]+_", names(combined_results), value=TRUE)

# Sort columns
total_cols <- sort(total_cols)
chr_cols <- sort(chr_cols)

final_cols <- c(total_cols, chr_cols)
combined_results <- combined_results[, final_cols]

# Save final combined results
final_output <- paste0(out_dir, "comprehensive_str_expansion_counts.tsv")
write.table(combined_results, final_output, sep="\t", quote=FALSE, row.names=TRUE)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Final combined results saved to:", final_output, "\n")
cat("Total samples:", nrow(combined_results), "\n")
cat("Total columns:", ncol(combined_results), "\n")
cat("Chromosomes processed:", length(existing_files), "\n")

# Print summary statistics for genome-wide totals
cat("\nGenome-wide expansion count summary:\n")
for (col in total_cols[1:min(10, length(total_cols))]) {
  values <- combined_results[, col]
  cat(col, ": mean =", round(mean(values, na.rm=TRUE), 1),
      ", median =", round(median(values, na.rm=TRUE), 1), "expansions\n")
}

cat("\nExample results (first 3 samples, first 5 columns):\n")
print(combined_results[1:3, 1:min(5, ncol(combined_results))])

cat("\nProcessing completed successfully!\n")
cat("Results represent COUNTS of STR expansion sites, not repeat unit lengths.\n")


##############
combined_results <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/burden_analysis/comprehensive_str_expansion_counts.tsv")

# Each column follows this pattern: chr[X]_ge[Y]_[Z]
# chr[X] = chromosome number (e.g., chr22)
# ge[Y] = "greater than or equal to Y repeat units" longer than reference
# [Z] = frequency cutoff category

#select the total
comp_burd <- combined_results %>%
  select(starts_with("total_")) %>%
  select(-contains("once")) 

comp_burd$FID <- row.names(comp_burd)

#to generate the bar plot
met_phen <- read.table("/home/projects/hearing_loss/clsaARHL_SA/phenotypes/metabolic.txt", header = TRUE)

#select relevant columns
met_phen <- met_phen %>%
  select(FID, better.rank)

#generate the ID column in this format (1_1)
met_phen$FID <- paste0(met_phen$FID, "_", met_phen$FID)

#merge with sum of allele length
comp_burd_met <- merge(comp_burd, met_phen, by = "FID")

burden_cols <- grep("^total_", names(comp_burd_met), value = TRUE)

assoc_results_met <- data.frame(
  burden = burden_cols,
  beta = NA, se = NA, p = NA
)

for (i in seq_along(burden_cols)) {
  f <- as.formula(paste("better.rank ~", burden_cols[i]))
  fit <- lm(f, data = comp_burd_met)
  s <- summary(fit)$coefficients[2, ]
  assoc_results_met[i, 2:4] <- c(s["Estimate"], s["Std. Error"], s["Pr(>|t|)"])
}

assoc_results_met

# FDR correction
assoc_results_met$FDR <- p.adjust(assoc_results_met$p, method = "fdr")
#scale but not centre
assoc_results_met$beta_scale <- scale(assoc_results_met$beta, center = FALSE)


#then make the bar plot
assoc_results_met_to_plot <- assoc_results_met %>%
  mutate(
    units = sub(".*_ge([0-9]+)_.*", "\\1", burden),
    frequency = sub(".*_ge[0-9]+_(.*)", "\\1", burden)
  )

#ensure that the frequency and units are factors with levels

assoc_results_met_to_plot$frequency <- factor(assoc_results_met_to_plot$frequency, levels = c("once", "le5", "le10", "le100", "no_cutoff"))
assoc_results_met_to_plot$units <- factor(assoc_results_met_to_plot$units, levels = c("1", "5", "10", "20"))
assoc_results_met_to_plot$unit_label <- factor(
  assoc_results_met_to_plot$units,
  levels = c(1, 5, 10, 20),
  labels = c("≥1", "≥5", "≥10", "≥20")
)

met_plot <- assoc_results_met_to_plot %>% 
ggplot(aes(x = unit_label, y = beta_scale, fill = frequency)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"), width = 0.8, color = "grey40") +
  scale_y_continuous(limits = c(0, 3)) +
  scale_fill_manual(
    name = "Frequency",
    values = c(
      "le5" = "#CFEAF2",
      "le10" = "#1E90FF",
      "le100" = "steelblue4",
      "no_cutoff" = "#1F3752"
    ),
    labels = c("≤5", "≤10", "≤100", "No limit")
  ) +
  labs(
    x = " ",
    y = "Burden of STRs on \nmetabolic hearing loss (beta from regression)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

####
#try the sensory

sen_phen <- read.table("/home/projects/hearing_loss/clsaARHL_SA/phenotypes/sensory.txt", header = TRUE)

#select relevant columns
sen_phen <- sen_phen %>%
  select(FID, better.rank)

#generate the ID column in this format (1_1)
sen_phen$FID <- paste0(sen_phen$FID, "_", sen_phen$FID)

#merge with sum of allele length
comp_burd_sen <- merge(comp_burd, sen_phen, by = "FID")

burden_cols <- grep("^total_", names(comp_burd_sen), value = TRUE)

assoc_results_sen <- data.frame(
  burden = burden_cols,
  beta = NA, se = NA, p = NA
)

for (i in seq_along(burden_cols)) {
  f <- as.formula(paste("better.rank ~", burden_cols[i]))
  fit <- lm(f, data = comp_burd_sen)
  s <- summary(fit)$coefficients[2, ]
  assoc_results_sen[i, 2:4] <- c(s["Estimate"], s["Std. Error"], s["Pr(>|t|)"])
}

assoc_results_sen

# FDR correction
assoc_results_sen$FDR <- p.adjust(assoc_results_sen$p, method = "fdr")
#scale but not centre
assoc_results_sen$beta_scale <- scale(assoc_results_sen$beta, center = FALSE)

#then make the bar plot
assoc_results_sen_to_plot <- assoc_results_sen %>%
  mutate(
    units = sub(".*_ge([0-9]+)_.*", "\\1", burden),
    frequency = sub(".*_ge[0-9]+_(.*)", "\\1", burden)
  )

#ensure that the frequency and units are factors with levels

assoc_results_sen_to_plot$frequency <- factor(assoc_results_sen_to_plot$frequency, levels = c("once", "le5", "le10", "le100", "no_cutoff"))
assoc_results_sen_to_plot$units <- factor(assoc_results_sen_to_plot$units, levels = c("1", "5", "10", "20"))
assoc_results_sen_to_plot$unit_label <- factor(
  assoc_results_sen_to_plot$units,
  levels = c(1, 5, 10, 20),
  labels = c("≥1", "≥5", "≥10", "≥20")
)

sen_plot <- ggplot(assoc_results_sen_to_plot, aes(x = unit_label, y = beta_scale, fill = frequency)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"), width = 0.8, color = "grey40") +
  scale_fill_manual(
    name = "Frequency",
    values = c(
      "le5" = "#E3C6C6",
      "le10" = "#C17B7B",
      "le100" = "darkred",
      "no_cutoff" = "#500000"
    ),
    labels = c("≤5", "≤10", "≤100", "No limit")
  ) +
  scale_x_discrete(position = "top") +  # Move x-axis to top
  scale_y_continuous(limits = c(NA, 0)) +  # Ensure y-axis goes up to 0
  labs(
    x = " ",
    y = "Burden of STRs on \nsensory hearing loss (beta from regression)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    axis.line.x.bottom = element_blank(),  # Remove the bottom x-axis line
    axis.line.x.top = element_line(color = "black")
  )


#merge results tables
assoc_results_merge <- merge(assoc_results_met, assoc_results_sen, by = "burden")

#combine the two plots
library(ggpubr)

combined_burden_plots <- ggarrange(met_plot, sen_plot, 
                                   ncol = 1, 
                                   common.legend = FALSE, 
                                   legend = "right")

# add a shared x label
combined_burden_plots <- annotate_figure(
  combined_burden_plots,
  bottom = text_grob("Number of repeat units longer than reference", 
                     size = 14)
)

tiff("/home/projects/hearing_loss/clsaARHL_SA/plots/combined_burden_plots.tiff",
     width = 8, height = 12, units = "in", res = 300)
print(combined_burden_plots)  # must use print() to render
dev.off()


#