#load libraries
#library(usethis)
#library(devtools)
library(data.table)
#install.packages("remotes")
#remotes::install_github("juliedwhite/miamiplot")
library(miamiplot)
library(tidyverse)
library(hrbrthemes)
library(broom)
library(gridExtra)
library(ggpubr)
library(dplyr)

metabolic_all <- fread(file = "/home/projects/hearing_loss/clsaARHL_SA/gwas/met_better_gwas_sorted_noNA_add.gz", header = TRUE)
sensory_all <- fread(file = "/home/projects/hearing_loss/clsaARHL_SA/gwas/sen_better_gwas_sorted_noNA_add.gz", header = TRUE)

#add phenotype column
metabolic_all <- mutate(metabolic_all, group = "metabolic")
sensory_all <- mutate(sensory_all, group = "sensory")


#select relevant columns and add group column
metabolic_all.1 <- metabolic_all |> 
  dplyr::select(CHR, BP, SNP, P, group) 

sensory_all.1 <- sensory_all|>
  dplyr::select(CHR, BP, SNP, P, group) 


data_all <- rbind(metabolic_all.1, sensory_all.1)

# # Filter the SNPs below the significance threshold of 1e-3
# data_below_threshold <- data_all %>%
#   filter(P < 5e-3)
# 
# # Filter the SNPs above or equal to the significance threshold of 1e-3
# data_above_threshold <- data_all %>%
#   filter(P >= 5e-3)
# 
# # Set a seed for reproducibility
# set.seed(42)
# 
# # Define the number of SNPs to select from those above the threshold
# num_above_threshold <- 9000  # Adjust this number as needed
# 
# # Randomly select SNPs from those above the threshold
# random_snps_above_threshold <- data_above_threshold %>%
#   sample_n(num_above_threshold, replace = FALSE)
# 
# # Combine the two datasets
# data_all <- bind_rows(data_below_threshold, random_snps_above_threshold)

#get significant snps 
sig_snps_met <- data_all |> 
  filter(group == "metabolic" & P < 1e-5)
sig_snps_met <- sig_snps_met$SNP

sig_snps_sen <- data_all |> 
  filter(group == "sensory" & P < 1e-5)
sig_snps_sen <- sig_snps_sen$SNP


#data_all$rs <- ifelse(data_all$SNP == "chr5:73780686:C:A" & data_all$group == "metabolic", "ARHGEF28",
#                                    ifelse(data_all$SNP == "chr22:50549676:G:A" & data_all$group == "sensory", "KLHDC7B", NA))

#reduce the number of variants for easy visualization
data_all_sub <- filter(data_all, P < 1e-2)

# #rsID <- c("ARHGEF28","KLHDC7B")
ggmiami_gwas2(data = data_all_sub, split_by = "group", split_at = "metabolic", chr = "CHR", pos = "BP", p="P",
              chr_colors = NULL, upper_ylab = "Metabolic",
              lower_ylab = "Sensory", upper_chr_colors = c("grey70", "grey40"),
              lower_chr_colors = c("grey70", "grey40"),suggestive_line = 1e-05,
              suggestive_line_color = "black",genome_line_color = "black",
              upper_highlight_col = "SNP", upper_highlight = sig_snps_met, upper_highlight_color = "steelblue4",
              lower_highlight_col = "SNP", lower_highlight = sig_snps_sen, lower_highlight_color = "darkred")

#####
#add STRs
comb_met_assoctiaTR <- read_tsv ("/home/projects/hearing_loss/clsaARHL_SA/str/results/combined_met_association_results.tsv")
comb_sen_assoctiaTR <- read_tsv ("/home/projects/hearing_loss/clsaARHL_SA/str/results/combined_sen_association_results.tsv")

#add phenotype column
comb_met_assoctiaTR <- mutate(comb_met_assoctiaTR, group = "metabolic")
comb_sen_assoctiaTR <- mutate(comb_sen_assoctiaTR, group = "sensory")


#select relevant columns and add group column
comb_met_assoctiaTR.1 <- comb_met_assoctiaTR |> 
  dplyr::select(chrom, pos, p_met_phenotype, group) |> 
  dplyr::rename(CHR=chrom,
         BP=pos,
         P=p_met_phenotype)

comb_sen_assoctiaTR.1 <- comb_sen_assoctiaTR|>
  dplyr::select(chrom, pos, p_sen_phenotype, group) |> 
  dplyr::rename(CHR=chrom,
                BP=pos,
                P=p_sen_phenotype)

str_data_all <- rbind(comb_met_assoctiaTR.1, comb_sen_assoctiaTR.1)
#remove NAs
str_data_all <- na.omit(str_data_all)
#fix the CHR column
str_data_all$CHR <- as.integer(gsub("^chr", "", str_data_all$CHR))
str_data_all <- as.data.frame(str_data_all)
#create SNP column
str_data_all$SNP <- paste(str_data_all$CHR, str_data_all$BP, sep=":")

#get significant snps 
sig_strs_met <- str_data_all |> 
  filter(group == "metabolic" & P < 1e-5)
sig_strs_met <- sig_strs_met$SNP

sig_strs_sen <- str_data_all |> 
  filter(group == "sensory" & P < 1e-5)
sig_strs_sen <- sig_strs_sen$SNP

#reduce the number of variants for easy visualization
str_data_all_sub <- filter(str_data_all, P < 1e-1)


ggmiami_gwas2(data = combined_data_sub, split_by = "group", split_at = "metabolic", chr = "CHR", pos = "BP", p="P",
              chr_colors = NULL, upper_ylab = "Metabolic",
              lower_ylab = "Sensory", upper_chr_colors = c("grey70", "grey40"),
              lower_chr_colors = c("grey70", "grey40"),suggestive_line = 1e-05, 
              suggestive_line_color = "black",genome_line_color = "black",
              upper_highlight_col = "SNP", upper_highlight = sig_met, upper_highlight_color = c("orange", "steelblue3"),
              lower_highlight_col = "SNP", lower_highlight = sig_sen, lower_highlight_color = c("orange", "darkred"))

#combine the datasets
combined_data_sub <- rbind(data_all_sub, str_data_all_sub)

#combine the significant strs and snos for each phenotype
sig_met <- c(sig_strs_met, sig_snps_met)
sig_sen <- c(sig_strs_sen, sig_snps_sen)

upper_highlight_snps <- c(sig_snps_met, sig_strs_met)
upper_highlight_colors <- c(rep("steelblue3", length(sig_snps_met)),
                            rep("orange", length(sig_strs_met)))
lower_highlight_snps <- c(sig_snps_sen, sig_strs_sen)
lower_highlight_colors <- c(rep("darkred", length(sig_snps_sen)),
                            rep("orange", length(sig_strs_sen)))

ggmiami_gwas2(data = combined_data_sub, split_by = "group", split_at = "metabolic", chr = "CHR", pos = "BP", p = "P",
              chr_colors = NULL,
              upper_ylab = "Metabolic",
              lower_ylab = "Sensory",
              upper_chr_colors = c("grey70", "grey40"),
              lower_chr_colors = c("grey70", "grey40"),
              suggestive_line = 1e-5,
              suggestive_line_color = "black",
              genome_line_color = "black",
              upper_highlight_col = "SNP",
              upper_highlight = upper_highlight_snps,
              upper_highlight_color = upper_highlight_colors,
              lower_highlight_col = "SNP",
              lower_highlight = lower_highlight_snps,
              lower_highlight_color = lower_highlight_colors)

#############################################
## ggmiami function
ggmiami_gwas2 <-function (data, split_by, split_at, chr = "chr", pos = "pos", 
                          p = "p", chr_colors = c("black", "grey"), upper_chr_colors = NULL, 
                          lower_chr_colors = NULL, upper_ylab = "-log10(p)", lower_ylab = "-log10(p)", 
                          genome_line = 5e-08, genome_line_color = "red", suggestive_line = 1e-05, 
                          suggestive_line_color = "blue", hits_label_col = NULL, hits_label = NULL, 
                          top_n_hits = 5, upper_labels_df = NULL, lower_labels_df = NULL, 
                          upper_highlight = NULL, upper_highlight_col = NULL, upper_highlight_color = "green", 
                          lower_highlight = NULL, lower_highlight_col = NULL, lower_highlight_color = "green") 
{
  plot_data <- prep_miami_data(data = data, split_by = split_by, 
                               split_at = split_at, chr = chr, pos = pos, p = p)
  # Create a vector for x-axis labels with chromosome 23 labeled as "X"
  custom_labels <- plot_data$axis$chr
  custom_labels[custom_labels == 23] <- "X"
  
  if (all(!is.null(chr_colors), any(!is.null(upper_chr_colors), 
                                    !is.null(lower_chr_colors)))) {
    stop("You have specified both chr_colors and upper_chr_colors and/or\n         lower_chr_colors. This package does not know how to use both\n         information simultaneously. Please only use one method for coloring:\n         either chr_colors, for making upper and lower plot have the same\n         colors, or upper_chr_colors + lower_chr_colors for specifying different\n         colors for upper and lower plot.")
  }
  if (upper_ylab == "-log10(p)") {
    upper_ylab <- expression("-log"[10] * "(p)")
  }
  else {
    upper_ylab <- bquote(atop(.(upper_ylab), "-log"[10] * 
                                "(p)"))
  }
  if (lower_ylab == "-log10(p)") {
    lower_ylab <- expression("-log"[10] * "(p)")
  }
  else {
    lower_ylab <- bquote(atop(.(lower_ylab), "-log"[10] * 
                                "(p)"))
  }
  upper_plot <- ggplot2::ggplot() + ggplot2::geom_point(data = plot_data$upper, 
                                                        aes(x = .data$rel_pos, y = .data$logged_p, color = as.factor(.data$chr)), 
                                                        size = 3) + ggplot2::scale_x_continuous(labels = custom_labels, 
                                                                                                breaks = plot_data$axis$chr_center, 
                                                                                                expand = ggplot2::expansion(mult = 0.01), 
                                                                                                guide = ggplot2::guide_axis(check.overlap = TRUE)) + 
    ggplot2::scale_y_continuous(limits = c(0, 12), 
                                expand = ggplot2::expansion(mult = c(0.02, 0))) + 
    ggplot2::labs(x = "", y = upper_ylab) + ggplot2::theme_classic() + 
    ggplot2::theme(legend.position = "none", axis.title.y = ggplot2::element_text(size = 20), 
                   axis.title.x = ggplot2::element_blank(), axis.text = element_text(size = 15),
                   plot.margin = ggplot2::margin(t = 10, b = 0, l = 10, r = 10))
  lower_plot <- ggplot2::ggplot() + ggplot2::geom_point(data = plot_data$lower, 
                                                        aes(x = .data$rel_pos, y = .data$logged_p, color = as.factor(.data$chr)), 
                                                        size = 3) + ggplot2::scale_x_continuous(breaks = plot_data$axis$chr_center,
                                                                                                position = "top", expand = ggplot2::expansion(mult = 0.01)) + 
    ggplot2::scale_y_reverse(limits = c(12, 0), 
                             expand = ggplot2::expansion(mult = c(0, 0.02))) + 
    ggplot2::labs(x = "", y = lower_ylab) + ggplot2::theme_classic() + 
    ggplot2::theme(legend.position = "none", axis.title.y = ggplot2::element_text(size = 20), 
                   axis.text.x = ggplot2::element_blank(), axis.text = element_text(size = 15),
                   axis.title.x = ggplot2::element_blank(), plot.margin = ggplot2::margin(t = 0, 
                                                                                          b = 10, l = 10, r = 10))
  if (all(!is.null(chr_colors), is.null(upper_chr_colors), 
          is.null(lower_chr_colors))) {
    if (length(chr_colors) == 2) {
      chr_colors <- rep(chr_colors, length.out = nrow(plot_data$axis))
    }
    else if (length(chr_colors) == nrow(plot_data$axis)) {
      chr_colors <- chr_colors
    }
    else {
      stop("The number of colors specified in {chr_colors} does not match the\n         number of chromosomes to be displayed.")
    }
    upper_plot <- upper_plot + ggplot2::scale_color_manual(values = chr_colors)
    lower_plot <- lower_plot + ggplot2::scale_color_manual(values = chr_colors)
  }
  else if (all(is.null(chr_colors), !is.null(upper_chr_colors), 
               !is.null(lower_chr_colors))) {
    if (length(upper_chr_colors) == 2) {
      upper_chr_colors <- rep(upper_chr_colors, length.out = nrow(plot_data$axis))
    }
    else if (length(upper_chr_colors) == nrow(plot_data$axis)) {
      upper_chr_colors <- upper_chr_colors
    }
    else {
      stop("The number of colors specified in {upper_chr_colors} does not match\n           the number of chromosomes to be displayed.")
    }
    if (length(lower_chr_colors) == 2) {
      lower_chr_colors <- rep(lower_chr_colors, length.out = nrow(plot_data$axis))
    }
    else if (length(lower_chr_colors) == nrow(plot_data$axis)) {
      lower_chr_colors <- lower_chr_colors
    }
    else {
      stop("The number of colors specified in {lower_chr_colors} does not match\n           the number of chromosomes to be displayed.")
    }
    upper_plot <- upper_plot + ggplot2::scale_color_manual(values = upper_chr_colors)
    lower_plot <- lower_plot + ggplot2::scale_color_manual(values = lower_chr_colors)
  }
  else if (all(is.null(chr_colors), any(is.null(upper_chr_colors), 
                                        is.null(lower_chr_colors)))) {
    stop("It looks like you've specified one of upper or lower chr colors\n         without specifying the other. This package needs both colors, unless\n         you want the upper and lower plot to have the same colors, which is\n         done using {chr_colors}.")
  }
  if (!is.null(suggestive_line)) {
    upper_plot <- upper_plot + ggplot2::geom_hline(yintercept = -log10(suggestive_line), 
                                                   color = suggestive_line_color, linetype = "dotted", 
                                                   linewidth = 0.5)
    lower_plot <- lower_plot + ggplot2::geom_hline(yintercept = -log10(suggestive_line), 
                                                   color = suggestive_line_color, linetype = "dotted", 
                                                   linewidth = 0.5)
  }
  if (!is.null(genome_line)) {
    upper_plot <- upper_plot + ggplot2::geom_hline(yintercept = -log10(genome_line), 
                                                   color = genome_line_color, linetype = "dashed", linewidth = 0.7)
    lower_plot <- lower_plot + ggplot2::geom_hline(yintercept = -log10(genome_line), 
                                                   color = genome_line_color, linetype = "dashed", linewidth = 0.7)
  }
  if (all(!is.null(hits_label_col), any(!is.null(upper_labels_df), 
                                        !is.null(lower_labels_df)))) {
    stop("You have specified both hits_label_col and a *_labels_df. This\n         package does not know how to use both information simultaneously.\n         Please only use one method for labelling: either hits_label_col (with\n         or without hits_label), or *_labels_df.")
  }
  if (all(!is.null(hits_label_col), is.null(upper_labels_df), 
          is.null(lower_labels_df))) {
    upper_labels_df <- make_miami_labels(data = plot_data$upper, 
                                         hits_label_col = hits_label_col, hits_label = hits_label, 
                                         top_n_hits = top_n_hits)
    lower_labels_df <- make_miami_labels(data = plot_data$lower, 
                                         hits_label_col = hits_label_col, hits_label = hits_label, 
                                         top_n_hits = top_n_hits)
    upper_plot <- upper_plot + ggrepel::geom_label_repel(data = upper_labels_df, 
                                                         aes(x = .data$rel_pos, y = .data$logged_p, label = .data$label), 
                                                         size = 5, segment.size = 0.2, point.padding = 0.3, 
                                                         ylim = c(plot_data$maxp/2, NA), min.segment.length = 0, 
                                                         force = 2, box.padding = 0.5, fontface = "italic")
    lower_plot <- lower_plot + ggrepel::geom_label_repel(data = lower_labels_df, 
                                                         aes(x = .data$rel_pos, y = .data$logged_p, label = .data$label), 
                                                         size = 5, segment.size = 0.2, point.padding = 0.3, 
                                                         ylim = c(NA, -(plot_data$maxp/2)), min.segment.length = 0, 
                                                         force = 2, box.padding = 0.5, fontface = "italic")
  }
  if (all(is.null(hits_label_col), !is.null(upper_labels_df))) {
    checkmate::assertNames(colnames(upper_labels_df), identical.to = c("rel_pos", 
                                                                       "logged_p", "label"))
    upper_plot <- upper_plot + ggrepel::geom_label_repel(data = upper_labels_df, 
                                                         aes(x = .data$rel_pos, y = .data$logged_p, label = .data$label), 
                                                         size = 2, segment.size = 0.2, point.padding = 0.3, 
                                                         ylim = c(plot_data$maxp/2, NA), min.segment.length = 0, 
                                                         force = 2, box.padding = 0.5, fontface = "italic",
                                                         max.overlaps = 50)
  }
  if (all(is.null(hits_label_col), !is.null(lower_labels_df))) {
    checkmate::assertNames(colnames(lower_labels_df), identical.to = c("rel_pos", 
                                                                       "logged_p", "label"))
    lower_plot <- lower_plot + ggrepel::geom_label_repel(data = lower_labels_df, 
                                                         aes(x = .data$rel_pos, y = .data$logged_p, label = .data$label), 
                                                         size = 2, segment.size = 0.2, point.padding = 0.3, 
                                                         ylim = c(NA, -(plot_data$maxp/2)), min.segment.length = 0, 
                                                         force = 2, box.padding = 0.5, fontface = "italic",
                                                         max.overlaps = 50)
  }
  if (all(!is.null(upper_highlight), !is.null(upper_highlight_col))) {
    upper_highlight_df <- highlight_miami(data = plot_data$upper, 
                                          highlight = upper_highlight, highlight_col = upper_highlight_col, 
                                          highlight_color = upper_highlight_color)
    upper_plot <- upper_plot + ggplot2::geom_point(data = upper_highlight_df, 
                                                   aes(x = .data$rel_pos, y = .data$logged_p), color = upper_highlight_df$color, 
                                                   size = 3)
  }
  if (all(!is.null(lower_highlight), !is.null(lower_highlight_col))) {
    lower_highlight_df <- highlight_miami(data = plot_data$lower, 
                                          highlight = lower_highlight, highlight_col = lower_highlight_col, 
                                          highlight_color = lower_highlight_color)
    lower_plot <- lower_plot + ggplot2::geom_point(data = lower_highlight_df, 
                                                   aes(x = .data$rel_pos, y = .data$logged_p), color = lower_highlight_df$color, 
                                                   size = 3)
  }
  gridExtra::grid.arrange(upper_plot, lower_plot, nrow = 2)
}


#####
