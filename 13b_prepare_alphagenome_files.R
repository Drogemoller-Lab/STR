#make alphagenome compatible .vcf

#read bi-allelic .bim file
#V1: chromosome
#V2: ID
#v3: 0
#V4: POS
#V5: ALT
#V6: REF

biallelic_bim <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/heritability/snp_str_comb_heritability_new/input/chrAll_str_biallelic_filtered.bim")

#now get the list of suggestive variants
comb_met_assoctiaTR <- read_tsv ("/home/projects/hearing_loss/clsaARHL_SA/str/results/combined_met_association_results.tsv")
comb_sen_assoctiaTR <- read_tsv ("/home/projects/hearing_loss/clsaARHL_SA/str/results/combined_sen_association_results.tsv")

#filter suggestive significant variants and the fine mapped variant
met_str_sig_sugest <- comb_met_assoctiaTR %>% 
  filter (p_met_phenotype < 1e-5 | pos == 73778077) %>% 
  select(chrom, pos)
#edit the chromosome column
biallelic_bim$V1 <- paste0("chr", biallelic_bim$V1)

#now extract the variants that match this file
biallelic_bim_filt <- merge(biallelic_bim, met_str_sig_sugest, by.x = c("V1", "V4"), by.y = c("chrom", "pos"))

#generate variant id column chr_pos_uniqnumber
# Create basic variant id
biallelic_bim_filt <- biallelic_bim_filt %>%
  mutate(variant_id_base = paste0(V1, "_", V4)) %>%
  group_by(variant_id_base) %>%
  mutate(variant_id = if (n() == 1) {
    variant_id_base
  } else {
    paste0(variant_id_base, "_", row_number())
  }) %>%
  ungroup() %>%
  select(-variant_id_base)

# select and rename
biallelic_bim_clean <- biallelic_bim_filt %>% 
  select(variant_id, V1, V4, V6, V5) %>% 
  rename(CHROM = V1,
         POS = V4,
         REF = V6,
         ALT = V5)

#save the file as .vcf
write.table(biallelic_bim_clean, "/home/projects/hearing_loss/clsaARHL_SA/str/AlphaGenome/met_suggestive.vcf", quote = F, row.names = F,  sep = "\t")
########


######
#for sensory phenotype
#filter suggestive significant variants and the fine mapped variant
sen_str_sig_sugest <- comb_sen_assoctiaTR %>% 
  filter (p_sen_phenotype < 1e-5 | pos == 73778077) %>% 
  select(chrom, pos)

#now extract the variants that match this file
sen_biallelic_bim_filt <- merge(biallelic_bim, sen_str_sig_sugest, by.x = c("V1", "V4"), by.y = c("chrom", "pos"))

#generate variant id column chr_pos_uniqnumber
# Create basic variant id
sen_biallelic_bim_filt <- sen_biallelic_bim_filt %>%
  mutate(variant_id_base = paste0(V1, "_", V4)) %>%
  group_by(variant_id_base) %>%
  mutate(variant_id = if (n() == 1) {
    variant_id_base
  } else {
    paste0(variant_id_base, "_", row_number())
  }) %>%
  ungroup() %>%
  select(-variant_id_base)

# select and rename
sen_biallelic_bim_clean <- sen_biallelic_bim_filt %>% 
  select(variant_id, V1, V4, V6, V5) %>% 
  rename(CHROM = V1,
         POS = V4,
         REF = V6,
         ALT = V5)

#save the file as .vcf
write.table(sen_biallelic_bim_clean, "/home/projects/hearing_loss/clsaARHL_SA/str/AlphaGenome/sen_suggestive.vcf", quote = F, row.names = F,  sep = "\t")
########
#read results
met_scores_batch <- read.csv("/home/projects/hearing_loss/clsaARHL_SA/str/AlphaGenome/met_STRs_batch_scores_100kb.csv")
sen_scores_batch <- read.csv("/home/projects/hearing_loss/clsaARHL_SA/str/AlphaGenome/sen_STRs_batch_scores_100kb.csv")

#get the annotation for the STRs of interest in the brain UBERON:0000955
met_scores_batch_brain <- met_scores_batch %>% 
  filter(ontology_curie == "UBERON:0000955")


met_scores_batch_brain_chr5_73782133 <- met_scores_batch_brain[grepl("chr5:73782133", met_scores_batch_brain[[1]]), ]
#try top 1%
met_scores_batch_brain_chr5_73782133_sig <- met_scores_batch_brain_chr5_73782133 %>% 
  filter(abs (quantile_score) > 0.99 & gene_name == "ARHGEF28")

#####################
#STR with xQTL annotation chr5:73754888

met_scores_batch_brain_chr5_73754888 <- met_scores_batch_brain[grepl("chr5:73754888", met_scores_batch_brain[[1]]), ]
#try top 1%
met_scores_batch_brain_chr5_73754888_sig <- met_scores_batch_brain_chr5_73754888 %>% 
  filter(abs (quantile_score) > 0.99 & gene_name == "ARHGEF28")

#STR with highest PIP chr5:73778077
met_scores_batch_brain_chr5_73778077 <- met_scores_batch_brain[grepl("chr5:73778077", met_scores_batch_brain[[1]]), ]
#try top 5%
met_scores_batch_brain_chr5_73778077_sig <- met_scores_batch_brain_chr5_73778077 %>% 
  filter(abs (quantile_score) > 0.99 & gene_name == "ARHGEF28")


####################
#STR with suggestive significance and xQTL chr16:31328793

met_scores_batch_brain_chr16_31328793 <- met_scores_batch_brain[grepl("chr16:31328793", met_scores_batch_brain[[1]]), ]
#try top 1%
met_scores_batch_brain_chr16_31328793_sig <- met_scores_batch_brain_chr16_31328793 %>% 
  filter(abs (quantile_score) > 0.99 & gene_name == "ITGAM")

####################
#for sensory phenotype, the STR with xQTL annotation
sen_scores_batch_brain <- sen_scores_batch %>% 
  filter(ontology_curie == "UBERON:0000955")


sen_scores_batch_brain_chr14_92705503 <- sen_scores_batch_brain[grepl("chr14:92705503", sen_scores_batch_brain[[1]]), ]
#try top 1%
sen_scores_batch_brain_chr14_92705503_sig <- sen_scores_batch_brain_chr14_92705503 %>% 
  filter(abs (quantile_score) > 0.99 & gene_name == "LGMN")



#########

