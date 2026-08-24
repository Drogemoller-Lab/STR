library(dplyr)
library(tidyr)
library(purrr)


comb_met_assoctiaTR <- read_tsv ("/home/projects/hearing_loss/clsaARHL_SA/str/results/combined_met_association_results.tsv")
comb_sen_assoctiaTR <- read_tsv ("/home/projects/hearing_loss/clsaARHL_SA/str/results/combined_sen_association_results.tsv")

#filter suggestive significant variants
met_str_sig <- filter (comb_met_assoctiaTR, p_met_phenotype < 1.49e-7)
met_str_sig_sugest <- filter (comb_met_assoctiaTR, p_met_phenotype < 1e-5)

#identify the number of loci
met_str_sig_sugest_group <- met_str_sig_sugest %>% 
  group_by(chrom) |> 
  arrange(p_met_phenotype) |>
  slice(1) |>  
  ungroup()

met_variants <- met_str_sig_sugest %>% 
  select(chrom,pos) %>% 
  rename (CHR=chrom,
          POS=pos)

# ####functional annotation of the STR uncovered by SuSiE
# "EnsTR:chr5:73778077"
# 
# #identify the number of loci
# chr5_73778077 <- comb_met_assoctiaTR %>% 
#   filter(chrom == "chr5" & pos == 73778077)
# 
# met_variants <- chr5_73778077 %>% 
#   select(chrom,pos) %>% 
#   rename (CHR=chrom,
#           POS=pos)


#filter suggestive significant variants
sen_str_sig <- filter (comb_sen_assoctiaTR, p_sen_phenotype < 1.49e-7)
sen_str_sig_sugest <- filter (comb_sen_assoctiaTR, p_sen_phenotype < 1e-5)
#identify the number of loci
sen_str_sig_sugest_group <- sen_str_sig_sugest %>% 
  group_by(chrom) |> 
  arrange(p_sen_phenotype) |>
  slice(1) |>  
  ungroup()

sen_variants <- sen_str_sig_sugest %>% 
  select(chrom,pos)%>% 
  rename (CHR=chrom,
          POS=pos)

#####

#read TR-xQTL files
#splice files
temporal_s <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/NYGCALS.Cortex_Temporal.TR-sQTL.txt")
frontal_s <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/NYGCALS.Cortex_Frontal.TR-sQTL.txt")
motor_s <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/NYGCALS.Cortex_Motor.TR-sQTL.txt")
cerebellum_s <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/NYGCALS.Cerebellum.TR-sQTL.txt")
cervical_s <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/NYGCALS.Spinal_Cord_Cervical.TR-sQTL.txt")
lumbar_s <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/NYGCALS.Spinal_Cord_Lumbar.TR-sQTL.txt")
dlpfc_s <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/ROSMAP.DLPFC.TR-sQTL.txt")
ipsc_s <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/AnswerALS.iPSC-MN.TR-sQTL.txt")

######

#####
# V2

find_matches_variant_df <- function(df, df_name, variants_df) {
  df %>%
    separate(V2, into = c("CHR", "start", "end", "str"), sep = "_", fill = "right") %>%
    mutate(start = as.numeric(start), end = as.numeric(end)) %>%
    right_join(variants_df, by = "CHR", relationship = "many-to-many") %>%
    filter(POS >= start, POS <= end) %>%
    mutate(dataset = df_name)
}



met_all_matches <- bind_rows(
  find_matches_variant_df(frontal_s, "frontal", met_variants),
  find_matches_variant_df(motor_s, "motor", met_variants),
  find_matches_variant_df(cerebellum_s, "cerebellum", met_variants),
  find_matches_variant_df(cervical_s, "cervical", met_variants),
  find_matches_variant_df(lumbar_s, "lumbar", met_variants),
  find_matches_variant_df(dlpfc_s, "dlpfc", met_variants),
  find_matches_variant_df(ipsc_s, "ipsc", met_variants),
  find_matches_variant_df(temporal_s, "temporal", met_variants)
)
print(met_all_matches)

#####
#for the sensory phenotype
sen_all_matches <- bind_rows(
  find_matches_variant_df(frontal_s, "frontal", sen_variants),
  find_matches_variant_df(motor_s, "motor", sen_variants),
  find_matches_variant_df(cerebellum_s, "cerebellum", sen_variants),
  find_matches_variant_df(cervical_s, "cervical", sen_variants),
  find_matches_variant_df(lumbar_s, "lumbar", sen_variants),
  find_matches_variant_df(dlpfc_s, "dlpfc", sen_variants),
  find_matches_variant_df(ipsc_s, "ipsc", sen_variants),
  find_matches_variant_df(temporal_s, "temporal", sen_variants)
)
print(sen_all_matches)


#try eQTL

#expression files
temporal_e <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/NYGCALS.Cortex_Temporal.TR-eQTL.txt")
frontal_e <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/NYGCALS.Cortex_Frontal.TR-eQTL.txt")
motor_e <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/NYGCALS.Cortex_Motor.TR-eQTL.txt")
cerebellum_e <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/NYGCALS.Cerebellum.TR-eQTL.txt")
cervical_e <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/NYGCALS.Spinal_Cord_Cervical.TR-eQTL.txt")
lumbar_e <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/NYGCALS.Spinal_Cord_Lumbar.TR-eQTL.txt")
dlpfc_e <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/ROSMAP.DLPFC.TR-eQTL.txt")
ipsc_e <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/AnswerALS.iPSC-MN.TR-eQTL.txt")

met_all_matches_e <- bind_rows(
  find_matches_variant_df(frontal_e, "frontal", met_variants),
  find_matches_variant_df(motor_e, "motor", met_variants),
  find_matches_variant_df(cerebellum_e, "cerebellum", met_variants),
  find_matches_variant_df(cervical_e, "cervical", met_variants),
  find_matches_variant_df(lumbar_e, "lumbar", met_variants),
  find_matches_variant_df(dlpfc_e, "dlpfc", met_variants),
  find_matches_variant_df(ipsc_e, "ipsc", met_variants),
  find_matches_variant_df(temporal_e, "temporal", met_variants)
)
print(met_all_matches_e)

sen_all_matches_e <- bind_rows(
    find_matches_variant_df(frontal_e, "frontal", sen_variants),
    find_matches_variant_df(motor_e, "motor", sen_variants),
    find_matches_variant_df(cerebellum_e, "cerebellum", sen_variants),
    find_matches_variant_df(cervical_e, "cervical", sen_variants),
    find_matches_variant_df(lumbar_e, "lumbar", sen_variants),
    find_matches_variant_df(dlpfc_e, "dlpfc", sen_variants),
    find_matches_variant_df(ipsc_e, "ipsc", sen_variants),
    find_matches_variant_df(temporal_e, "temporal", sen_variants)
)
print(sen_all_matches_e)

#####
#for 3' APA
temporal_3a <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/NYGCALS.Cortex_Temporal.TR-3aQTL.txt")
frontal_3a <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/NYGCALS.Cortex_Frontal.TR-3aQTL.txt")
motor_3a <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/NYGCALS.Cortex_Motor.TR-3aQTL.txt")
cerebellum_3a <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/NYGCALS.Cerebellum.TR-3aQTL.txt")
cervical_3a <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/NYGCALS.Spinal_Cord_Cervical.TR-3aQTL.txt")
lumbar_3a <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/NYGCALS.Spinal_Cord_Lumbar.TR-3aQTL.txt")
dlpfc_3a <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/ROSMAP.DLPFC.TR-3aQTL.txt")
ipsc_3a <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/AnswerALS.iPSC-MN.TR-3aQTL.txt")
  

met_all_matches_3a <- bind_rows(
  find_matches_variant_df(frontal_3a, "frontal", met_variants),
  find_matches_variant_df(motor_3a, "motor", met_variants),
  find_matches_variant_df(cerebellum_3a, "cerebellum", met_variants),
  find_matches_variant_df(cervical_3a, "cervical", met_variants),
  find_matches_variant_df(lumbar_3a, "lumbar", met_variants),
  find_matches_variant_df(dlpfc_3a, "dlpfc", met_variants),
  find_matches_variant_df(ipsc_3a, "ipsc", met_variants),
  find_matches_variant_df(temporal_3a, "temporal", met_variants)
)
print(met_all_matches_3a)

sen_all_matches_3a <- bind_rows(
  find_matches_variant_df(frontal_3a, "frontal", sen_variants),
  find_matches_variant_df(motor_3a, "motor", sen_variants),
  find_matches_variant_df(cerebellum_3a, "cerebellum", sen_variants),
  find_matches_variant_df(cervical_3a, "cervical", sen_variants),
  find_matches_variant_df(lumbar_3a, "lumbar", sen_variants),
  find_matches_variant_df(dlpfc_3a, "dlpfc", sen_variants),
  find_matches_variant_df(ipsc_3a, "ipsc", sen_variants),
  find_matches_variant_df(temporal_3a, "temporal", sen_variants)
)
print(sen_all_matches_3a)


###
#for methylation QTL
dlpfc_m <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/ROSMAP.DLPFC.TR-mQTL.txt")

met_all_matches_m <- bind_rows(find_matches_variant_df(dlpfc_m, "dlpfc", met_variants))
print(met_all_matches_m)

sen_all_matches_m <- bind_rows(find_matches_variant_df(dlpfc_m, "dlpfc", sen_variants))
print(sen_all_matches_m)

#for ATAC-seq
ipsc_ca <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/AnswerALS.iPSC-MN.TR-caQTL.txt")
met_all_matches_ca <- bind_rows(find_matches_variant_df(ipsc_ca, "ipsc", met_variants))
print(met_all_matches_ca)

sen_all_matches_ca <- bind_rows(find_matches_variant_df(ipsc_ca, "ipsc", sen_variants))
print(sen_all_matches_ca)

#for hQTL
dlpfc_h <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/ROSMAP.DLPFC.TR-hQTL.txt")

met_all_matches_h <- bind_rows(find_matches_variant_df(dlpfc_h, "dlpfc", met_variants))
print(met_all_matches_h)

sen_all_matches_h <- bind_rows(find_matches_variant_df(dlpfc_h, "dlpfc", sen_variants))
print(sen_all_matches_h)

#for protein QTL
ipsc_p <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/AnswerALS.iPSC-MN.TR-pQTL.txt")
met_all_matches_p <- bind_rows(find_matches_variant_df(ipsc_p, "ipsc", met_variants))
print(met_all_matches_p)

sen_all_matches_p <- bind_rows(find_matches_variant_df(ipsc_p, "ipsc", sen_variants))
print(sen_all_matches_p)



#####
#Try GTEx
gtex_amygdala <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/Brain_Amygdala.TR-eQTL.txt.gz")
gtex_anterior_cortex <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/Brain_Anterior_cingulate_cortex_BA24.TR-eQTL.txt.gz")
gtex_caudate <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/Brain_Caudate_basal_ganglia.TR-eQTL.txt.gz")
gtex_cerebellar_hemisphere <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/Brain_Cerebellar_Hemisphere.TR-eQTL.txt.gz")
gtex_cerebellum <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/Brain_Cerebellum.TR-eQTL.txt.gz")
gtex_cortex <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/Brain_Cortex.TR-eQTL.txt.gz")
gtex_hippocampus <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/Brain_Hippocampus.TR-eQTL.txt.gz")
gtex_nucleus <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/Brain_Nucleus_accumbens_basal_ganglia.TR-eQTL.txt.gz")
gtex_putamen <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/Brain_Putamen_basal_ganglia.TR-eQTL.txt.gz")
gtex_spinal_cord <- read.table("/home/projects/hearing_loss/clsaARHL_SA/str/tr_xqtls/Brain_Spinal_cord_cervical_c-1.TR-eQTL.txt.gz")


# Create a named list of all GTEx brain region datasets
brain_datasets <- list(
  gtex_amygdala = gtex_amygdala,
  gtex_anterior_cortex = gtex_anterior_cortex,
  gtex_caudate = gtex_caudate,
  gtex_cerebellar_hemisphere = gtex_cerebellar_hemisphere,
  gtex_cerebellum = gtex_cerebellum,
  gtex_cortex = gtex_cortex,
  gtex_hippocampus = gtex_hippocampus,
  gtex_nucleus = gtex_nucleus,
  gtex_putamen = gtex_putamen,
  gtex_spinal_cord = gtex_spinal_cord
)

# Apply your function across all datasets and combine results
met_all_matches_gtex <- map_dfr(
  names(brain_datasets),
  ~ find_matches_variant_df(brain_datasets[[.x]], .x, met_variants)
)

sen_all_matches_gtex <- map_dfr(
  names(brain_datasets),
  ~ find_matches_variant_df(brain_datasets[[.x]], .x, sen_variants)
)



#try the same function but with locus instead of the repeat position

# V1

find_matches_variant_region <- function(df, df_name, variants_df) {
  df %>%
    separate(V1, into = c("CHR", "start", "end", "str"), sep = ":", fill = "right") %>%
    mutate(start = as.numeric(start), end = as.numeric(end)) %>%
    right_join(variants_df, by = "CHR", relationship = "many-to-many") %>%
    filter(POS >= start, POS <= end) %>%
    mutate(dataset = df_name)
}

met_all_matches <- bind_rows(
  find_matches_variant_region(frontal_s, "frontal", met_variants),
  find_matches_variant_region(motor_s, "motor", met_variants),
  find_matches_variant_region(cerebellum_s, "cerebellum", met_variants),
  find_matches_variant_region(cervical_s, "cervical", met_variants),
  find_matches_variant_region(lumbar_s, "lumbar", met_variants),
  find_matches_variant_region(dlpfc_s, "dlpfc", met_variants),
  find_matches_variant_region(ipsc_s, "ipsc", met_variants),
  find_matches_variant_region(temporal_s, "temporal", met_variants)
)
print(met_all_matches)

met_all_matches_e <- bind_rows(
  find_matches_variant_region(frontal_e, "frontal", met_variants),
  find_matches_variant_region(motor_e, "motor", met_variants),
  find_matches_variant_region(cerebellum_e, "cerebellum", met_variants),
  find_matches_variant_region(cervical_e, "cervical", met_variants),
  find_matches_variant_region(lumbar_e, "lumbar", met_variants),
  find_matches_variant_region(dlpfc_e, "dlpfc", met_variants),
  find_matches_variant_region(ipsc_e, "ipsc", met_variants),
  find_matches_variant_region(temporal_e, "temporal", met_variants)
)
print(met_all_matches_e)


met_all_matches_3a <- bind_rows(
  find_matches_variant_region(frontal_3a, "frontal", met_variants),
  find_matches_variant_region(motor_3a, "motor", met_variants),
  find_matches_variant_region(cerebellum_3a, "cerebellum", met_variants),
  find_matches_variant_region(cervical_3a, "cervical", met_variants),
  find_matches_variant_region(lumbar_3a, "lumbar", met_variants),
  find_matches_variant_region(dlpfc_3a, "dlpfc", met_variants),
  find_matches_variant_region(ipsc_3a, "ipsc", met_variants),
  find_matches_variant_region(temporal_3a, "temporal", met_variants)
)
print(met_all_matches_3a)


###
#for methylation QTL
met_all_matches_m <- bind_rows(find_matches_variant_region(dlpfc_m, "dlpfc", met_variants))
print(met_all_matches_m)

#for ATAC-seq
met_all_matches_ca <- bind_rows(find_matches_variant_region(ipsc_ca, "ipsc", met_variants))
print(met_all_matches_ca)

#for hQTL
met_all_matches_h <- bind_rows(find_matches_variant_region(dlpfc_h, "dlpfc", met_variants))
print(met_all_matches_h)

#for protein QTL
met_all_matches_p <- bind_rows(find_matches_variant_region(ipsc_p, "ipsc", met_variants))
print(met_all_matches_p)


