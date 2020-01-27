library(data.table)
library(tidyverse)
library(cegwas)
library(cegwas2)
library(broom)
library(cowplot)

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
setwd("..")


# # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # # # # # # # # # # COMPLETE DOSE RESPONSES
raw_df <- fread("Data/Heritability/PA_doseresponse_herit.csv")

long_df <- raw_df %>%
  gather(strain_n, ratio, -DATE, -PA_conc) 

complete_dr <- data.frame(date = long_df$DATE, pa_conc = long_df$PA_conc, strain = long_df$strain_n, pheno = long_df$ratio)%>%
  separate(strain, into = c("strain", "replicate"), sep = "_")

write.table(complete_dr, file = "Final_Tables/Supplemental_Table_COARSE-DOSE-RESPONSE-PHENOTYPES.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

# # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # # # # # # # # # # FINE_SCALE DOSE RESPONSES
pheno <- fread("Data/20150219_heritability_v2/cleaned_phenotypes.csv") %>%
  data.frame()%>%
  gather(strain, pheno, -Dose, -day)

pheno$strain <- gsub(pattern = "strain.", replacement = "", pheno$strain)
pheno$strain <- gsub(pattern = "\\.1$|\\.2$", replacement = "", pheno$strain)

# # Fix strain names

pheno2 <- fread("Data/20150219_heritability_v2/strain_lookup.csv")
colnames(pheno2) <- c("strain", "strain_name")

pheno$strain <- as.numeric(pheno$strain)

pheno3 <- pheno %>%
  separate(Dose, into = c("Concentration", "value"), sep = " ",remove = FALSE)%>%
  mutate(n_conc = as.numeric(Concentration))%>%
  arrange(n_conc) %>%
  left_join(.,pheno2, by="strain") %>%
  dplyr::select(-Concentration, -value, -strain, -n_conc)

write.table(pheno3, file = "Final_Tables/Supplemental_Table_FINESCALE-DOSE-RESPONSE-PHENOTYPES.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

# # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # # # # # # # # # # GWAS PHENOTYPES

p <- data.table::fread("Data/20161128_new_gwas_replicates/20170103_set1-3.csv")%>%
  dplyr::select(-number, -average)%>%
  tidyr::gather(t_replicate, ratio, -strain, -replicate,-set) %>%
  dplyr::filter(strain!="") %>%
  dplyr::arrange(strain)

write.table(p, file = "Final_Tables/Supplemental_Table_GWAS_PHENOTYPES.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

# # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # # # # # # # # # # GWAS GENOTYPES

gm <- readr::read_tsv(glue::glue("Processed_Data/Genotype_Matrix.tsv"))

write_tsv(gm, "Final_Tables/Supplemental_Table_GWAS_GENOTYPES.tsv.gz", col_names = T)

# # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # # # # # # # # # # NIL GENOTYPES
df1 <- readr::read_tsv("Data/201810_NILs/gt_hmm_fill.tsv") %>%
  dplyr::group_by(sample) %>%
  dplyr::filter(chrom != "MtDNA") %>%
  dplyr::mutate(gt = ifelse(gt == 1, "DL238", "BRC20067")) %>%
  dplyr::mutate(low_sites = ifelse(sites < 100, TRUE, FALSE)) %>%
  dplyr::mutate(nil_name = ifelse(grepl("_",sample), strsplit(sample,split = "_")[[1]][2], sample))

df2 <- readr::read_tsv("Data/NIL_genotypes/gt_hmm_fill.tsv") %>%
  dplyr::group_by(sample) %>%
  dplyr::filter(chrom != "MtDNA") %>%
  dplyr::mutate(gt = ifelse(gt == 1, "DL238", "BRC20067")) %>%
  dplyr::mutate(low_sites = ifelse(sites < 100, TRUE, FALSE)) %>%
  dplyr::mutate(nil_name = ifelse(grepl("_",sample), strsplit(sample,split = "_")[[1]][2], sample))%>%
  dplyr::bind_rows(.,df1)%>%
  dplyr::filter(chrom == "V")%>%
  dplyr::distinct(sample,start,end,.keep_all=T) %>%
  dplyr::mutate(nil_name1 = factor(nil_name, levels = c("DL238",
                                                        "B3", "B6", "B9","B2","A3","B10",
                                                        "B8","C6","A6","BRC20067"),
                                   labels = c("DL238",
                                              "B3", "B6", "B9","B2","A3","B10",
                                              "B8","C6","A6","BRC20067")))


nil_name_conv <- readr::read_csv("Data/NIL_name_conversion.csv") %>%
  dplyr::rename(nil_name=Strain) %>%
  dplyr::left_join(df2,., by = "nil_name") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(new_name = ifelse(is.na(new_name), nil_name, new_name))


write.table(nil_name_conv, file = "Final_Tables/Supplemental_Table_NIL-GENOTYPES.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

# # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # # # # # # # # # # NIL PHENOTYPES

prop_pheno <- data.table::fread("Data/20180625_nil_phenotypes/20180625_nil_phenos.csv")%>%
  na.omit()
prop_pheno$strain <- gsub("-","_",prop_pheno$strain)
prop_pheno$day <- as.factor(prop_pheno$day)

prop_pheno_n4 <- prop_pheno %>%
  dplyr::filter(day!=4)

write.table(prop_pheno_n4, file = "Final_Tables/Supplemental_Table_NIL-PHENOTYPES.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

# # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # ## # # # # # # # # # # # # # # # CRISPR PHENOTYPES

phen <- data.table::fread("Data/glct_causality/glct_deletion_expt_full.csv")

write.table(phen, file = "Final_Tables/Supplemental_Table_CRISPR-PHENOTYPES.tsv", sep = "\t", quote = F, row.names = F, col.names = T)
