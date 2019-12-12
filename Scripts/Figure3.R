library(data.table)
library(tidyverse)
library(cegwas)
library(cegwas2)
library(broom)
library(cowplot)

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
setwd("..")


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

# phenotypes

prop_pheno <- data.table::fread("Data/20180625_nil_phenotypes/20180625_nil_phenos.csv")%>%
  na.omit()
prop_pheno$strain <- gsub("-","_",prop_pheno$strain)
prop_pheno$day <- as.factor(prop_pheno$day)

prop_pheno_n4 <- prop_pheno %>%
  dplyr::filter(day!=4)

prop_pheno_n4$day_regressed <- residuals(glm(survival ~ day, data = prop_pheno_n4))

strain_order <- prop_pheno_n4 %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(m_p = mean(day_regressed))%>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(m_p)) %>%
  dplyr::distinct(strain)


dplyr::ungroup() %>%
  dplyr::arrange(resid) %>%
  dplyr::mutate(norm_pheno_temp = ifelse(resid == min(resid), 0, 1))%>%
  dplyr::mutate(delta_pheno = ifelse(norm_pheno_temp == 0, 0, abs(dplyr::lag(resid) - resid)))%>%
  dplyr::mutate(norm_pheno = cumsum(delta_pheno)) %>%
  dplyr::mutate(final_pheno = norm_pheno/max(norm_pheno))

tidy_pheno <- prop_pheno_n4 %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(m_p = mean(day_regressed))%>%
  dplyr::ungroup() %>%
  dplyr::filter(strain != "A6_R2") %>%
  dplyr::mutate(strain = gsub("_R1","",strain)) %>%
  dplyr::mutate(nil_name = factor(strain, levels = c("DL238",
                                                     "B3", "B6", "B9","B2","A3","B10",
                                                     "B8","C6","A6","BRC20067"),
                                  labels = c("DL238",
                                             "B3", "B6", "B9","B2","A3","B10",
                                             "B8","C6","A6","BRC20067"))) %>%
  dplyr::left_join(.,nil_name_conv,by="nil_name") %>%
  dplyr::ungroup() %>%
  dplyr::arrange(day_regressed) %>%
  dplyr::mutate(norm_pheno_temp = ifelse(day_regressed == min(day_regressed), 0, 1))%>%
  dplyr::mutate(delta_pheno = ifelse(norm_pheno_temp == 0, 0, abs(dplyr::lag(day_regressed) - day_regressed)))%>%
  dplyr::mutate(norm_pheno = cumsum(delta_pheno)) %>%
  dplyr::mutate(final_pheno = norm_pheno/max(norm_pheno)) %>%
  dplyr::arrange(desc(m_p)) %>%
  dplyr::mutate(final_name = factor(new_name, levels = unique(new_name), labels = unique(new_name)))

g_pt <- tidy_pheno %>%
  dplyr::filter(chrom =="V") %>%
  ggplot(.,  aes(x = start/1e6, xend = end/1e6, y = final_name, yend = final_name, color = gt)) +
  geom_segment( size = 3) +
  geom_vline(aes(xintercept=3.213649), linetype = 2, alpha = 0.5)+
  geom_vline(aes(xintercept=4.284434), linetype = 2, alpha = 0.5)+
  scale_color_manual(values = c("hotpink3", "cadetblue3")) +
  theme_classic(12) +
  theme(strip.background = element_blank(),
        legend.position = "None")+
  coord_cartesian(xlim=c(3, 5)) +
  labs(x = "Genomic position (Mb)", y = "Strain")


p_pt <- tidy_pheno %>%
  dplyr::distinct(final_name, final_pheno) %>%
  ggplot()+
  aes(x = final_name, y = final_pheno, fill = final_name)+
  ggbeeswarm::geom_beeswarm(cex = 0.7, size =0.5, priority = "density")+
  geom_boxplot(alpha = 0.8, outlier.colour = NA)+
  theme_classic(12)+
  scale_fill_manual(values=c("BRC20067"="hotpink3", "DL238"  = "cadetblue3",rep("gray60",10)))+
  coord_flip()+
  labs(y = "Normalized L1 survival") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        legend.position = "none")

cowplot::plot_grid(g_pt, 
                   p_pt,
                   
                   # align = "h", 
                   labels = "AUTO",
                   # rel_widths = c(.5,1),
                   label_x = c(0,-.05),
                   label_size = 14, scale = 0.95)

ggsave("Plots/SVG_PLOTS/Figure3.pdf", height = 3.5, width = 6.5, units = "in")
ggsave("Plots/SVG_PLOTS/Figure3.png", height = 3.5, width = 6.5, units = "in", dpi = 300)
ggsave("Plots/SVG_PLOTS/Figure3.svg", height = 3.5, width = 6.5, units = "in")
