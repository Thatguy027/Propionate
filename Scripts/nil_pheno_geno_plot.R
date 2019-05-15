library(tidyverse)

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
  dplyr::left_join(.,df2,by="nil_name") %>%
  dplyr::arrange(desc(m_p)) %>%
  dplyr::mutate(nil_name1 = factor(nil_name, levels = c("DL238",
                                                        "B3", "B6", "B9","B2","A3","B10",
                                                        "B8","C6","A6","BRC20067"),
                                   labels = c("DL238",
                                              "B3", "B6", "B9","B2","A3","B10",
                                              "B8","C6","A6","BRC20067")))

g_pt <- tidy_pheno %>%
  dplyr::filter(chrom =="V") %>%
  ggplot(.,  aes(x = start/1e6, xend = end/1e6, y = nil_name1, yend = nil_name1, color = gt)) +
  geom_segment( size = 3) +
  scale_color_manual(values = c("hotpink3", "cadetblue3")) +
  theme_classic(20) +
  theme(strip.background = element_blank(),
        legend.position = "None")+
  # xlim(c(3,5))+
  labs(x = "Genomic Position (Mb)", y = "Strain")

p_pt <- tidy_pheno %>%
  dplyr::distinct(nil_name1, day_regressed) %>%
  ggplot()+
  aes(x = nil_name1, y = day_regressed)+
  geom_hline(yintercept = -3.66, color = "red", linetype = 2)+
  ggbeeswarm::geom_beeswarm(cex = 0.7, priority = "density")+
  geom_boxplot(alpha = 0.8)+
  theme_classic(20)+
  coord_flip()+
  labs(y = "L1 Survival") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank())

cowplot::plot_grid(g_pt, 
                   p_pt,
                   
                   # align = "h", 
                   labels = "AUTO",
                   # rel_widths = c(.5,1),
                   label_x = c(0,-.05),
                   label_size = 20)

ggsave("Plots/NIL_pheno_geno_plot.pdf", height = 8, width = 12)
ggsave("Plots/NIL_pheno_geno_plot.png", height = 8, width = 12, dpi = 300)
