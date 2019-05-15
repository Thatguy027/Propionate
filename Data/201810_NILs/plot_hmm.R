library(tidyverse)

setwd("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/201810_NILs/")

df <- readr::read_tsv("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/201810_NILs/gt_hmm_fill.tsv") %>%
  dplyr::group_by(sample) %>%
  dplyr::filter(chrom != "MtDNA") %>%
  dplyr::mutate(gt = ifelse(gt == 1, "DL238", "BRC20067")) %>%
  dplyr::mutate(low_sites = ifelse(sites < 100, TRUE, FALSE)) 

df$index <- dplyr::group_indices(df)

strain_index <- df$sample
names(strain_index) <- df$index + 0.5

ggplot(df,  aes(xmin = start, xmax = end, ymin = index, ymax = index + 1, fill = gt)) +
  geom_rect(aes(alpha = low_sites)) +
  scale_alpha_discrete(range = c(1.0, 0.65)) +
  scale_fill_manual(values = c("hotpink3", "cadetblue3")) +
  facet_grid(.~chrom, scales="free", space="free") +
  scale_x_continuous(labels = function(x) { x/1e6 }, expand = c(0,0)) +
  scale_y_continuous(breaks = unique(df$index) + 0.5, labels = function(x) { strain_index[as.character(x)] }, expand = c(0,0)) + 
  theme(strip.background = element_blank(),
        legend.position = "None")

ggsave("gt_hmm.png", height = 8, width = 12)

df%>%
  dplyr::filter(chrom =="V") %>%
  ggplot(.,  aes(xmin = start, xmax = end, ymin = index, ymax = index + 1, fill = gt)) +
  geom_rect(aes(alpha = low_sites)) +
  scale_alpha_discrete(range = c(1.0, 0.65)) +
  scale_fill_manual(values = c("hotpink3", "cadetblue3")) +
  facet_grid(.~chrom, scales="free", space="free") +
  scale_x_continuous(labels = function(x) { x/1e6 }, expand = c(0,0)) +
  scale_y_continuous(breaks = unique(df$index) + 0.5, labels = function(x) { strain_index[as.character(x)] }, expand = c(0,0)) + 
  theme(strip.background = element_blank(),
        legend.position = "None")

ggsave("gt_hmm_chrV.png", height = 8, width = 12)

df%>%
  dplyr::filter(chrom =="V", !sample %in% c("DL238", "BRC20067")) %>%
  ggplot(.,  aes(xmin = start, xmax = end, ymin = index, ymax = index + 1, fill = gt)) +
  geom_rect(aes(alpha = low_sites)) +
  scale_alpha_discrete(range = c(1.0, 0.65)) +
  scale_fill_manual(values = c("hotpink3", "cadetblue3")) +
  facet_grid(.~chrom, scales="free", space="free") +
  scale_x_continuous(labels = function(x) { x/1e6 }, expand = c(0,0)) +
  scale_y_continuous(breaks = unique(df$index) + 0.5, labels = function(x) { strain_index[as.character(x)] }, expand = c(0,0)) + 
  theme(strip.background = element_blank(),
        legend.position = "None")+xlim(c(3e6,5e6))

ggsave("gt_hmm_chrV_zoom.png", height = 8, width = 12)
