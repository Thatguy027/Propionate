# generate data
# Rscript --vanilla Scripts/Interval_Popgen.R I 12274204 12488791 ../Ce330_annotated.vcf.gz ../WS245_exons.gff Propionate ../ce330_strains.txt

# negative fay an wu h = excess of high frequency derived alleles
# Rscript --vanilla Scripts/Interval_Popgen.R I 1 15072434 ../Ce330_annotated.vcf.gz ../WS245_exons.gff Propionate ../ce330_strains.txt

library(tidyverse)

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
setwd("..")

source("Scripts/Figure_Themes.R")

load("Propionate_I_1-15072434_Diversity_Statistics.Rda")

qtl_start <- 12374204
qtl_end <- 12388791
glct3 <- c(12385766, 12388791)
glct1 <- c(12337968, 12339239)
glct2 <- c(12435726, 12437285)
glct4 <- c(12317321, 12318735)
glct5 <- c(12453389, 12454764)

td_df <- d2_df %>%
  dplyr::filter(statistic %in% c("Fay.Wu.H")) %>%
  dplyr::group_by(statistic) %>%
  dplyr::mutate(scaled_value = scale(value)) %>%
  dplyr::mutate(q10 = quantile(value, 0.01, na.rm = T)) %>%
  dplyr::mutate(outlier = ifelse(value < q10, "Low", "not"))

chrom1 <- td_df %>%
  ggplot()+
  aes(x = startWindow/1e6, y = value, color = outlier)+
  geom_vline(aes(xintercept = glct3[1]/1e6), color = "#E68FAC")+
  geom_vline(aes(xintercept = glct3[2]/1e6), color = "#E68FAC")+
  geom_vline(aes(xintercept = glct1[1]/1e6), color = "gray80")+
  geom_vline(aes(xintercept = glct1[2]/1e6), color = "gray80")+
  geom_vline(aes(xintercept = glct2[1]/1e6), color = "gray60")+
  geom_vline(aes(xintercept = glct2[2]/1e6), color = "gray60")+
  geom_vline(aes(xintercept = glct4[1]/1e6), color = "gray40")+
  geom_vline(aes(xintercept = glct4[2]/1e6), color = "gray40")+
  geom_vline(aes(xintercept = glct5[1]/1e6), color = "gray20")+
  geom_vline(aes(xintercept = glct5[2]/1e6), color = "gray20")+
  geom_point(size = point_size, alpha = point_alpha)+
  scale_color_manual(values = c(highlight_color, "#222222"))+
  # facet_grid(statistic~., scales = "free")+
  # xlim(12,13)+
  theme_bw(18) +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  labs(x = "Genomic Position (Mb)",
       y = "Fay and Wu's H") 


region <-  td_df %>%
  ggplot()+
  aes(x = startWindow/1e6, y = value, color = outlier)+
  geom_vline(aes(xintercept = glct3[1]/1e6), color = "#E68FAC")+
  geom_vline(aes(xintercept = glct3[2]/1e6), color = "#E68FAC")+
  geom_vline(aes(xintercept = glct1[1]/1e6), color = "gray80")+
  geom_vline(aes(xintercept = glct1[2]/1e6), color = "gray80")+
  geom_vline(aes(xintercept = glct2[1]/1e6), color = "gray60")+
  geom_vline(aes(xintercept = glct2[2]/1e6), color = "gray60")+
  geom_vline(aes(xintercept = glct4[1]/1e6), color = "gray40")+
  geom_vline(aes(xintercept = glct4[2]/1e6), color = "gray40")+
  geom_vline(aes(xintercept = glct5[1]/1e6), color = "gray20")+
  geom_vline(aes(xintercept = glct5[2]/1e6), color = "gray20")+
  geom_point(size = point_size, alpha = point_alpha)+
  scale_color_manual(values = c(highlight_color, "#222222"))+
  # facet_grid(statistic~., scales = "free")+
  xlim(12,13)+
  theme_bw(18) +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  labs(x = "Genomic Position (Mb)",
       y = "Fay and Wu's H") 

cowplot::plot_grid(chrom1 + theme(axis.title.x = element_blank()),
                   region,
                   ncol =1, 
                   label_size = 20, 
                   align = "v", axis = "l",
                   labels = "AUTO")

ggsave(filename = "Plots/glct3_popgen.png", height = 8, width = 12, dpi = 400)
ggsave(filename = "Plots/glct3_popgen.pdf", height = 8, width = 12, dpi = 400)


td_df <- d2_df %>%
  dplyr::filter(statistic %in% c("nuc.diversity.within")) %>%
  dplyr::group_by(statistic) %>%
  dplyr::mutate(scaled_value = scale(value)) %>%
  dplyr::mutate(q10 = quantile(value, .99, na.rm = T)) %>%
  dplyr::mutate(outlier = ifelse(value > q10, "Low", "not"))

chrom1 <- td_df %>%
  ggplot()+
  aes(x = startWindow/1e6, y = value/10, color = outlier)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  annotation_logticks(sides = "l", 
                      color = "gray60", 
                      short = unit(0.05, "cm"), 
                      mid = unit(0.1, "cm"), 
                      long = unit(0.15, "cm"), 
                      size = 0.5)+
  geom_vline(aes(xintercept = glct3[1]/1e6), color = "#E68FAC")+
  geom_vline(aes(xintercept = glct3[2]/1e6), color = "#E68FAC")+
  geom_vline(aes(xintercept = glct1[1]/1e6), color = "gray80")+
  geom_vline(aes(xintercept = glct1[2]/1e6), color = "gray80")+
  geom_vline(aes(xintercept = glct2[1]/1e6), color = "gray60")+
  geom_vline(aes(xintercept = glct2[2]/1e6), color = "gray60")+
  geom_vline(aes(xintercept = glct4[1]/1e6), color = "gray40")+
  geom_vline(aes(xintercept = glct4[2]/1e6), color = "gray40")+
  geom_vline(aes(xintercept = glct5[1]/1e6), color = "gray20")+
  geom_vline(aes(xintercept = glct5[2]/1e6), color = "gray20")+
  geom_point(size = point_size, alpha = point_alpha)+
  scale_color_manual(values = c(highlight_color, "#222222"))+
  # facet_grid(statistic~., scales = "free")+
  # xlim(12,13)+
  theme_bw(18) +
  theme(axis.text.x = ggplot2::element_text(size = 12),
        axis.text.y = ggplot2::element_text(size = 12, hjust = 1),
        axis.title.x = ggplot2::element_text(size = 14, face = "bold", color = "black", vjust = -0.3),
        axis.title.y = ggplot2::element_text(size = 14, face = "bold", color = "black", vjust = -0.3),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.position = "none") +
  labs(x = "Genomic Position (Mb)",
       y = expression(pi~"(per kb)")) 

chrom1

region <- td_df %>%
  ggplot()+
  aes(x = startWindow/1e6, y = value/10, color = outlier)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  annotation_logticks(sides = "l", 
                      color = "gray60", 
                      short = unit(0.05, "cm"), 
                      mid = unit(0.1, "cm"), 
                      long = unit(0.15, "cm"), 
                      size = 0.5)+
  geom_vline(aes(xintercept = glct3[1]/1e6), color = "#E68FAC")+
  geom_vline(aes(xintercept = glct3[2]/1e6), color = "#E68FAC")+
  geom_vline(aes(xintercept = glct1[1]/1e6), color = "gray80")+
  geom_vline(aes(xintercept = glct1[2]/1e6), color = "gray80")+
  geom_vline(aes(xintercept = glct2[1]/1e6), color = "gray60")+
  geom_vline(aes(xintercept = glct2[2]/1e6), color = "gray60")+
  geom_vline(aes(xintercept = glct4[1]/1e6), color = "gray40")+
  geom_vline(aes(xintercept = glct4[2]/1e6), color = "gray40")+
  geom_vline(aes(xintercept = glct5[1]/1e6), color = "gray20")+
  geom_vline(aes(xintercept = glct5[2]/1e6), color = "gray20")+
  geom_point(size = point_size, alpha = point_alpha)+
  scale_color_manual(values = c(highlight_color, "#222222"))+
  # facet_grid(statistic~., scales = "free")+
  xlim(12,13)+
  theme_bw(18) +
  theme(axis.text.x = ggplot2::element_text(size = 12),
        axis.text.y = ggplot2::element_text(size = 12, hjust = 1),
        axis.title.x = ggplot2::element_text(size = 14, face = "bold", color = "black", vjust = -0.3),
        axis.title.y = ggplot2::element_text(size = 14, face = "bold", color = "black", vjust = -0.3),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 14, face = "bold"),
        legend.position = "none") +
  labs(x = "Genomic Position (Mb)",
       y = expression(pi~"(per kb)")) 

cowplot::plot_grid(chrom1 + theme(axis.title.x = element_blank()),
                   region,
                   ncol =1, 
                   label_size = 20, 
                   align = "v", axis = "l",
                   labels = "AUTO")

ggsave(filename = "Plots/glct3_popgen.png", height = 8, width = 12, dpi = 400)
ggsave(filename = "Plots/glct3_popgen.pdf", height = 8, width = 12, dpi = 400)
