library(data.table)
library(tidyverse)
library(cegwas)
library(cegwas2)
library(broom)
library(cowplot)

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
setwd("..")

source("Scripts/Figure_Functions.R")

# figure 2
c2_prmaps <- data.table::fread("Processed_Data/GWAS_processed_mapping_cegwas2.tsv")
pr_resid <- data.table::fread("Processed_Data/GWAS_PRpheno_replicates.tsv")

pr_resid$strain <- factor(pr_resid$strain,levels = pr_resid$strain,labels = pr_resid$strain, ordered = T)

pr_resid_pldf <- pr_resid %>%
  dplyr::ungroup() %>%
  dplyr::arrange(resid) %>%
  dplyr::mutate(norm_pheno_temp = ifelse(resid == min(resid), 0, 1))%>%
  dplyr::mutate(delta_pheno = ifelse(norm_pheno_temp == 0, 0, abs(dplyr::lag(resid) - resid)))%>%
  dplyr::mutate(norm_pheno = cumsum(delta_pheno)) %>%
  dplyr::mutate(final_pheno = norm_pheno/max(norm_pheno)) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(mph = mean(final_pheno, na.rm = T),
                sph = sd(final_pheno, na.rm = T)) %>%
  dplyr::distinct(strain, mph, sph,.keep_all=T) %>%
  dplyr::ungroup()%>%
  dplyr::mutate(st_colors = ifelse(strain == "N2", "N2",
                                   ifelse(strain == "BRC20067", "BRC20067",
                                          ifelse(strain == "DL238", "DL238", "not"))))

pheno_plot <-ggplot(pr_resid_pldf) +
  aes(y = mph, x = strain, fill = st_colors) +
  geom_bar(stat="identity", color="black", position=position_dodge(), size = 0.25) +
  geom_errorbar(aes(ymin=mph, ymax=mph+sph), width=.4,
                position=position_dodge(.9)) +
  scale_fill_manual(values = c("DL238"  = "cadetblue3", "N2"  = "orange", "BRC20067" = "hotpink3", "not" = "gray70")) +
  theme_classic(12) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank(), axis.line = element_line(size = 0.25)) +
  geom_rect(aes(ymin = 0.75, ymax = 0.75+0.015*2, xmin = "N2", xmax = "LSJ1"), fill = "orange", color = "black")+
  geom_rect(aes(ymin = 0.75+0.015*5, ymax = 0.75+0.015*7, xmin = "N2", xmax = "LSJ1"), fill = "hotpink3", color = "black")+
  geom_rect(aes(ymin = 0.75+0.015*10, ymax = 0.75+0.015*12, xmin = "N2", xmax = "LSJ1"), fill = "cadetblue3", color = "black")+
  geom_text(aes(x="JU1212", y=(0.75+0.75+0.015*2)/2, label= "N2"), size=4, hjust = 0, check_overlap = TRUE) +
  geom_text(aes(x="JU1212", y=((0.75+0.015*5)+( 0.75+0.015*7))/2, label= "BRC20067"), size=4, hjust = 0, check_overlap = TRUE) +
  geom_text(aes(x="JU1212", y=((0.75+0.015*10)+( 0.75+0.015*12))/2, label= "DL238"), size=4, hjust = 0, check_overlap = TRUE)+
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  labs(x = "Strain", y = "Normalized\nL1 survival")


rrblup_map <- c2_prmaps %>%
  dplyr::filter(CHROM != "MtDNA") %>%
  dplyr::distinct(marker, .keep_all = T) %>%
  dplyr::mutate(EIGEN_CUTOFF = -log10(.05/772)) %>%
  dplyr::mutate(EIGEN_SIG = ifelse(log10p > BF, "1", 
                                   ifelse(log10p > EIGEN_CUTOFF, "2", "0")) )%>%
  ggplot2::ggplot(.) +
  ggplot2::aes(x = POS/1e6, y = log10p) +
  ggplot2::scale_color_manual(values = c("0" = "black", 
                                         "1" = "red",
                                         "2" = "hotpink3")) +
  ggplot2::scale_alpha_manual(values = c("0" = 0.5, 
                                         "1" = 1,
                                         "2" = 1)) +
  ggplot2::geom_rect(ggplot2::aes(xmin = startPOS/1e6, 
                                  xmax = endPOS/1e6, 
                                  ymin = 0, 
                                  ymax = Inf, 
                                  fill = "blue"), 
                     color = "blue",fill = "cyan",
                     linetype = 2,
                     size = 0.25,
                     alpha=.3)+
  ggplot2::geom_hline(ggplot2::aes(yintercept = BF),
                      color = "gray60", 
                      alpha = .75,  
                      size = 0.5) +
  ggplot2::geom_hline(ggplot2::aes(yintercept = EIGEN_CUTOFF),
                      color = "gray60", 
                      alpha = .75,  
                      size = 0.5,
                      linetype = 2) +
  ggplot2::geom_point( size = 0.5, ggplot2::aes(color= factor(EIGEN_SIG), alpha = factor(EIGEN_SIG)) ) +
  ggplot2::facet_grid( . ~ CHROM, scales = "free_x" , space = "free_x") +
  ggplot2::theme_bw(12) +
  scale_y_continuous(limits = c(0,10), expand = c(0, 0)) +
  ggplot2::theme(strip.background = element_blank(),
                 legend.position = "none",
                 panel.grid.minor = element_blank(),
                 panel.grid.major = element_blank()) +
  ggplot2::labs(x = "Genomic position (Mb)",
                y = expression(-log[10](italic(p))))

skat_maps <- data.table::fread("Processed_Data/GWAS_Skat.assoc") %>%
  dplyr::rowwise()%>%
  dplyr::mutate(CHROM = strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][1],
                POS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][1]),
                endPOS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][2]))%>%
  dplyr::mutate(size = abs(POS-endPOS)) %>%
  dplyr::filter(CHROM!="MtDNA")%>%
  dplyr::ungroup()%>%
  dplyr::filter(NumVar > 1)%>%
  dplyr::mutate(significant = ifelse(Pvalue < .05/n(), TRUE,FALSE ))%>%
  dplyr::arrange(Pvalue)

skat_plot <- skat_maps%>%
  ggplot()+
  aes(x = POS/1e6, y = -log10(Pvalue), alpha = 0.5, color = significant)+
  geom_point(size = 0.5)+
  scale_color_manual(values=c("black","red"))+
  facet_grid(.~CHROM, scales = "free", space = "free")+
  ggplot2::theme_bw(12) +
  ggplot2::theme(strip.background = element_blank(),
                 legend.position = "none",
                 panel.grid.minor = element_blank(),
                 panel.grid.major = element_blank()) +
  ggplot2::geom_hline(ggplot2::aes(yintercept = -log10(0.05/nrow(skat_maps))),
                      color = "gray60", 
                      alpha = .75,  
                      size = 0.5) +
  scale_y_continuous(limits = c(0,10), expand = c(0, 0)) +
  labs(x = "Genomic position (Mb)", 
       y = expression(-log[10](italic(p))))

# figure 2 cowplot
figure2<-cowplot::plot_grid(pheno_plot,
                            rrblup_map+ theme(axis.text.x = element_blank(),
                                              axis.title.x = element_blank()),
                            skat_plot + theme(strip.text.x = element_blank()),
                            ncol =1, 
                            label_size = 14, 
                            rel_heights = c(0.8,1,1),
                            align = "v", 
                            axis = "l",
                            labels = "AUTO")


ggsave(filename = "Plots/SVG_PLOTS/Figure2_up.svg", height = 8, width = 6.5, units = "in")
ggsave(filename = "Plots/SVG_PLOTS/Figure2_up.pdf", height = 8, width = 6.5, units = "in")
ggsave(filename = "Plots/SVG_PLOTS/Figure2_up.png", height = 8, width = 6.5, units = "in")
