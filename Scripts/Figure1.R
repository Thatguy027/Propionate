library(data.table)
library(tidyverse)
library(cegwas)
library(cegwas2)
library(broom)
library(cowplot)
library(ggpubr)

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
setwd("..")

source("Scripts/Figure_Functions.R")

# figure 1
complete_dr_means <- data.table::fread( "Processed_Data/processed_coarsescale_dr.tsv")
herits_df <- data.table::fread( "Processed_Data/processed_coarsescale_h2.tsv")
fine_scale_dr_mean <- data.table::fread( "Processed_Data/processed_finescale_dr.tsv")
sampled_h2_df <- data.table::fread( "Processed_Data/subsambled_finescale_h2.tsv")
pl50 <- data.table::fread( "Data/propionate_LD50.csv")

label_df <- dplyr::filter(herits_df, pa_conc == 100)

complete_drplt <- ggplot()+
  geom_smooth(aes(x = pa_conc, y = dose_mean, group = strain), color = "black", fill = "gray90",
              data = complete_dr_means %>% dplyr::filter(!(strain %in% c("N2", "DL238"))),
              span=.8, size = 0.5) +
  geom_smooth(aes(x = pa_conc, y = dose_mean, color = strain), 
              data = complete_dr_means %>% dplyr::filter(strain %in% c("N2", "DL238", "BRC20067")),
              span=.8, size = 0.5) +
  scale_x_continuous(limits = c(0, 150), breaks = unique(complete_dr_means$pa_conc)) +
  scale_y_continuous(breaks = c(0,25,50,75,100)) +
  coord_cartesian(ylim=c(-10, 120)) +
  scale_color_manual(values = c("cadetblue3", "orange")) +
  geom_hline(aes(yintercept = 50), linetype = 2, color = "red", alpha = 0.7) +
  # geom_vline(aes(xintercept = 85), linetype = 2, color = "orange", alpha = 0.7) +
  geom_segment(aes(x = 85, xend=85, y=-20,yend=50), linetype = 2, color = "orange", alpha = 0.7) +
  geom_segment(aes(x = 105, xend=105, y=-20,yend=50), linetype = 2, color = "cadetblue3", alpha = 0.7) +
  theme_classic(14) +
  labs(x = "Propionate concentration (mM)", y = "L1 survival (%)", color = "Strain")

complete_h2 <- ggplot(herits_df)+
  aes(x = pa_conc, y = H2)+
  geom_point(size = 2)+
  theme_classic(14)+
  annotate("text", x = label_df$pa_conc, y = signif(label_df$H2,2)+.05, label = signif(label_df$H2,2), size = 4)+
  scale_x_continuous(limits = c(0, 150), breaks = unique(complete_dr_means$pa_conc)) +
  scale_y_continuous(limits = c(0,1),breaks = c(0,0.25,0.50,0.75,1)) +
  labs(x = "Propionate concentration (mM)", y = "Broad-sense heritability")

fine_scale_dr <- ggplot(fine_scale_dr_mean)+
  geom_smooth(aes(x = n_conc, y = dose_mean, group = strain_name), color = "black", fill = "gray80",
              data = fine_scale_dr_mean %>% dplyr::filter(!(strain %in% c("N2", "DL238", "BRC20067"))),
              span=1, size = 0.5) +
  geom_smooth(aes(x = n_conc, y = dose_mean, color = strain_name),
              data = fine_scale_dr_mean %>% dplyr::filter(strain_name %in% c("N2", "DL238", "BRC20067")),
              span=1, size = 0.5) +
  scale_y_continuous(breaks = c(0,25,50,75,100)) +
  coord_cartesian(ylim=c(-10, 120)) +
  scale_color_manual(values = c("cadetblue3", "orange")) +
  geom_vline(aes(xintercept = 100), linetype = 2, color = "red", alpha = 0.7) +
  labs(x = "Propionate concentration (mM)", y = "L1 survival (%)", color = "Strain")


finescale_h2 <- sampled_h2_df %>%
  dplyr::group_by(Concentration) %>%
  dplyr::summarise(mh2 = mean(H2),
                   sh2 = sd(H2)) %>%
  ggplot() +
  geom_point(aes(x = Concentration, y = mh2),size = 0.5)+
  geom_errorbar(aes(ymin=mh2-sh2, ymax=mh2+sh2,x = Concentration), width=.4,
                position=position_dodge(.9)) +
  scale_y_continuous(limits = c(0,1),breaks = c(0,0.25,0.50,0.75,1)) +
  labs(x = "Propionate concentration (mM)", y = "Broad-sense heritability")
  
 
finescale_h2 <- sampled_h2_df %>%
  dplyr::group_by(Concentration) %>%
  dplyr::summarise(mh2 = median(H2),
                   q2h2 = stats::quantile(H2, probs=0.25),
                   q7h2 =stats::quantile(H2, probs=0.75),
                   iqrh2 =IQR(H2)) %>%
  ggplot() +
  geom_segment(aes(x = Concentration-1, xend = Concentration+1, y=mh2, yend=mh2)) +
  geom_segment(aes(x = Concentration-1, xend = Concentration+1, y=q2h2, yend=q2h2)) +
  geom_segment(aes(x = Concentration-1, xend = Concentration+1, y=q7h2, yend=q7h2)) +
  geom_segment(aes(x = Concentration-1, xend = Concentration-1, y=q2h2, yend=q7h2)) +
  geom_segment(aes(x = Concentration+1, xend = Concentration+1, y=q2h2, yend=q7h2)) +
  geom_segment(aes(x = Concentration, xend = Concentration, y=q2h2, yend=q2h2-1.5*iqrh2)) +
  geom_segment(aes(x = Concentration, xend = Concentration, y=q7h2, yend=q7h2+1.5*iqrh2)) +
  # geom_boxplot(aes(x = Concentration, y = H2),size = 0.5)+
  # geom_errorbar(aes(ymin=mh2-sh2, ymax=mh2+sh2,x = Concentration), width=.4,
  #               position=position_dodge(.9)) +
  scale_y_continuous(limits = c(0,1),breaks = c(0,0.25,0.50,0.75,1)) +
  labs(x = "Propionate concentration (mM)", y = "Broad-sense heritability")
 
# ld50 

l50.lm <- lm(ld50 ~ strain, data = pl50)
l50.av <- aov(l50.lm)
summary(l50.av)
tukey.test <- TukeyHSD(l50.av)

plot(tukey.test)

tukey.test$strain %>%
  data.frame() %>%
  dplyr::mutate(s_comp = row.names(.)) %>%
  dplyr::filter(p.adj < 0.05)

l_compare <- list(c("N2", "DL238"), c("N2", "EG4725"))

# finescale_h2 <- ggplot(sampled_h2_df)+
#   aes(x = n_conc, y = H2)+
#   geom_boxplot(outlier.colour = NA)+
#   ggbeeswarm::geom_beeswarm(size = 0.5)+
#   scale_y_continuous(limits = c(0,1),breaks = c(0,0.25,0.50,0.75,1)) +
#   labs(x = "Propionate concentration (mM)", y = "Broad-sense heritability")

ld50_plot <- pl50 %>%
  dplyr::mutate(strain = factor(strain, levels = c("N2", "CB4856", "CX11314", "DL238", "ED3017", "EG4725", "JT11398", "JU258", "JU775", "LKC34", "MY16", "MY23"))) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(m50 = mean(ld50), 
                s50 = sd(ld50)) %>%
  dplyr::ungroup()%>%
  ggplot()+
  aes(x = strain, y = ld50) +
  geom_bar(aes(y = m50, fill = strain), stat = "identity", color = "black", size = 0.5) +
  geom_point(size = 0.5)+
  scale_fill_manual(values = c("N2" = "orange", "CB4856" = "gray50", "CX11314"= "gray50", 
                               "DL238"= "cadetblue3", "ED3017"= "gray50", "EG4725"= "gray50",
                               "JT11398"= "gray50", "JU258"= "gray50", "JU775"= "gray50", 
                               "LKC34"= "gray50", "MY16"= "gray50", "MY23"= "gray50"))+
  scale_y_continuous(limits = c(0,125),breaks = c(0,20,40, 60, 80, 100, 120), expand = c(0, 0)) +
  theme_classic(12) +
  # stat_compare_means(method = "anova", label.y = 100, label.x = "JU775", size =3,label.sep = "\n", hjust =0)+        # Add global annova p-value
  stat_compare_means(comparisons = l_compare, label.y = c(110, 115),
                     label = "p.signif", method = "t.test",
                     ref.group = "N2", hide.ns = TRUE, size = 5, color = "gray50") +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = expression(paste("Propionate LD"[50], " (mM)")))

left_plots <- cowplot::plot_grid(ggdraw() + draw_image("Plots/SVG_PLOTS/Propionate_Pathway.jpg"),
                                 ld50_plot,
                                 complete_h2+ theme_classic(12),
                                 
                                 ncol = 1, label_size = 14,
                                 labels = c("A","C","E"),scale = c(0.9,1,1))

right_plots <- cowplot::plot_grid(lemon::reposition_legend(complete_drplt+ theme_classic(12) , 'top right'),
                                  fine_scale_dr + 
                                    theme_classic(12)+
                                    xlim(79,121)+
                                    theme(legend.position = "none"),
                                  finescale_h2+ 
                                    theme_classic(12), 
                                  ncol = 1, label_size = 14,
                                  labels = c("B", "D", "F"), align = "hv")


cowplot::plot_grid(left_plots,
                   right_plots,
                   ncol = 2, label_size = 14,
                   labels = NA, align = "hv")

ggsave(filename = "Plots/SVG_PLOTS/Figure1_new.pdf", height = 10, width = 6.5, units = "in")
ggsave(filename = "Plots/SVG_PLOTS/Figure1_new.svg", height = 10, width = 6.5, units = "in")
ggsave(filename = "Plots/SVG_PLOTS/Figure1_new.png", height = 10, width = 6.5, units = "in")


# figure s1

samp.all.df <- data.table::fread("Processed_Data/Power_samples.tsv")

ggplot(samp.all.df)+
  aes(x = d, y = m.sd, 
      color = factor(n),
      fill = factor(n))+
  geom_line(size= 0.25) +
  theme_classic(12) +
  geom_hline(aes(yintercept = 0.8), color = "red", linetype = 2, size = 0.25)+
  geom_ribbon(aes(ymin=m.sd-sd.sd,ymax=m.sd+sd.sd), alpha= 0.25, linetype=4, size = 0 )+
  scale_color_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                     name = "Sample size")+
  scale_fill_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                    name = "Sample size")+
  xlim(0,.75) +
  labs(x = "Difference between means", y = "Power")

ggsave(filename = "Plots/SVG_PLOTS/FigureS1.pdf", height = 4, width = 6)
ggsave(filename = "Plots/SVG_PLOTS/FigureS1.png", height = 4, width = 6, dpi = 300)
ggsave(filename = "Plots/SVG_PLOTS/FigureS1.svg", height = 4, width = 6)


# figure s2

p <- data.table::fread("Data/20161128_new_gwas_replicates/20170103_set1-3.csv")%>%
  dplyr::select(-number, -average)%>%
  tidyr::gather(t_replicate, ratio, -strain, -replicate,-set) %>%
  dplyr::filter(strain != "")

replicated_strains <- p %>%
  dplyr::group_by(strain) %>%
  dplyr::summarise(ct = n()) %>%
  dplyr::arrange(desc(ct)) %>%
  head() 

other_strains <- p %>%
  dplyr::filter(!strain %in% replicated_strains$strain) %>%
  dplyr::distinct(strain)

plot_df <- p %>%
  dplyr::mutate(strain_f = factor(strain, levels = c(replicated_strains$strain, other_strains$strain))) %>%
  dplyr::group_by(strain, set) %>%
  dplyr::mutate(md_s = median(ratio, na.rm = T)) %>%
  dplyr::mutate(rep_str = ifelse(strain %in%replicated_strains$strain, "yes","no"))


batch1 <- plot_df %>%
  dplyr::filter(set == 1) %>%
  ggplot()+
  aes(x = strain_f, y = ratio)+
  geom_point(size = 0.25) +
  geom_point(aes(y = md_s, color = rep_str), shape = 45, size = 10, 
             data = plot_df %>%
               dplyr::filter(set == 1) %>%
               distinct(strain, md_s, .keep_all=T)) +
  scale_color_manual(values=c("yes" = "red", "no" = "blue"))+
  theme_classic(12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.title.x = element_blank(), 
        legend.position = "none")+
  scale_y_continuous(limits = c(0,100), expand = c(0.01,0)) +
  labs(y = "Batch 1\nSurvival rate (%)")

batch2 <- plot_df %>%
  dplyr::filter(set == 2) %>%
  ggplot()+
  aes(x = strain_f, y = ratio)+
  geom_point(size = 0.25) +
  geom_point(aes(y = md_s, color = rep_str), shape = 45, size = 10, 
             data = plot_df %>%
               dplyr::filter(set == 2) %>%
               distinct(strain, md_s, .keep_all=T)) +
  scale_color_manual(values=c("yes" = "red", "no" = "blue"))+
  theme_classic(12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.title.x = element_blank(), 
        legend.position = "none")+
  scale_y_continuous(limits = c(0,100), expand = c(0.01,0)) +
  labs(y = "Batch 2\nSurvival rate (%)")

batch3 <- plot_df %>%
  dplyr::filter(set == 3) %>%
  ggplot()+
  aes(x = strain_f, y = ratio)+
  geom_point(size = 0.25) +
  geom_point(aes(y = md_s, color = rep_str), shape = 45, size = 10, 
             data = plot_df %>%
               dplyr::filter(set == 3) %>%
               distinct(strain, md_s, .keep_all=T)) +
  scale_color_manual(values=c("yes" = "red", "no" = "blue"))+
  theme_classic(12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.title.x = element_blank(), 
        legend.position = "none")+
  scale_y_continuous(limits = c(0,100), expand = c(0.01,0)) +
  labs(y = "Batch 3\nSurvival rate (%)")

bottom_batches <- cowplot::plot_grid(batch1,batch2,batch3,
                                     align = "h",
                                     labels = c("C",NA,NA), 
                                     label_size = 14,ncol = 1)

top_left <- plot_df %>%
  dplyr::filter(strain %in% replicated_strains$strain) %>%
  ggplot()+
  aes(x = strain_f, y = ratio)+
  geom_boxplot(outlier.colour = NA)+
  ggbeeswarm::geom_quasirandom(size = 0.25)+
  theme_classic(12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank(), 
        legend.position = "none")+
  scale_y_continuous(limits = c(0,100), expand = c(0.01,0)) +
  labs(y = "Survival rate (%)")


top_plot <- cowplot::plot_grid(top_left,
                               ggdraw() + draw_image("Plots/SVG_PLOTS/L1_Survival_Assay.jpg"),
                                     labels = c("A","B"), rel_widths =c(1.5,2.5), scale = c(1,0.9),
                                     label_size = 14,ncol = 2)

cowplot::plot_grid(top_plot,bottom_batches,
                   labels = c(NA,NA), rel_heights =c(1.25,3), align = "v",
                   label_size = 14,ncol = 1)

ggsave(filename = "Plots/SVG_PLOTS/FigureS2.pdf", height = 10, width = 7)
ggsave(filename = "Plots/SVG_PLOTS/FigureS2.png", height = 10, width = 7, dpi = 300)
ggsave(filename = "Plots/SVG_PLOTS/FigureS2.svg", height = 10, width = 7)
