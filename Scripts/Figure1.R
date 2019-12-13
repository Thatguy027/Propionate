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
  geom_point(size = 0.5)+
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

ld50_plot <-pl50 %>%
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
  coord_flip() +
  theme(legend.position = "none", 
        axis.title.y = element_blank()) +
  labs(y = expression(paste("Propionate LD"[50], " (mM)")))



# finescale_h2 <- ggplot(sampled_h2_df)+
#   aes(x = n_conc, y = H2)+
#   geom_boxplot(outlier.colour = NA)+
#   ggbeeswarm::geom_beeswarm(size = 0.5)+
#   scale_y_continuous(limits = c(0,1),breaks = c(0,0.25,0.50,0.75,1)) +
#   labs(x = "Propionate concentration (mM)", y = "Broad-sense heritability")

top_row <- cowplot::plot_grid(NULL,
                              ld50_plot,
                              ncol = 2, label_size = 14,
                              labels = c("A","B"), rel_heights = c(1,2))

coarse_cow <- cowplot::plot_grid(lemon::reposition_legend(complete_drplt+ theme_classic(12)+theme(axis.text.x = element_blank(), axis.title.x = element_blank()) , 'top right'),
                                 complete_h2+ theme_classic(12),
                                 ncol = 1, label_size = 14,
                                 labels = c("C", "D"), align = "hv")

fine_cow <- cowplot::plot_grid(fine_scale_dr + 
                                 theme_classic(12)+
                                 xlim(79,121)+
                                 theme(axis.text.x = element_blank(), 
                                       axis.title.x = element_blank(), legend.position = "none"),
                               finescale_h2 + 
                                 theme_classic(12) ,
                               ncol = 1, label_size = 14,
                               labels = c("E", "F"), align = "hv")

bottom <- cowplot::plot_grid(coarse_cow,
                             fine_cow,
                             ncol = 2, label_size = 14,
                             labels = NA)


cowplot::plot_grid(top_row,
                   bottom,
                   ncol = 1, label_size = 14,
                   labels = NA, rel_heights = c(1,2))

ggsave(filename = "Plots/SVG_PLOTS/Figure1_up.pdf", height = 10, width = 6.5, units = "in")
ggsave(filename = "Plots/SVG_PLOTS/Figure1_up.svg", height = 10, width = 6.5, units = "in")
ggsave(filename = "Plots/SVG_PLOTS/Figure1_up.png", height = 10, width = 6.5, units = "in")


# figure s1

samp.all.df <- data.table::fread("Processed_Data/Power_samples.tsv")

ggplot(samp.all.df)+
  aes(x = d, y = m.sd, 
      color = factor(n),
      fill = factor(n))+
  geom_line() +
  theme_classic(18) +
  geom_hline(aes(yintercept = 0.8), color = "red")+
  geom_ribbon(aes(ymin=m.sd-sd.sd,ymax=m.sd+sd.sd), alpha= 0.25, linetype=4, size = 0 )+
  scale_color_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                     name = "Sample size")+
  scale_fill_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                    name = "Sample size")+
  xlim(0,.75) +
  labs(x = "Difference between means", y = "Power")

ggsave(filename = "Plots/Power_Calculation.pdf", height = 6, width = 10)
ggsave(filename = "Plots/Power_Calculation.png", height = 6, width = 10, dpi = 300)
ggsave(filename = "Plots/Power_Calculation.svg", height = 6, width = 10)
