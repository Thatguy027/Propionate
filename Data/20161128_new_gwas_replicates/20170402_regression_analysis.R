library(data.table)
library(dplyr)
library(tidyr)
library(cegwas)
library(ggplot2)
library(broom)
library(cowplot)

p <- data.table::fread("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170103_set1-3.csv")%>%
  dplyr::select(-number, -average)%>%
  tidyr::gather(t_replicate, ratio, -strain, -replicate,-set)

# pull out strains that were repeated multiple days, take median of replicates (mean & median produces same QTL)
rep_strains <- p %>%
  dplyr::filter(strain %in% c("CX11307", "DL238", "EG4725", "JU775", "JU830", "MY23"))%>%
  dplyr::group_by(strain,set,replicate)%>%
  dplyr::mutate(m_rat = median(ratio, na.rm = T))%>%
  dplyr::select( replicate, set, strain, m_rat)%>%
  dplyr::ungroup()%>%
  dplyr::group_by(strain)

# fit linear model to see how set number affects survival
set_fit <- lm(m_rat ~ set, data = rep_strains)
coeffs <- coef(set_fit)
# visualize model
plot(rep_strains$set, rep_strains$m_rat)
abline(8.280132, 3.803452, col = "red")

# correct phenotypes of all GWAS strains based on the parameters identified above
pr_p <- p %>%
  dplyr::group_by(strain,set,replicate)%>%
  dplyr::mutate(m_rat = median(ratio, na.rm = T))%>%
  dplyr::select( replicate, set, strain, m_rat)%>%
  dplyr::ungroup()%>%
  dplyr::group_by(strain)

# adjust
adjusted_survival <- pr_p %>%
  dplyr::mutate(set_adj = coeffs[1] + coeffs[2]*m_rat  )

# visualize to make sure that the uncorrected vs. corrected values correspond to the fit line identified in the model above. 
plot(adjusted_survival$m_rat, adjusted_survival$set_adj)
abline(8.280132, 3.803452, col = "red")


adjusted_survival_pr <-  adjusted_survival%>%
  dplyr::select(strain, replicate, set_adj)%>%
  dplyr::filter(strain != "")%>%
  dplyr::group_by(strain,replicate)%>%
  dplyr::mutate(m_rat = mean(set_adj, na.rm = T))%>%
  dplyr::select( replicate, strain, m_rat)%>%
  dplyr::ungroup()%>%
  dplyr::distinct(strain, m_rat,replicate,.keep_all = T)%>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(mph = median(m_rat, na.rm = T),
                sph = sd(m_rat, na.rm = T))%>%
  dplyr::filter(m_rat <= 1.5*sph+mph,
                m_rat >= mph-1.5*sph)


mean_adjusted <- adjusted_survival_pr%>%
  dplyr::summarise(phenotype = mean(m_rat, na.rm =T))

pr_pheno <- process_pheno(mean_adjusted)
maps <- gwas_mappings(pr_pheno)

pr_maps <- process_mappings(maps,pr_pheno, BF = 4)

save(pr_maps, file ="~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170402_regression_analysis_output/3set_pr_maps.Rda")



manplot(pr_maps)
ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170402_regression_analysis_output/manplot.pdf",
       width = 10,
       height = 5)
pxg_plot(pr_maps,color_strains = c("CX11307", "DL238", "EG4725", "JU775"))

ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170402_regression_analysis_output/pxg.pdf",
       width = 6,
       height = 5)

genes <- process_correlations(variant_correlation(pr_maps,condition_trait = F,quantile_cutoff_high = .5, variant_severity = "ALL"))

save(genes, file ="~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170402_regression_analysis_output/genes.Rda")


pr_genes <- genes%>%
  dplyr::filter(strain == "CB4856")%>%
  dplyr::distinct(CHROM,POS,aa_change, .keep_all = T)

save(pr_genes, file ="~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170402_regression_analysis_output/pr_genes.Rda")

ggplot(pr_genes)+
  aes(x = POS/1e6, y = abs_spearman_cor, fill = GT)+
  geom_point(shape=21)+
  scale_fill_manual(values = c("ALT" = "cyan", "REF" = "gray"))+
  scale_color_manual(values = c("ALT" = "cyan", "REF" = "gray"))+
  ggplot2::theme_bw() + 
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16, face = "bold", color = "black"),
                 axis.text.y = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                 axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                 axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.x = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                 strip.text.y = ggplot2::element_text(size = 16,  face = "bold", color = "black"), 
                 plot.title = ggplot2::element_text(size = 0,  face = "bold", vjust = 1),
                 panel.background = ggplot2::element_rect(color = "black",  size = 1.2), 
                 strip.background = ggplot2::element_rect(color = "black", size = 1.2),
                 legend.position = 'none')+
  labs(x = "Genomic Position (Mb)", y = expression(bold(paste("Spearman's ", bolditalic("rho")))))




hist <-ggplot(p_reg)+
  aes(x = phenotype)+
  geom_histogram(color = "black", fill = "gray50")+
  # scale_fill_manual(values = c("ALT" = "cyan", "REF" = "gray"))+
  # scale_color_manual(values = c("ALT" = "cyan", "REF" = "gray"))+
  ggplot2::theme_bw() + 
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 14, face = "bold", color = "black"),
                 axis.text.y = ggplot2::element_text(size = 14, face = "bold", color = "black"), 
                 axis.title.x = ggplot2::element_text(size = 16, face = "bold", color = "black", vjust = -0.3), 
                 axis.title.y = ggplot2::element_text(size = 18, color = "black"), 
                 strip.text.x = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                 strip.text.y = ggplot2::element_text(size = 16,  face = "bold", color = "black"), 
                 plot.title = ggplot2::element_text(size = 0,  face = "bold", vjust = 1),
                 panel.background = ggplot2::element_rect(color = "black",  size = 1.2), 
                 strip.background = ggplot2::element_rect(color = "black", size = 1.2),
                 legend.position = 'none')+
  labs(x = "Regressed Propionate Survival", y = "Count")

load( file ="~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/3set_pr_maps.Rda")

mp <- manplot(pr_maps)

mp <- mp[[1]] +  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 14, face = "bold", color = "black"),
                                axis.text.y = ggplot2::element_text(size = 14, face = "bold", color = "black"), 
                                axis.title.x = ggplot2::element_text(size = 16, face = "bold", color = "black", vjust = -0.3), 
                                axis.title.y = ggplot2::element_text(size = 18, face = "bold", color = "black"), 
                                strip.text.x = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                                strip.text.y = ggplot2::element_text(size = 16,  face = "bold", color = "black"), 
                                plot.title = ggplot2::element_text(size = 0,  face = "bold", vjust = 1),
                                panel.background = ggplot2::element_rect(color = "black",  size = 1.2), 
                                strip.background = ggplot2::element_rect(color = "black", size = 1.2),
                                legend.position = 'none')




ggdraw() +
  draw_plot(hist, 0, 0, .2667, .925) +
  draw_plot(mp, .2667, 0, 1-.2667, 1) +
  draw_plot_label(c("A", "B"), c(0, .2667), c(1, 1), size = 16)

cowplot::ggsave(filename = "~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/grant_figure.png",
                plot = last_plot(),
                width = 15,
                height = 4)
