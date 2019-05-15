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

p_pr <-  p%>%
  dplyr::filter(strain != "")%>%
  dplyr::group_by(strain,replicate,set)%>%
  dplyr::mutate(m_rat = mean(ratio, na.rm = T))%>%
  dplyr::select( replicate, set, strain, m_rat)%>%
  dplyr::ungroup()%>%
  dplyr::distinct(set, strain,replicate, m_rat,.keep_all = T)%>%
  dplyr::group_by(strain)%>%
  dplyr::mutate(mph = median(m_rat, na.rm = T),
                sph = sd(m_rat, na.rm = T))%>%
  dplyr::filter(m_rat <= 1.5*sph+mph,
                m_rat >= mph-1.5*sph)

pr_resid <- p_pr%>%
  ungroup()%>%
  dplyr::do(augment(lm(m_rat ~ set + replicate , data =.)))%>%
  dplyr::rename(resid = .resid)%>%
  data.frame()%>%
  dplyr::left_join(p_pr, ., by = c("m_rat","set","replicate"))%>%
  dplyr::distinct(m_rat,set,replicate,strain,.keep_all = T) %>%
  dplyr::ungroup()%>%
  dplyr::select(strain, resid)%>%
  dplyr::group_by(strain)%>%
  dplyr::group_by(strain)%>%
  dplyr::mutate(mph = median(resid, na.rm = T),
                sph = sd(resid, na.rm = T))%>%
  dplyr::filter(resid <= 1.5*sph+mph,
                resid >= mph-1.5*sph)%>%
  dplyr::mutate(m.p = mean(resid))%>%
  dplyr::arrange(desc(m.p))

pr_resid$strain <- factor(pr_resid$strain,levels = pr_resid$strain,labels = pr_resid$strain, ordered = T)

ggplot(pr_resid) + 
  aes(y = resid, x =strain , fill = strain)+
  geom_boxplot()

p_reg <- pr_resid%>%
  dplyr::summarise(phenotype = mean(resid, na.rm = T))

pr_pheno <- process_pheno(p_reg)
maps <- gwas_mappings(pr_pheno)

pr_maps <- process_mappings(maps,pr_pheno)

save(pr_maps, file ="~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/3set_pr_maps.Rda")

inter_sum <- interval_summary("V:3158728-4323750")

save(inter_sum, file ="~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/interval_summary.Rda")


manplot(pr_maps)
pxg_plot(pr_maps,color_strains = c("CB4856", "JU258","ED3049","N2"))

genes <- process_correlations(variant_correlation(pr_maps,condition_trait = F))

save(genes, file ="~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/genes.Rda")


pr_genes <- genes%>%
  dplyr::filter(strain == "CB4856")%>%
  dplyr::distinct(CHROM,POS,aa_change, .keep_all = T)

save(pr_genes, file ="~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/pr_genes.Rda")

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
