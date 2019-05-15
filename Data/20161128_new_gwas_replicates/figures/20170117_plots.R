
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

raw <- p%>%
  dplyr::filter(strain %in% c("MY23","JU830","JU775","EG4725","DL238","CX11307"))%>%
  dplyr::group_by(strain,set)%>%
  dplyr::mutate(means = mean(ratio, na.rm = T))



raw <- p%>%
  dplyr::filter(strain %in% c("MY23","JU830","JU775","EG4725","DL238","CX11307"))%>%
  dplyr::group_by(strain,replicate, set)%>%
  dplyr::mutate(means = mean(ratio, na.rm = T))

ggplot(mixed_p) + 
  aes(y = log(ratio), x =strain )+
  geom_boxplot(fill = "gray", outlier.colour = NA)+  
  # geom_jitter(aes(color=factor(replicate)))+
  # geom_point(aes(x=strain, y = means,fill=factor(set)), shape  =23, size = 3, alpha = .5)+
  theme_bw()+ggplot2::theme(axis.text.x = ggplot2::element_text(size = 14, face = "bold", color = "black"),
                                             axis.text.y = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                                             axis.title.x = ggplot2::element_text(size = 0, face = "bold", color = "black", vjust = -0.3), 
                                             axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                                             strip.text.x = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                                             strip.text.y = ggplot2::element_text(size = 16,  face = "bold", color = "black"), 
                                             plot.title = ggplot2::element_text(size = 0,  face = "bold", vjust = 1),
                                             panel.background = ggplot2::element_rect(color = "black",  size = 1.2), 
                                             strip.background = ggplot2::element_rect(color = "black", size = 1.2))+
  labs(x="Strain", y="Raw Phenotype")

ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/figures/control-strain-boxplots.pdf",
       height = 4,
       width = 8)

p_pr <-  p%>%
  dplyr::filter(strain != "")%>%
  dplyr::group_by(strain,replicate,set)%>%
  dplyr::mutate(m_rat = mean(ratio, na.rm = T))%>%
  dplyr::select( replicate, set, strain, m_rat)%>%
  dplyr::ungroup()%>%
  dplyr::distinct(set, strain, m_rat,.keep_all = T)%>%
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


# raw <- pr_resid%>%
#   dplyr::filter(strain %in% c("MY23","JU830","JU775","EG4725","DL238","CX11307"))%>%
#   dplyr::group_by(strain)%>%
#   dplyr::mutate(means = mean(resid, na.rm = T))
# 
# 
# ggplot(raw) + 
#   aes(y = resid, x =strain )+
#   geom_boxplot(fill = "gray", outlier.colour = NA)+  
#   geom_jitter()+
#   # geom_point(aes(x=strain, y = means), shape  =23, size = 3, alpha = .5)+
#   theme_bw()+ggplot2::theme(axis.text.x = ggplot2::element_text(size = 14, face = "bold", color = "black"),
#                             axis.text.y = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
#                             axis.title.x = ggplot2::element_text(size = 0, face = "bold", color = "black", vjust = -0.3), 
#                             axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
#                             strip.text.x = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
#                             strip.text.y = ggplot2::element_text(size = 16,  face = "bold", color = "black"), 
#                             plot.title = ggplot2::element_text(size = 0,  face = "bold", vjust = 1),
#                             panel.background = ggplot2::element_rect(color = "black",  size = 1.2), 
#                             strip.background = ggplot2::element_rect(color = "black", size = 1.2))+
#   labs(x="Strain", y="Raw Phenotype")
# 
# ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/figures/control-strain-boxplots.pdf",
#        height = 4,
#        width = 8)

pr_resid$strain <- factor(pr_resid$strain,levels = pr_resid$strain,labels = pr_resid$strain, ordered = T)

ggplot(pr_resid) + 
  aes(y = resid, x =strain , fill = strain)+
  geom_boxplot()+  theme_bw()+ggplot2::theme(axis.text.x = ggplot2::element_text(size = 0, face = "bold", color = "black"),
                                  axis.text.y = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                                  axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                                  axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                                  strip.text.x = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                                  strip.text.y = ggplot2::element_text(size = 16,  face = "bold", color = "black"), 
                                  plot.title = ggplot2::element_text(size = 0,  face = "bold", vjust = 1),
                                  panel.background = ggplot2::element_rect(color = "black",  size = 1.2), 
                                  strip.background = ggplot2::element_rect(color = "black", size = 1.2),
                                  legend.position = 'none')+
  labs(x="Strain", y="Processed Phenotype")

hist <-ggplot(pr_resid)+
  aes(x = resid)+
  geom_histogram(color = "black", fill = "gray50")+
  # scale_fill_manual(values = c("ALT" = "cyan", "REF" = "gray"))+
  # scale_color_manual(values = c("ALT" = "cyan", "REF" = "gray"))+
  ggplot2::theme_bw() + 
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16, face = "bold", color = "black"),
                 axis.text.y = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                 axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                 axis.title.y = ggplot2::element_text(size = 20, color = "black"), 
                 strip.text.x = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                 strip.text.y = ggplot2::element_text(size = 16,  face = "bold", color = "black"), 
                 plot.title = ggplot2::element_text(size = 0,  face = "bold", vjust = 1),
                 panel.background = ggplot2::element_rect(color = "black",  size = 1.2), 
                 strip.background = ggplot2::element_rect(color = "black", size = 1.2),
                 legend.position = 'none')+
  labs(x = "Propionate Survival", y = "Count")

ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/figures/phenotype_histogram.pdf",
       height = 4,
       width = 6)

load( file ="~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/3set_pr_maps.Rda")


mp <- manplot(pr_maps)

mp <- mp[[1]] +  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16, face = "bold", color = "black"),
                                axis.text.y = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                                axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                                axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                                strip.text.x = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                                strip.text.y = ggplot2::element_text(size = 16,  face = "bold", color = "black"), 
                                plot.title = ggplot2::element_text(size = 0,  face = "bold", vjust = 1),
                                panel.background = ggplot2::element_rect(color = "black",  size = 1.2), 
                                strip.background = ggplot2::element_rect(color = "black", size = 1.2),
                                legend.position = 'none')

mp

ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/figures/manplot.pdf",
       height = 4,
       width = 10)


pxg <- pxg_plot(pr_maps,color_strains = c("CB4856","N2"))

pxg[[1]] +  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16, face = "bold", color = "black"),
                          axis.text.y = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                          axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                          axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                          strip.text.x = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                          strip.text.y = ggplot2::element_text(size = 16,  face = "bold", color = "black"), 
                          plot.title = ggplot2::element_text(size = 0,  face = "bold", vjust = 1),
                          panel.background = ggplot2::element_rect(color = "black",  size = 1.2), 
                          strip.background = ggplot2::element_rect(color = "black", size = 1.2),
                          legend.position = 'none')

ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/figures/peak_pxg_split.pdf",
       height = 4,
       width = 6)

load("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/genes.Rda")

gene_table <- dplyr::select(genes, CHROM,POS,REF,ALT,nt_change, aa_change, gene_name, gene_id,molecular_name, strain, GT, trait, pheno_value, abs_spearman_cor)

write.table(gene_table, file = "~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/genes_verbose.csv",
            quote=FALSE, sep=",", row.names = F)

gene_table <- dplyr::select(genes, CHROM,POS,REF,ALT,nt_change, aa_change, gene_name, gene_id,molecular_name, strain, GT, trait, pheno_value, abs_spearman_cor)%>%
  dplyr::filter(strain == "CB4856")

write.table(gene_table, file = "~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/genes_concise.csv",
            quote=FALSE, sep=",", row.names = F)

genes %>% 
  # dplyr::filter(!grepl("stream|syn",effect), strain == "CB4856")%>%
  dplyr::filter(strain == "CB4856")%>%
  ggplot(.)+
  aes(x = POS/1e6, y = abs_spearman_cor, fill = GT)+
  geom_point(shape=21, size =2, alpha = .7)+
  scale_fill_manual(values = c("ALT" = "cyan", "REF" = "gray50"))+
  scale_color_manual(values = c("ALT" = "cyan", "REF" = "gray50"))+
  xlim(3.8, 4.06)+
  ggplot2::theme_bw() + 
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 14, face = "bold", color = "black"),
                 axis.text.y = ggplot2::element_text(size = 14, face = "bold", color = "black"), 
                 axis.title.x = ggplot2::element_text(size = 16, face = "bold", color = "black", vjust = -0.3), 
                 axis.title.y = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                 strip.text.x = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                 strip.text.y = ggplot2::element_text(size = 16,  face = "bold", color = "black"), 
                 plot.title = ggplot2::element_text(size = 0,  face = "bold", vjust = 1),
                 panel.background = ggplot2::element_rect(color = "black",  size = 1.2), 
                 strip.background = ggplot2::element_rect(color = "black", size = 1.2),
                 legend.position = 'none')+
  labs(x = "Genomic Position (Mb)", y = expression( bold(paste("Spearman's ",bolditalic("rho")))))

ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/figures/fine_map.pdf",
       height = 4,
       width = 6)