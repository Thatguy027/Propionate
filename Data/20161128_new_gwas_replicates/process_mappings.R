library(data.table)
library(dplyr)
library(tidyr)
library(cegwas)
library(ggplot2)
library(broom)

p <- data.table::fread("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170103_set1-3.csv")%>%
  dplyr::select(-number, -average)%>%
  tidyr::gather(t_replicate, ratio, -strain, -replicate,-set)


p_pr <- p%>%
  dplyr::group_by(strain,t_replicate,set)%>%
  dplyr::mutate(m_rat = mean(ratio))

ggplot(p_pr)+
  aes(x= strain, y = ratio)+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 0, face = "bold", color = "black"),
                 axis.text.y = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                 axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                 axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.x = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                 strip.text.y = ggplot2::element_text(size = 16,  face = "bold", color = "black"), 
                 plot.title = ggplot2::element_text(size = 0,  face = "bold", vjust = 1),
                 panel.background = ggplot2::element_rect(color = "black",  size = 1.2), 
                 strip.background = ggplot2::element_rect(color = "black", size = 1.2),
                 legend.position = 'none')


ctr_strains <- c("CX11307" ,
                 "JU830" ,
                 "EG4725",
                 "JU775",
                 "MY23",
                 "DL238")

p <- data.table::fread("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170103_set1-3.csv")%>%
  dplyr::select(-number, -average)%>%
  tidyr::gather(t_replicate, ratio, -strain, -replicate, -set)

p_ctr <- dplyr::filter(p,strain %in% ctr_strains)%>%
  dplyr::group_by(strain,t_replicate,set)%>%
  dplyr::mutate(m_rat = mean(ratio))

ggplot(p_ctr)+
  aes(x= strain, y = m_rat)+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(aes(color = factor(set)))+
  theme_bw()+
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 0, face = "bold", color = "black"),
                 axis.text.y = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                 axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                 axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.x = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                 strip.text.y = ggplot2::element_text(size = 16,  face = "bold", color = "black"), 
                 plot.title = ggplot2::element_text(size = 0,  face = "bold", vjust = 1),
                 panel.background = ggplot2::element_rect(color = "black",  size = 1.2), 
                 strip.background = ggplot2::element_rect(color = "black", size = 1.2),
                 legend.position = 'none')




p_pr <-  p%>%
  dplyr::filter(strain != "")%>%
  dplyr::group_by(strain,replicate,set)%>%
  dplyr::mutate(m_rat = mean(ratio, na.rm = T))%>%
  mutate(med1 = m_rat+ 1.5*IQR(ratio, na.rm = T))%>%
  dplyr::group_by(strain)%>%
  dplyr::mutate(mph = median(ratio, na.rm = T),
                sph = sd(ratio, na.rm = T))%>%
  dplyr::filter(ratio <= 1.5*sph+mph,
                ratio >= mph-1.5*sph)

ggplot(p_pr)+
  aes(x= strain, y = ratio)+
  geom_boxplot()+
  geom_jitter(aes(color = factor(set)))+
  theme_bw()+
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 0, face = "bold", color = "black"),
                 axis.text.y = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                 axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                 axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.x = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                 strip.text.y = ggplot2::element_text(size = 16,  face = "bold", color = "black"), 
                 plot.title = ggplot2::element_text(size = 0,  face = "bold", vjust = 1),
                 panel.background = ggplot2::element_rect(color = "black",  size = 1.2), 
                 strip.background = ggplot2::element_rect(color = "black", size = 1.2),
                 legend.position = 'none')


p_reg <- dplyr::select(p_pr, replicate, set, strain, m_rat)%>%
  dplyr::distinct(set, strain, m_rat,.keep_all = T)%>%
  dplyr::do(augment(lm(m_rat ~ set + replicate, data =.)))%>%
  dplyr::rename(resid = .resid)%>%
  data.frame()

ggplot(p_reg)+
  aes(x= strain, y = resid)+
  geom_boxplot()+
  geom_jitter(aes(color = factor(set)))+
  theme_bw()+
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 0, face = "bold", color = "black"),
                 axis.text.y = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                 axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                 axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.x = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                 strip.text.y = ggplot2::element_text(size = 16,  face = "bold", color = "black"), 
                 plot.title = ggplot2::element_text(size = 0,  face = "bold", vjust = 1),
                 panel.background = ggplot2::element_rect(color = "black",  size = 1.2), 
                 strip.background = ggplot2::element_rect(color = "black", size = 1.2),
                 legend.position = 'none')

p_reg <- dplyr::select(p_pr, replicate, set, strain, m_rat)%>%
  dplyr::distinct(set, strain, m_rat,.keep_all = T)%>%
  dplyr::do(augment(lm(m_rat ~ replicate , data =.)))%>%
  dplyr::rename(resid = .resid)%>%
  data.frame()

ggplot(p_reg)+
  aes(x= strain, y = resid)+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 0, face = "bold", color = "black"),
                 axis.text.y = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                 axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                 axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.x = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                 strip.text.y = ggplot2::element_text(size = 16,  face = "bold", color = "black"), 
                 plot.title = ggplot2::element_text(size = 0,  face = "bold", vjust = 1),
                 panel.background = ggplot2::element_rect(color = "black",  size = 1.2), 
                 strip.background = ggplot2::element_rect(color = "black", size = 1.2),
                 legend.position = 'none')


p_reg <- dplyr::select(p_pr, replicate, set, strain, m_rat)%>%
  dplyr::ungroup()%>%
  dplyr::distinct(set, strain, m_rat,.keep_all = T)%>%
  dplyr::do(augment(lm(m_rat ~ set + replicate , data =.)))%>%
  dplyr::rename(resid = .resid)%>%
  data.frame()%>%
  dplyr::left_join(p_pr, ., by = c("m_rat","set","replicate"))%>%
  dplyr::distinct(m_rat,set,replicate,strain,.keep_all = T)

ggplot(p_reg)+
  aes(x= strain, y = resid)+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 0, face = "bold", color = "black"),
                 axis.text.y = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                 axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                 axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.x = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                 strip.text.y = ggplot2::element_text(size = 16,  face = "bold", color = "black"), 
                 plot.title = ggplot2::element_text(size = 0,  face = "bold", vjust = 1),
                 panel.background = ggplot2::element_rect(color = "black",  size = 1.2), 
                 strip.background = ggplot2::element_rect(color = "black", size = 1.2),
                 legend.position = 'none')

p_means <- p_reg%>%
  dplyr::ungroup()%>%
  dplyr::select(strain, resid)%>%
  dplyr::group_by(strain)%>%
  dplyr::summarise(ratio = log(mean(resid)))

pr_pheno <- process_pheno(p_means)
maps <- gwas_mappings(pr_pheno)

pr_maps <- process_mappings(maps,pr_pheno)

save(pr_maps, file = "~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/processed_mappings.Rda")

manplot(pr_maps)

pxg_plot(pr_maps,color_strains = NA)

genes <- process_correlations(variant_correlation(pr_maps,condition_trait = F))

pr_genes <- genes%>%
  # dplyr::filter(strain == "N2")%>%
  dplyr::distinct(CHROM,POS,aa_change, .keep_all = T)

save(genes, file = "~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/processed_genes.Rda")


fine_pl<-ggplot(pr_genes)+
  aes(x = POS/1e6, y = abs_spearman_cor, fill = GT)+
  geom_point(shape=21)+
  scale_fill_manual(values = c("ALT" = "cyan", "REF" = "cyan"))+
  scale_color_manual(values = c("ALT" = "cyan", "REF" = "cyan"))+
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


fine_pl