library(cegwas)
library(lmer)
library(lmerTest)
library(broom)

pxg_plot1<-function (plot_df, loc = NA, use_base = F, color_strains = c("N2", 
                                                                        "CB4856")) 
{
  plot_traits <- unique(plot_df$trait)
  plots <- lapply(plot_traits, function(x) {
    plot_peak <- plot_df %>% na.omit() %>% dplyr::filter(trait == 
                                                           x) %>% dplyr::distinct(strain, value, peakPOS, .keep_all = T) %>% 
      dplyr::select(strain, value, CHROM, POS = peakPOS, 
                    -allele) %>% dplyr::mutate(chr_pos = paste(CHROM, 
                                                               POS, sep = ":"))
    if (is.na(loc)) {
      loc <- plot_peak %>% dplyr::select(CHROM, POS) %>% 
        dplyr::distinct(.keep_all = T) %>% dplyr::transmute(chr_pos = paste0(CHROM, 
                                                                             ":", POS))
    }
    to_plot <- snpeff(loc[1], severity = "ALL", elements = "ALL") %>% 
      dplyr::select(strain, CHROM, POS, GT, REF, ALT) %>% 
      dplyr::distinct(.keep_all = T) %>% dplyr::mutate(chr_pos = paste0(CHROM, 
                                                                        ":", POS))
    to_plot <- dplyr::left_join(to_plot, plot_peak) %>% dplyr::mutate(chr_pos = paste(CHROM, 
                                                                                      POS, sep = ":"))
    to_plot <- dplyr::filter(to_plot, !is.na(value)) %>% 
      dplyr::distinct(strain, value, POS, .keep_all = T) %>% 
      dplyr::filter(!is.na(GT)) %>% dplyr::group_by(GT) %>% 
      dplyr::mutate(GT2 = ifelse(use_base, ifelse(GT == 
                                                    "REF", REF, ALT), GT)) %>% dplyr::ungroup() %>% 
      dplyr::mutate(GT = GT2)
    if (!unique(is.na(color_strains))) {
      to_plot <- to_plot %>% dplyr::mutate(colors = ifelse(strain %in% 
                                                             color_strains, strain, "aa_ignore")) %>% dplyr::arrange(colors)
    }
    else {
      to_plot <- to_plot %>% dplyr::mutate(colors = "aa_ignore")
    }
    to_plot %>% ggplot2::ggplot(., ggplot2::aes(x = GT, y = value,  fill = GT)) + 
      ggplot2::scale_fill_brewer(palette = "Set1") + 
      ggplot2::geom_boxplot(outlier.colour = NA) + ggplot2::theme_bw() + 
      ggplot2::geom_jitter(ggplot2::aes(color = colors, size = colors)) +
      ggplot2::scale_color_manual(values = c("black", "#2474FF", "#EE8F03", "red","pink","cyan", "green"), 
                                  labels = c("Other",  unique(to_plot$colors)[2], unique(to_plot$colors)[3], unique(to_plot$colors)[4], unique(to_plot$colors)[5], unique(to_plot$colors)[6], unique(to_plot$colors)[7])) + 
      ggplot2::scale_size_manual(values = c(2, 12, 12, 12,12,12, 12), labels = c("Other", unique(to_plot$colors)[2], 
                                                                unique(to_plot$colors)[3], unique(to_plot$colors)[4], 
                                                                unique(to_plot$colors)[5],unique(to_plot$colors)[6],unique(to_plot$colors)[7])) + ggplot2::facet_grid(. ~ 
                                                                                                                    chr_pos, scales = "free") + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 24, 
                                                                                                                                                                                                   face = "bold", color = "black"), axis.text.y = ggplot2::element_text(size = 24, 
                                                                                                                                                                                                                                                                        face = "bold", color = "black"), axis.title.x = ggplot2::element_text(size = 24, 
                                                                                                                                                                                                                                                                                                                                              face = "bold", color = "black", vjust = -0.3), axis.title.y = ggplot2::element_text(size = 24, 
                                                                                                                                                                                                                                                                                                                                                                                                                                  face = "bold", color = "black"), strip.text.x = ggplot2::element_text(size = 24, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        face = "bold", color = "black"), strip.text.y = ggplot2::element_text(size = 16, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              face = "bold", color = "black"), plot.title = ggplot2::element_text(size = 24, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  face = "bold", vjust = 1), panel.background = ggplot2::element_rect(color = "black", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      size = 1.2), strip.background = ggplot2::element_rect(color = "black", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            size = 1.2)) + ggplot2::labs(y = "Phenotype", x = "Genotype", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         title = x)
  })
  plots
}

load("/Users/stefan/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170403_simple_set-and-ratioRegression_maps.Rda")

pxg_plot1(pr_maps,color_strains = c("BRC20067", "DL238", "EG4725", "JU775","JU830"))[[1]]+
  labs(y = "L1 survival")

ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/Validation/pxg_plot.pdf",
       width = 10,
       height = 6)

manplot(pr_maps)

ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/Validation/manplot.pdf",
       width = 12,
       height = 4)

# # # replicate data
p <- data.table::fread("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170103_set1-3.csv")%>%
  dplyr::select(strain, ratio = average, replicate, set)

pr_resid <- p%>%
  ungroup()%>%
  dplyr::do(augment(lm(ratio ~ set + replicate, data =.)))%>%
  dplyr::rename(resid = .resid)%>%
  data.frame()%>%
  dplyr::left_join(p, ., by = c("ratio","set","replicate"))%>%
  dplyr::distinct(ratio,set,replicate,strain,.keep_all = T) %>%
  dplyr::ungroup()%>%
  dplyr::select(strain, resid)%>%
  dplyr::group_by(strain)%>%
  dplyr::mutate(mph = median(resid, na.rm = T),
                sph = sd(resid, na.rm = T))%>%
  dplyr::filter(resid <= 2*sph+mph,
                resid >= mph-2*sph)%>%
  dplyr::mutate(m.p = mean(resid))%>%
  dplyr::arrange(desc(m.p))%>%
  dplyr::filter(strain %in% c("BRC20067", "DL238", "EG4725", "JU775","JU830"))

pr_resid$strain <- factor(pr_resid$strain,levels = pr_resid$strain,labels = pr_resid$strain, ordered = T)

ggplot(pr_resid) + 
  aes(y = resid, x =strain , fill = strain)+
  geom_boxplot()+
  geom_jitter(alpha = .5, width = .5)+
  scale_fill_manual(values = c( "#2474FF", "#EE8F03","cyan","pink", "red", "green"))+
  theme_bw()+
  theme(axis.text.x = element_text(size=12, face="bold", color="black"),
        axis.text.y = element_text(size=12, face="bold", color="black"),
        axis.title.x = element_text(size=14, face="bold", color="black"),
        axis.title.y = element_text(size=14, face="bold", color="black"),
        strip.text.x = element_text(size=20,face="bold", color="black"),
        strip.text.y = element_text(size=12, angle =0, face="bold", color="black"),
        plot.title = element_text(size=24, face="bold"),
        legend.position = "none")+
  labs(y = "Regressed L1 survival", y = "Strain")

ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/Validation/replicate_pheno_plot.pdf",
       width = 8,
       height = 4)

p%>%
  dplyr::filter(strain %in% c("BRC20067", "DL238", "EG4725", "JU775","JU830"))%>%
  dplyr::mutate(strain1 = factor(strain,levels = c( "DL238", "EG4725","BRC20067", "JU830", "JU775"),labels = c( "DL238", "EG4725","BRC20067", "JU830", "JU775"), ordered = T))%>%
  ggplot(.) + 
  aes(y = ratio, x =strain1 , fill = strain1)+
  geom_boxplot()+
  geom_jitter(alpha = .5, width = .5)+
  scale_fill_manual(values = c( "#2474FF", "#EE8F03","cyan","pink", "red", "green"))+
  theme_bw()+
  theme(axis.text.x = element_text(size=12, face="bold", color="black"),
        axis.text.y = element_text(size=12, face="bold", color="black"),
        axis.title.x = element_text(size=14, face="bold", color="black"),
        axis.title.y = element_text(size=14, face="bold", color="black"),
        strip.text.x = element_text(size=20,face="bold", color="black"),
        strip.text.y = element_text(size=12, angle =0, face="bold", color="black"),
        plot.title = element_text(size=24, face="bold"),
        legend.position = "none")+
  labs(y = "L1 survival", y = "Strain")

ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/Validation/non_regressed_replicate_pheno_plot.pdf",
       width = 8,
       height = 4)
