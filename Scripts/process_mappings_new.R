library(data.table)
library(tidyverse)
library(cegwas)
library(cegwas2)
library(broom)
library(cowplot)

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
setwd("..")

source("Scripts/Figure_Functions.R")

p <- data.table::fread("Data/20161128_new_gwas_replicates/20170103_set1-3.csv")%>%
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

pheno_plot <- ggplot(pr_resid_pldf) +
  aes(y = mph, x = strain, fill = st_colors) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=mph, ymax=mph+sph), width=.4,
                position=position_dodge(.9)) +
  scale_fill_manual(values = c("DL238"  = "cadetblue3", "N2"  = "orange", "BRC20067" = "hotpink3", "not" = "gray70")) +
  theme_classic(20) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  labs(x = "Strain", y = "Normalized\nL1 Survival")

ggsave(pheno_plot, filename = "Plots/GWAS_Propionate_Phenotypes.pdf", height = 4, width = 10)
ggsave(pheno_plot, filename = "Plots/GWAS_Propionate_Phenotypes.png", height = 4, width = 10, dpi = 300)
ggsave(pheno_plot, filename = "Plots/GWAS_Propionate_Phenotypes.svg", height = 4, width = 10)

p_reg <- pr_resid%>%
  dplyr::summarise(phenotype = mean(resid, na.rm = T))

pr_pheno <- process_pheno(p_reg)
maps <- gwas_mappings(pr_pheno)

pr_maps <- process_mappings(maps,pr_pheno)

save(pr_maps, file ="Processed_Data/GWAS_PRmaps.Rda")

# save pheno for cegwas2 pipeline

pa_pheno <- na.omit(pr_maps) %>%
  dplyr::select(strain, value)

write.table(pa_pheno, file ="Processed_Data/GWAS_PRpheno.tsv", col.names = T, row.names = F, quote = F, sep = "\t")

# load cegwas2 processed mapping

c2_prmaps <- data.table::fread("Processed_Data/GWAS_processed_mapping_cegwas2.tsv")

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
                     color = "blue",fill = "cyan",linetype = 2, 
                     alpha=.3)+
  ggplot2::geom_hline(ggplot2::aes(yintercept = BF),
                      color = "gray60", 
                      alpha = .75,  
                      size = 1) +
  ggplot2::geom_hline(ggplot2::aes(yintercept = EIGEN_CUTOFF),
                      color = "gray60", 
                      alpha = .75,  
                      size = 1,
                      linetype = 2) +
  ggplot2::geom_point( ggplot2::aes(color= factor(EIGEN_SIG), alpha = factor(EIGEN_SIG)) ) +
  ggplot2::facet_grid( . ~ CHROM, scales = "free_x" , space = "free_x") +
  ggplot2::theme_bw(20) +
  scale_y_continuous(limits = c(0,10), expand = c(0, 0)) +
  ggplot2::theme(strip.background = element_blank(),
                 legend.position = "none",
                 panel.grid.minor = element_blank(),
                 panel.grid.major = element_blank()) +
  ggplot2::labs(x = "Genomic Position (Mb)",
                y = expression(-log[10](italic(p))))

ggsave("Plots/GWAS_Propionate_rrblup.pdf", 
       plot = rrblup_map,
       height = 4,
       width = 12)

ggsave("Plots/GWAS_Propionate_rrblup.png", 
       plot = rrblup_map,
       height = 4,
       width = 12,
       dpi = 300)

ggsave("Plots/GWAS_Propionate_rrblup.svg", 
       plot = rrblup_map,
       height = 4,
       width = 12)

# # # # # # # # # # # # # # # BURDEN MAPPING
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
  geom_point()+
  scale_color_manual(values=c("black","red"))+
  facet_grid(.~CHROM, scales = "free", space = "free")+
  ggplot2::theme_bw(18) +
  ggplot2::theme(strip.background = element_blank(),
                 legend.position = "none",
                 panel.grid.minor = element_blank(),
                 panel.grid.major = element_blank()) +
  scale_y_continuous(limits = c(0,10), expand = c(0, 0)) +
  labs(x = "Genomic Position (Mb)", 
       y = expression(-log[10](italic(p))))

# figure 2 cowplot
cowplot::plot_grid(pheno_plot,
                   rrblup_map+ theme(axis.text.x = element_blank(),
                                     axis.title.x = element_blank()),
                   skat_plot + theme(strip.text.x = element_blank()),
                   ncol =1, 
                   label_size = 20, 
                   rel_heights = c(0.8,1,1),
                   align = "v", axis = "l",
                   labels = "AUTO")

ggsave("Plots/Figure2.png", 
       height = 12,
       width = 12,
       dpi = 300)

ggsave("Plots/Figure2.svg", 
       height = 10,
       width = 12,
       dpi = 300)

# # # # # # # # # # # # # # PLOT QTL LD

gm <- readr::read_tsv(glue::glue("Processed_Data/Genotype_Matrix.tsv"))

LD_output <- Plot_Peak_LD(c2_prmaps, gm)

LD_output[[1]] + 
  theme_bw(18) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

ggsave(filename = "Plots/GWAS_Peak_LD.pdf", height = 8, width = 12)
ggsave(filename = "Plots/GWAS_Peak_LD.png", height = 8, width = 12, dpi = 300)
ggsave(filename = "Plots/GWAS_Peak_LD.svg", height = 8, width = 12)

# # # # # # # # # # # # # # all gene pxg

genes <- query_vcf(c("I:12374204-12388791"))

norm_pheno <- na.omit(c2_prmaps) %>%
  dplyr::distinct(strain, value) %>%
  dplyr::ungroup()%>%
  dplyr::arrange(value) %>%
  dplyr::mutate(norm_pheno_temp = ifelse(value == min(value), 0, 1))%>%
  dplyr::mutate(delta_pheno = ifelse(norm_pheno_temp == 0, 0, abs(dplyr::lag(value) - value)))%>%
  dplyr::mutate(norm_pheno = cumsum(delta_pheno)) %>%
  dplyr::mutate(final_pheno = norm_pheno/max(norm_pheno)) %>%
  dplyr::select(strain, final_pheno)

gene_pr <- dplyr::filter(genes, SAMPLE %in% unique(c2_prmaps$strain)) %>%
  dplyr::select(strain = SAMPLE, REF, ALT, gene_id, gene_name, allele = a1, impact, effect, aa_change) %>%
  dplyr::mutate(TGT = ifelse(allele == ALT, "ALT", "REF")) %>%
  dplyr::left_join(., norm_pheno, by = "strain") %>%
  dplyr::group_by(strain, gene_id) %>%
  dplyr::mutate(vt_ct = sum(grepl("ALT",TGT)))

gene_ref <- dplyr::filter(gene_pr,  vt_ct == 0) %>%
  dplyr::distinct(strain, final_pheno, gene_id, .keep_all = T) %>%
  dplyr::mutate(plot_gt = "REF")

gene_high <- dplyr::filter(gene_pr, impact == "HIGH", TGT == "ALT") %>%
  dplyr::group_by(gene_id) %>%
  dplyr::arrange(impact) %>%
  dplyr::filter(impact == "HIGH") %>%
  dplyr::mutate(plot_gt = effect)

gene_alt <- dplyr::filter(gene_pr, TGT == "ALT") %>%
  # dplyr::filter(!strain %in% gene_high$strain | strain %in% c("BRC20067", "DL238")) %>%
  dplyr::group_by(strain, gene_id) %>%
  dplyr::mutate(plot_gt = paste(unique(aa_change), collapse = " ")) %>%
  dplyr::mutate(plot_gt = gsub("p\\.","", plot_gt)) %>%
  dplyr::distinct(strain, final_pheno, gene_id,plot_gt, .keep_all =T)

sig_genes <- skat_maps%>%
  dplyr::filter(significant == T, grepl("I:",RANGE))%>% dplyr::pull(Gene)

gene_all <- dplyr::bind_rows(list(gene_ref,gene_high,gene_alt)) %>%
  dplyr::filter(gene_id %in% sig_genes)



countSpaces <- function(s) { sapply(gregexpr(" ", s), function(p) { sum(p>=0) } ) }

gene_plots <- list()
for(g in 1:length(unique(gene_all$gene_id))){
  
  temp_all <- gene_all %>% 
    dplyr::filter(gene_id == unique(gene_all$gene_id)[g]) %>%
    dplyr::ungroup()
  
  temp_all$ids <- temp_all %>%
    dplyr::group_by(plot_gt) %>%
    group_indices(., plot_gt) 
  
  temp_all <- temp_all %>%
    dplyr::rowwise() %>%
    dplyr::mutate(alt_ct = countSpaces(plot_gt)) %>%
    dplyr::mutate(alt_ct = ifelse(plot_gt!="REF" & !plot_gt %in% gene_high$plot_gt, paste0("ALT_CT:",alt_ct+1," Hap:", ids), plot_gt)) %>%
    dplyr::mutate(alt_ct = gsub("\\&","\n",alt_ct)) %>%
    dplyr::group_by(alt_ct, gene_id) %>%
    dplyr::mutate(mean_gt = median(final_pheno))
  
  temp_df <- temp_all %>%
    dplyr::ungroup()%>%
    dplyr::arrange(desc(mean_gt)) %>%
    dplyr::mutate(pt_gt = factor(alt_ct, levels = unique(c(alt_ct[!alt_ct %in% "REF"], "REF")))) %>%
    dplyr::mutate(pt_col = factor(strain, 
                                  levels = c("BRC20067", "DL238", unique(strain)[!unique(strain) %in% c("BRC20067", "DL238")]) ,
                                  labels = c("BRC20067", "DL238", rep("other", length(unique(strain)[!unique(strain) %in% c("BRC20067", "DL238")]))))) %>%
    dplyr::distinct(strain, final_pheno, alt_ct,pt_col,.keep_all = T)
  
  gene_plots[[g]] <- ggplot(temp_df)+
    ggbeeswarm::geom_beeswarm(aes(x =pt_gt, y = final_pheno, 
                                  fill = pt_col, shape = pt_col, size = pt_col))+
    scale_fill_manual(values = c("hotpink3", "cadetblue3", "gray50"))+
    scale_size_manual(values = c(5, 5, 2))+
    scale_shape_manual(values = c(23, 23, 21))+
    geom_point(aes(x = pt_gt, y = mean_gt), 
               data = dplyr::distinct(temp_df, alt_ct, gene_id, mean_gt,.keep_all=T),
               shape = 23, 
               size = 3,
               fill = "red")+
    theme_classic(18) +
    facet_wrap(~gene_name, scale = "free_y") +
    labs(x = "Genotype",
         y = "Normalized L1 Survival") +
    theme(legend.position = "none",strip.text = element_text(face = "italic"))+
    coord_flip()
}

length(gene_plots)

all_gene_plot <- cowplot::plot_grid(gene_plots[[1]]+ theme(axis.title.x = element_blank()) ,
                   gene_plots[[2]]+ theme(axis.title.y = element_blank(),axis.title.x = element_blank()),
                   gene_plots[[3]],
                   gene_plots[[4]]+theme(axis.title.y = element_blank()), 
                   ncol = 2, 
                   label_size = 20, align = "vh", labels = "AUTO")
  
ggsave(plot = all_gene_plot, filename = "Plots/chr1_gene_geno_pheno.pdf", height = 12, width = 20)
ggsave(plot = all_gene_plot, filename = "Plots/chr1_gene_geno_pheno.png", height = 12, width = 20, dpi = 300)
ggsave(plot = all_gene_plot, filename = "Plots/chr1_gene_geno_pheno.svg", height = 12, width = 20)

# unique haplotypes that contain stop gain in glct-e
glct_stop_hap <- dplyr::filter(temp_df, plot_gt!="REF", grepl("Gly16*", plot_gt)) %>% dplyr::pull(alt_ct) %>% unique()

dplyr::filter(temp_df, alt_ct %in% glct_stop_hap)  %>% dplyr::pull(plot_gt) %>% unique()

# # # # # # # # # # # # # # GLCT-3

glct3 <- query_vcf("glct-3")

pa_pheno <- na.omit(c2_prmaps) %>%
  dplyr::distinct(strain, value)

glct3_pr <- dplyr::filter(glct3, SAMPLE %in% unique(c2_prmaps$strain)) %>%
  dplyr::select(strain = SAMPLE, REF, ALT, allele = a1, impact, effect, aa_change) %>%
  dplyr::mutate(TGT = ifelse(allele == ALT, "ALT", "REF")) %>%
  dplyr::left_join(., pa_pheno, by = "strain") %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(vt_ct = length(unique(TGT)))

glct_ref <- dplyr::filter(glct3_pr, TGT == "REF", vt_ct == 1) %>%
  dplyr::distinct(strain, value, .keep_all = T) %>%
  dplyr::mutate(plot_gt = "REF")

glct_high <- dplyr::filter(glct3_pr, TGT == "ALT") %>%
  dplyr::arrange(impact) %>%
  dplyr::filter(impact == "HIGH") %>%
  dplyr::mutate(plot_gt = "Stop Gained")

glct_alt <- dplyr::filter(glct3_pr, TGT == "ALT") %>%
  dplyr::arrange(impact) %>%
  dplyr::filter(!strain %in% glct_high$strain) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(plot_gt = paste(unique(aa_change), collapse = " ")) %>%
  dplyr::mutate(plot_gt = gsub("p\\.","", plot_gt)) %>%
  dplyr::distinct(strain, value, .keep_all =T)

glct_all <- dplyr::bind_rows(list(glct_ref,glct_high,glct_alt))

plot_df <- glct_all %>%
  dplyr::ungroup()%>%
  dplyr::arrange(value) %>%
  dplyr::mutate(norm_pheno_temp = ifelse(value == min(value), 0, 1))%>%
  dplyr::mutate(delta_pheno = ifelse(norm_pheno_temp == 0, 0, abs(dplyr::lag(value) - value)))%>%
  dplyr::mutate(norm_pheno = cumsum(delta_pheno)) %>%
  dplyr::mutate(final_pheno = norm_pheno/max(norm_pheno)) %>%
  dplyr::group_by(plot_gt) %>%
  dplyr::mutate(mean_gt = median(final_pheno)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(plot_gt = factor(plot_gt, levels = c("REF", "Ile35Phe Leu36Ser", "Asp53Glu", "Stop Gained"))) %>%
  dplyr::mutate(pt_col = ifelse(strain == "BRC20067", "br",
                                ifelse(strain == "DL238", "dl", "other")))

glct_var_pt <- ggplot(plot_df)+
  ggbeeswarm::geom_beeswarm(aes(x = plot_gt, y = final_pheno, 
                                fill = pt_col, shape = pt_col, size = pt_col))+
  scale_fill_manual(values = c("hotpink3", "cadetblue3", "gray50"))+
  scale_size_manual(values = c(5, 5, 2))+
  scale_shape_manual(values = c(23, 23, 21))+
  geom_point(aes(x = plot_gt, y = mean_gt), 
             data = dplyr::distinct(plot_df, plot_gt,mean_gt),
             shape = 23, 
             size = 3,
             fill = "red")+
  theme_classic(20) +
  labs(x = expression(paste(italic("glct-3"), " Genotype")),
       y = "Normalized L1 Survival") +
  theme(legend.position = "none")

ggsave(plot = glct_var_pt, filename = "Plots/glct3_geno_pheno.pdf", height = 6, width = 10)
ggsave(plot = glct_var_pt, filename = "Plots/glct3_geno_pheno.png", height = 6, width = 10, dpi = 300)
ggsave(plot = glct_var_pt, filename = "Plots/glct3_geno_pheno.svg", height = 6, width = 10)
# 
# glct_3_crispr <- data.table::fread("Data/glct_deletion_expt.csv") %>%
#   dplyr::select(BRC20067, DL238, glct3, Day) %>%
#   tidyr::gather(Strain, pheno, -Day) 
# 
# glct_plot <- glct_3_crispr %>%
#   ggplot()+
#   aes(x = Strain, y = pheno, fill = Strain) +
#   geom_boxplot() +
#   ggbeeswarm::geom_beeswarm() +
#   scale_fill_manual(values = c("hotpink3", "cadetblue3", "gray70"))+
#   labs(y = "L1 Survival") +
#   theme_classic(20) +
#   scale_x_discrete(labels=c("BRC20067" = "BRC20067", "DL238" = "DL238", "glct3" = expression("BRC20067 "~ Delta~italic("glct-3"))))+
#   theme(axis.title.x = element_blank(),
#         legend.position = "none")
# 
# ggsave(plot = glct_plot, filename = "Plots/glct3_crispr.pdf", height = 6, width = 10)
# ggsave(plot = glct_plot, filename = "Plots/glct3_crispr.png", height = 6, width = 10, dpi = 300)
# 
# plots <- align_plots(glct_var_pt, skat_plot, align = 'v', axis = 'l')
# 
# bottom_plot <- cowplot::plot_grid(plots[[1]],
#                                   glct_plot,
#                    ncol =2, 
#                    label_size = 20, 
#                    align = "h", rel_widths = c(1,0.8),
#                    labels = c('B','C'))
# 
# cowplot::plot_grid(plots[[2]],
#                    bottom_plot,
#                    ncol = 1,
#                    label_size = 20, 
#                    labels = c("A",""))
# 
# 
# ggsave(filename = "Plots/Figure4.pdf", height = 10, width = 14)
# ggsave(filename = "Plots/Figure4.png", height = 10, width = 14, dpi = 300)
# 
# pr_genes <- genes%>%
#   dplyr::filter(SAMPLEW == "CB4856")%>%
#   dplyr::distinct(CHROM,POS,aa_change, .keep_all = T)
# 
# save(pr_genes, file ="~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/pr_genes.Rda")
# 
# ggplot(pr_genes)+
#   aes(x = POS/1e6, y = abs_spearman_cor, fill = GT)+
#   geom_point(shape=21)+
#   scale_fill_manual(values = c("ALT" = "cyan", "REF" = "gray"))+
#   scale_color_manual(values = c("ALT" = "cyan", "REF" = "gray"))+
#   ggplot2::theme_bw() + 
#   ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16, face = "bold", color = "black"),
#                  axis.text.y = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
#                  axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
#                  axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
#                  strip.text.x = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
#                  strip.text.y = ggplot2::element_text(size = 16,  face = "bold", color = "black"), 
#                  plot.title = ggplot2::element_text(size = 0,  face = "bold", vjust = 1),
#                  panel.background = ggplot2::element_rect(color = "black",  size = 1.2), 
#                  strip.background = ggplot2::element_rect(color = "black", size = 1.2),
#                  legend.position = 'none')+
#   labs(x = "Genomic Position (Mb)", y = expression(bold(paste("Spearman's ", bolditalic("rho")))))
# 
# 
# 
# 
# hist <-ggplot(p_reg)+
#   aes(x = phenotype)+
#   geom_histogram(color = "black", fill = "gray50")+
#   # scale_fill_manual(values = c("ALT" = "cyan", "REF" = "gray"))+
#   # scale_color_manual(values = c("ALT" = "cyan", "REF" = "gray"))+
#   ggplot2::theme_bw() + 
#   ggplot2::theme(axis.text.x = ggplot2::element_text(size = 14, face = "bold", color = "black"),
#                  axis.text.y = ggplot2::element_text(size = 14, face = "bold", color = "black"), 
#                  axis.title.x = ggplot2::element_text(size = 16, face = "bold", color = "black", vjust = -0.3), 
#                  axis.title.y = ggplot2::element_text(size = 18, color = "black"), 
#                  strip.text.x = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
#                  strip.text.y = ggplot2::element_text(size = 16,  face = "bold", color = "black"), 
#                  plot.title = ggplot2::element_text(size = 0,  face = "bold", vjust = 1),
#                  panel.background = ggplot2::element_rect(color = "black",  size = 1.2), 
#                  strip.background = ggplot2::element_rect(color = "black", size = 1.2),
#                  legend.position = 'none')+
#   labs(x = "Regressed Propionate Survival", y = "Count")
# 
# load( file ="~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/3set_pr_maps.Rda")
# 
# mp <- manplot(pr_maps)
# 
# mp <- mp[[1]] +  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 14, face = "bold", color = "black"),
#                                 axis.text.y = ggplot2::element_text(size = 14, face = "bold", color = "black"), 
#                                 axis.title.x = ggplot2::element_text(size = 16, face = "bold", color = "black", vjust = -0.3), 
#                                 axis.title.y = ggplot2::element_text(size = 18, face = "bold", color = "black"), 
#                                 strip.text.x = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
#                                 strip.text.y = ggplot2::element_text(size = 16,  face = "bold", color = "black"), 
#                                 plot.title = ggplot2::element_text(size = 0,  face = "bold", vjust = 1),
#                                 panel.background = ggplot2::element_rect(color = "black",  size = 1.2), 
#                                 strip.background = ggplot2::element_rect(color = "black", size = 1.2),
#                                 legend.position = 'none')
# 
# 
# 
# 
# ggdraw() +
#   draw_plot(hist, 0, 0, .2667, .925) +
#   draw_plot(mp, .2667, 0, 1-.2667, 1) +
#   draw_plot_label(c("A", "B"), c(0, .2667), c(1, 1), size = 16)
# 
# cowplot::ggsave(filename = "~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/grant_figure.png",
#                 plot = last_plot(),
#                 width = 15,
#                 height = 4)


