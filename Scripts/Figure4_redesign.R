library(data.table)
library(tidyverse)
library(cegwas)
library(cegwas2)
library(broom)
library(cowplot)
library(ggtree)
library(ggbeeswarm)
library(ggpubr)

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
setwd("..")

get_vcf2 <- function() {
  # Function for fetching the path the VCF
  path <- glue::glue("~/Dropbox/Andersenlab/Reagents/WormReagents/_SEQ/WI/WI-{cendr_dataset_release}/vcf/WI.{cendr_dataset_release}.snpeff.vcf.gz") # nolint
  if (file.exists(path)) {
    message("Using local VCF")
  } else {
    message("Using remote VCF")
    path <- glue::glue("http://storage.googleapis.com/elegansvariation.org/releases/{cendr_dataset_release}/variation/WI.{cendr_dataset_release}.soft-filter.vcf.gz") # nolint
  }
  path
}

#####################################################################################################################################################################
phen <- data.table::fread("Data/glct_causality/glct_deletion_expt_full.csv")

reg_pheno <- phen %>%
  dplyr::do(broom::augment(lm(Survival ~ Day, data =.)))%>%
  dplyr::select(-Survival) %>%
  dplyr::rename(Survival = .resid) %>%
  dplyr::mutate(Strain = phen$Strain) %>%
  dplyr::select(Strain, Survival) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Survival) %>%
  dplyr::mutate(norm_pheno_temp = ifelse(Survival == min(Survival), 0, 1))%>%
  dplyr::mutate(delta_pheno = ifelse(norm_pheno_temp == 0, 0, abs(dplyr::lag(Survival) - Survival)))%>%
  dplyr::mutate(norm_pheno = cumsum(delta_pheno)) %>%
  dplyr::mutate(final_pheno = norm_pheno/max(norm_pheno)) %>%
  dplyr::select(-Survival) %>%
  dplyr::rename(Survival = final_pheno)

l_compare <- list(c("BRC20067", "BRC20067_Del"), c("BRC20067", "BRC20067_Swap"), c("BRC20067", "DL238"))

pheno_plot <- reg_pheno%>%
  ggplot()+
  aes(x = Strain, y = Survival, fill = Strain)+
  geom_boxplot(alpha = 0.8, outlier.colour = NA)+
  scale_fill_manual(values=c("hotpink3","gray60","gray60","cadetblue3"))+
  scale_x_discrete(labels=c("BRC20067" = "BRC20067", 
                            "DL238" = "DL238", 
                            "BRC20067_Del" = expression(atop("BRC20067","Phe19fs")),
                            "BRC20067_Swap" = expression(atop("BRC20067","Gly16*"))))+
  ggbeeswarm::geom_beeswarm(cex = 0.7, priority = "density",  size = 0.5, alpha = 0.5)+
  theme_classic(12)+
  labs(y = "Normalized L1 survival") +
  scale_y_continuous(limits = c(0,1), expand = c(0.02, 0)) +
  theme(axis.title.x = element_blank(), 
        legend.position = "none") +
  stat_compare_means(comparisons = l_compare, label.y = c(0.75, 0.85, 0.05),tip.length =c(0.01,0.01,-0.01),
                     label = "p.signif", method = "t.test",
                     ref.group = "BRC20067", hide.ns = TRUE, size = 5, color = "gray50")

# effect sizes
p_effect <- sjstats::anova_stats(car::Anova(aov(Survival ~ Strain, data = dplyr::filter(reg_pheno, Strain %in% c("BRC20067", "DL238")))))[1,11]
d_effect <- sjstats::anova_stats(car::Anova(aov(Survival ~ Strain, data = dplyr::filter(reg_pheno, Strain %in% c("BRC20067", "BRC20067_Del")))))[1,11]
s_effect <- sjstats::anova_stats(car::Anova(aov(Survival ~ Strain, data = dplyr::filter(reg_pheno, Strain %in% c("BRC20067", "BRC20067_Swap")))))[1,11]
# BRC - del
d_effect/p_effect
# BRC - swap
s_effect/p_effect

#####################################################################################################################################################################
# load world map and remove antartica 
world <- map_data("world")
world <- world[world$region != "Antarctica",] # intercourse antarctica

# download data from CeNDR
isolation_info <- readr::read_tsv("https://elegansvariation.org/strain/strain_data.tsv")

cendr_dataset_release <- "20180527"

g3 <- cegwas2::query_vcf(c("glct-3"), vcf = get_vcf2()) 

alt_strains <-g3%>%
  dplyr::rowwise() %>%
  dplyr::filter( impact == "HIGH" & (a1==ALT || a2 == ALT || grepl(a1, ALT) ) | SAMPLE == "ECA369" | SAMPLE == "ECA347" | SAMPLE == "MY10") %>% 
  dplyr::filter(SAMPLE!="ECA701" & query == "glct-3") %>%
  dplyr::select(strain = SAMPLE, gene_name, aa_change, REF, ALT, a1, a2) %>%
  dplyr::mutate(a_col = ifelse((a1 == "T" | a2 == "T") & gene_name == "glct-3", "ALT", 
                               ifelse(a1 == ALT | a2 == ALT,"ALT","REF")) ) %>%
  dplyr::filter(a_col != "REF") %>%
  dplyr::distinct(strain, .keep_all=T)



strains_330 <- isolation_info%>%
  dplyr::filter(reference_strain == "True")%>%
  dplyr::select(strain = isotype, long = longitude, lat = latitude) %>%
  dplyr::left_join(., alt_strains, by = "strain")  


map <- ggplot()+ 
  geom_map(data=world, map=world,
           aes(x=long, y=lat, map_id=region),
           color="black", fill="white", size=0.25)+
  geom_point(data = strains_330 %>% dplyr::filter(is.na(aa_change)) , 
             aes(x=as.numeric(long), y=as.numeric(lat)), 
             fill = "hotpink3",
             shape = 21, 
             size = 1, stroke = 0.1) +
  geom_point(data = strains_330 %>% dplyr::filter(!is.na(aa_change)), 
             aes(x=as.numeric(long), y=as.numeric(lat)), 
             fill = "cadetblue3",
             shape = 21, 
             size = 1.5, stroke = 0.1) +
  scale_fill_manual(values = c("cadetblue3","hotpink3"))+
  theme_map()+
  scale_y_continuous(expand=c(0.01,0.01))+
  scale_x_continuous(expand=c(0.01,0.01))

ggsave(plot = map, "Plots/SVG_PLOTS/FigureSXX_MAP.pdf", height = 4, width = 6.5, units = "in")
ggsave(plot = map, "Plots/SVG_PLOTS/FigureSXX_MAP.png", height = 4, width = 6.5, units = "in", dpi = 300)

##################################################################################################################################################################### GWAS pheno plot, by glct allele
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
  dplyr::mutate(st_colors = ifelse(strain %in% alt_strains$strain, "glct3", "ref")) 


gwas_split <- ggplot(pr_resid_pldf) +
  aes(y = mph, x = st_colors) +
  geom_boxplot(outlier.colour = NA) +
  geom_beeswarm(data = pr_resid_pldf %>% dplyr::filter(!strain %in% c("DL238", "BRC20067")), size = 0.5, alpha = 0.5) +
  geom_point(data = pr_resid_pldf %>% dplyr::filter(strain %in% c("DL238", "BRC20067")), shape = 21, size = 3, aes(fill = strain)) +
  scale_fill_manual(values = c("DL238"  = "cadetblue3", "N2"  = "orange", "BRC20067" = "hotpink3", "not" = "gray70")) +
  theme_classic(12) +
  scale_x_discrete(labels=c("ref" = "GLCT-3", 
                            "glct3" = expression(atop("GLCT-3","Gly16*"))))+
  # scale_x_discrete(labels=c("ref" = "GLCT-3", 
  #                           "glct3" = expression(""~ Delta~"GLCT-3"))) + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), legend.title = element_blank()) +
  scale_y_continuous(limits = c(0,1), expand = c(0.02, 0)) +
  labs(x = "Strain", y = "Normalized L1 survival") 

##################################################################################################################################################################### GWAS fine mapping
finemap <- data.table::fread("Data/chr1_finemap/value.I.11078010.13929469_prLD_df.tsv")

load("/Users/Stefan/github_repos/Propionate/vcf/gene_ref_flat.Rda")

gene_df <- gene_ref_flat %>%
  dplyr::filter(wbgene %in% c("WBGene00011781")) %>%
  dplyr::select(gene_id = wbgene, strand, txstart, txend, feature_id = gene) %>%
  dplyr::arrange(txstart, feature_id)%>%
  dplyr::distinct(gene_id, feature_id, .keep_all = TRUE) %>%
  dplyr::distinct(gene_id, peak_marker, CHROM, strand, txstart, txend, start_pos, end_pos) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(gene_name = ifelse(gene_id == "WBGene00011779", "bgnt-1.7",
                                   ifelse(gene_id == "WBGene00001642", "gly-17",
                                          ifelse(gene_id == "WBGene00001641", "gly-16",
                                                 ifelse(gene_id == "WBGene00011780", "T15D6.5",
                                                        ifelse(gene_id == "WBGene00003667", "nhr-77",
                                                               ifelse(gene_id == "WBGene00011781", "glct-3",NA)))))),
                plot_pos = mean(c(txstart,txend)))



fmap_plt <- ggplot(finemap) +
  geom_rect(aes(xmin = ifelse(strand == "+", txstart/1e6, txend/1e6),
                xmax = ifelse(strand == "+", txend/1e6, txstart/1e6),
                ymin = 0,
                ymax = max(finemap$value)), size = 1, data = gene_df, alpha = 0.25) +
  ggplot2::geom_point( aes(y = value, x = POS/1e6), size = 0.5 ) +
  # geom_segment(aes(x = ifelse(strand == "+", txstart/1e6, txend/1e6),
  #                  xend = ifelse(strand == "+", txend/1e6, txstart/1e6),
  #                  y = -0.02,
  #                  yend =  -0.02),
  #              arrow = arrow(length = unit(5, "points")), size = 1, data = gene_df)  +
  ggplot2::theme_bw(12) +
  xlim(12.3,12.5) +
  scale_y_continuous(limits = c(-0.02,5.1), expand = c(0, 0.1))+
  ggplot2::labs(x = "Genomic position (Mb)",
                y = expression(-log[10](italic(p))))
fmap_plt

#####################################################################################################################################################################
# load tree
tree <- ape::read.tree(glue::glue("Data/whole_genome_tree/330_genome.raxml.bestTree"))

# highlight branches for strains of interest
branch_strains <- list(CONNECT = alt_strains$strain)

tree_pt_h <- ggtree::groupOTU(tree, branch_strains)

ggtree(tree_pt_h,
                    branch.length="rate", 
                    aes(color=group), size = 0.25) + 
  scale_color_manual(values=c("hotpink3", "cadetblue3"), 
                     name = "GLCT-3\nallele", 
                     labels=c("REF", "Gly16*")) + 
  theme(legend.position="right")+
  theme_tree2() + 
  scale_y_continuous(expand=c(0.01,0.01))

ggsave( "Plots/SVG_PLOTS/FigureSXX_GENOME_TREE.pdf", height = 10, width = 6.5, units = "in")
ggsave( "Plots/SVG_PLOTS/FigureSXX_GENOME_TREE.png", height = 10, width = 6.5, units = "in", dpi = 300)

##################################################################################################################################################################### gene_tree

# glct3_gene <- cegwas2::query_vcf("WBGene00011781", vcf = "~/UCLA/Genomics_Data/VCFs/328_Soft_csq_annotated.vcf.gz")
# 
# glct_gene_matrix <- glct3_gene %>%
#   dplyr::filter(FILTER=="PASS") %>%
#   dplyr::select(POS, SAMPLE, ALT, a1,a2) %>%
#   dplyr::mutate(nGT = ifelse(a1 == ALT | a2 == ALT, 1, 0)) %>%
#   dplyr::select(POS, SAMPLE, nGT) %>%
#   dplyr::distinct(POS, SAMPLE, .keep_all=T) %>%
#   tidyr::spread(POS, nGT)
# 
# glct_gene_matrix$altct <- rowSums(glct_gene_matrix[,2:ncol(glct_gene_matrix)], na.rm = T)
# 
# glct_gene_matrix <- glct_gene_matrix %>% dplyr::filter(altct > 0 | SAMPLE %in% c("BRC20067", "N2"))
# 
# row.names(glct_gene_matrix) <- glct_gene_matrix$SAMPLE
# 
# glct_gene_tree = ape::nj(ape::dist.gene(glct_gene_matrix, pairwise.deletion = T,method = "percentage"))
# 
# 
# ggtree(glct_gene_tree) +
#   geom_tiplab(size = 2) +
#   theme_tree()+
#   theme(axis.line.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_blank()) 



top_panel <- cowplot::plot_grid(fmap_plt,
                                labels = c("A"))

mid_panel <- cowplot::plot_grid(gwas_split + theme(legend.position = "none"),
                                pheno_plot + theme(axis.title.y = element_blank(), axis.text.y = element_blank()), 
                                nrow = 1, 
                                rel_widths = c(2,4),
                                label_size = 14,
                                labels = c("B","C"), align = "hv")

fig4 <- cowplot::plot_grid(top_panel,
                           mid_panel, 
                           nrow = 2, 
                           label_size = 14,
                           align = "hv",
                           axis = "left", rel_heights = c(0.7,1))



ggsave(plot = fig4, "Plots/SVG_PLOTS/Figure4.pdf", height = 5, width = 6.5, units = "in")
ggsave(plot = fig4, "Plots/SVG_PLOTS/Figure4.png", height = 5, width = 6.5, units = "in", dpi = 300)
