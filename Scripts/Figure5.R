# generate data
# Rscript --vanilla Scripts/Interval_Popgen.R I 12274204 12488791 ../Ce330_annotated.vcf.gz ../WS245_exons.gff Propionate ../ce330_strains.txt

# negative fay an wu h = excess of high frequency derived alleles
# Rscript --vanilla Scripts/Interval_Popgen.R I 1 15072434 ../Ce330_annotated.vcf.gz ../WS245_exons.gff Propionate ../ce330_strains.txt

# glct6
# Rscript --vanilla Scripts/Interval_Popgen.R IV 3684358 3894044 ../Ce330_annotated.vcf.gz ../WS245_exons.gff Propionate ../ce330_strains.txt

# all chrom 4
# Rscript --vanilla Scripts/Interval_Popgen.R IV 1 17493829 ../Ce330_annotated.vcf.gz ../WS245_exons.gff Propionate ../ce330_strains.txt


library(tidyverse)
library(ggtree)

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
setwd("..")

source("Scripts/Figure_Themes.R")

load("Propionate_I_1-15072434_Diversity_Statistics.Rda")

qtl_start <- 12374204
qtl_end <- 12388791
glct3 <- c(12385766, 12388791)
glct1 <- c(12337968, 12339239)
glct2 <- c(12435726, 12437285)
glct4 <- c(12317321, 12318735)
glct5 <- c(12453389, 12454764)
glct6 <- c(3784237, 3794152)

td_df <- d2_df %>%
  dplyr::filter(statistic %in% c("Pi")) %>%
  dplyr::group_by(statistic) %>%
  dplyr::mutate(scaled_value = scale(value)) %>%
  dplyr::mutate(q10 = quantile(value, 0.99, na.rm = T)) %>%
  dplyr::mutate(outlier = ifelse(value > q10, "Low", "not"))

chrom1 <- td_df %>%
  ggplot()+
  geom_point(size = 0.5)+
  aes(x = startWindow/1e6, y = value, color = outlier)+
  geom_vline(aes(xintercept = glct3[1]/1e6), color = "#E68FAC")+
  geom_vline(aes(xintercept = glct3[2]/1e6), color = "#E68FAC")+
  geom_vline(aes(xintercept = glct1[1]/1e6), color = "gray80")+
  geom_vline(aes(xintercept = glct1[2]/1e6), color = "gray80")+
  geom_vline(aes(xintercept = glct2[1]/1e6), color = "gray60")+
  geom_vline(aes(xintercept = glct2[2]/1e6), color = "gray60")+
  geom_vline(aes(xintercept = glct4[1]/1e6), color = "gray40")+
  geom_vline(aes(xintercept = glct4[2]/1e6), color = "gray40")+
  geom_vline(aes(xintercept = glct5[1]/1e6), color = "gray20")+
  geom_vline(aes(xintercept = glct5[2]/1e6), color = "gray20")+
  scale_color_manual(values = c(highlight_color, "#222222"))+
  # facet_grid(statistic~., scales = "free")+
  # xlim(12,13)+
  theme_bw(18) +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  labs(x = "Genomic Position (Mb)",
       y = (expression(italic(theta["w"]))))

region <-  td_df %>%
  ggplot()+
  aes(x = startWindow/1e6, y = value, color = outlier)+
  geom_point(size = 0.5)+
  geom_vline(aes(xintercept = glct3[1]/1e6), color = "#E68FAC", size = 0.5)+
  # geom_vline(aes(xintercept = glct3[2]/1e6), color = "#E68FAC")+
  geom_vline(aes(xintercept = glct1[1]/1e6), color = "gray80", size = 0.5)+
  # geom_vline(aes(xintercept = glct1[2]/1e6), color = "gray80")+
  geom_vline(aes(xintercept = glct2[1]/1e6), color = "gray60", size = 0.5)+
  # geom_vline(aes(xintercept = glct2[2]/1e6), color = "gray60")+
  geom_vline(aes(xintercept = glct4[1]/1e6), color = "gray40", size = 0.5)+
  # geom_vline(aes(xintercept = glct4[2]/1e6), color = "gray40")+
  geom_vline(aes(xintercept = glct5[1]/1e6), color = "gray20", size = 0.5)+
  # geom_vline(aes(xintercept = glct5[2]/1e6), color = "gray20")+
  scale_color_manual(values = c(highlight_color, "#222222"))+
  # facet_grid(statistic~., scales = "free")+
  xlim(12,13)+
  theme_bw(12) +
  theme(legend.position = "none") +
  labs(x = "Genomic Position (Mb)",
       y = (expression(italic(theta["w"]))))


load("Propionate_IV_1-17493829_Diversity_Statistics.Rda")

glct6_r <- d2_df %>%
  dplyr::filter(statistic %in% c("Pi")) %>%
  dplyr::group_by(statistic) %>%
  dplyr::mutate(scaled_value = scale(value)) %>%
  dplyr::mutate(q10 = quantile(value, 0.99, na.rm = T)) %>%
  dplyr::mutate(outlier = ifelse(value > q10, "Low", "not")) %>%
  ggplot()+
  aes(x = startWindow/1e6, y = value, color = outlier)+
  geom_point(size = 0.5)+
  geom_vline(aes(xintercept = glct6[1]/1e6), color = "gray20", size = 0.5)+
  # geom_vline(aes(xintercept = glct6[2]/1e6), color = "gray20")+
  scale_color_manual(values = c(highlight_color, "#222222"))+
  # facet_grid(statistic~., scales = "free")+
  theme_bw(12) +
  xlim(3.5,4.5)+
  theme(legend.position = "none") +
  labs(x = "Genomic Position (Mb)",
       y = (expression(italic(theta["w"]))))


#####################################################################################################################################################################
# load nucleotide tree

n_tree <- treeio::read.raxml("Data/phylo/Elegans_only/RAxML_bipartitionsBranchLabels.elegans_partition")
# 
# original <- n_tree@phylo$tip.label
# clean_names <- c("glct-5", "glct-6", "glct-6", "glct-6","glct-1","glct-2","glct-3","glct-4")
# rename_df <- data.frame(label = original, newlab = clean_names)
# 
# n_tree_clean <- treeio::rename_taxa(n_tree, rename_df, label, newlab)

n_t <- ggtree(n_tree) + 
  scale_color_viridis_c(direction = -1) +
  theme(legend.position="right")+
  theme_classic(12) +
  scale_x_continuous(limits = c(0,2),breaks = c(0,0.5,1, 1.5, 2)) +
  geom_tiplab(aes(label=paste0('~italic(', label, ')')), parse=T, hjust = 0, size = 2) +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) 


# cp <- collapse(n_t, node=11)
# cp + geom_cladelabel(node=11, "A")
# 
# n_tree_final <- cp + geom_label2(aes(subset=(node == 11), label=paste0('~italic(', "glct-6", ')')), parse = T, label.size = NA) +
#   xlim(0,1.1)+ geom_label2(aes(label=bootstrap, subset=(node != 11)))+
#   theme_classic(12) +
#   theme(axis.line.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_blank()) 
                 
# load protein tree
p_tree <- treeio::read.raxml("Data/phylo/Protein/RAxML_bipartitionsBranchLabels.elegans_partition_paralogs")

# original <- p_tree@phylo$tip.label
# clean_names <- c("GLCT-2", "GLCT-3", "GLCT-6", "GLCT-4","GLCT-5","GLCT-1")
# rename_df <- data.frame(label = original, newlab = clean_names)
# 
# p_tree_clean <- treeio::rename_taxa(p_tree, rename_df, label, newlab)

p_t <- ggtree(p_tree) + 
  scale_color_viridis_c(direction = -1) +
  theme(legend.position="right")+
  theme_classic(12) +
  geom_tiplab(size = 2) +
  xlim(0,2.5)+ 
  # geom_label2(aes(label=bootstrap), size = 2) +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) 

# top_plots <- cowplot::plot_grid(n_tree_final, p_t, labels = "AUTO", label_size = 14)
# 
# cowplot::plot_grid(top_plots, 
#                    chrom1 + theme(axis.title.x = element_blank()),
#                    region, 
#                    ncol = 1,
#                    labels = c(NA,"C","D"), label_size = 14)
# 
# 
# ggsave(filename = "Plots/SVG_PLOTS/Figure5.pdf", height = 10, width = 6.5, units = "in")
# ggsave(filename = "Plots/SVG_PLOTS/Figure5.pdf", height = 10, width = 6.5, units = "in")
# ggsave(filename = "Plots/SVG_PLOTS/Figure5.pdf", height = 10, width = 6.5, units = "in")


fig5 <- cowplot::plot_grid(n_t ,
                   p_t,
                   region + scale_y_continuous(limits = c(0,160),breaks = c(0,50,100, 150), expand = c(.01, 0)),
                   glct6_r +
                     scale_y_continuous(limits = c(0,160),breaks = c(0,50,100, 150), expand = c(.01, 0)) +
                     theme(axis.title.y = element_blank(),
                           axis.text.y = element_blank(),
                           plot.margin = margin(.1, .1, .1, .1, "cm")) ,
                   ncol = 2,
                   labels = "AUTO", label_size = 14, align = "hv")


ggsave(plot = fig5, filename = "Plots/SVG_PLOTS/Figure5.pdf", height = 6, width = 6.5, units = "in")
ggsave(plot = fig5,filename = "Plots/SVG_PLOTS/Figure5.png", height = 6, width = 6.5, units = "in", dpi = 300)
ggsave(plot = fig5,filename = "Plots/SVG_PLOTS/Figure5.svg", height = 6, width = 6.5, units = "in")


test_dnds <- PopGenome::readData("Data/phylo/Elegans_only/fasta", format = "phylip")


cowplot::plot_grid(n_t ,
                   p_t,
                   region + scale_y_continuous(limits = c(0,20),breaks = c(0,5,10, 20), expand = c(.01, 0)),
                   glct6_r +
                     scale_y_continuous(limits = c(0,20),breaks = c(0,5,10, 20), expand = c(.01, 0)) +
                     theme(axis.title.y = element_blank(),
                           axis.text.y = element_blank(),
                           plot.margin = margin(.1, .1, .1, .1, "cm")) ,
                   ncol = 2,
                   labels = "AUTO", label_size = 14, align = "hv")

  # compare trees
p_tree <- treeio::read.raxml("Data/phylo/Protein/RAxML_bipartitionsBranchLabels.elegans_partition_paralogs")

n_tree <- treeio::read.raxml("Data/phylo/Elegans_only/RAxML_bipartitionsBranchLabels.elegans_partition")

ape::comparePhylo(n_tree@phylo, p_tree@phylo, plot = TRUE, force.rooted = F,use.edge.length = TRUE)
