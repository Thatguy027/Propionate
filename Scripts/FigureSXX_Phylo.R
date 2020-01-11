library(tidyverse)
library(ggtree)

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
setwd("..")

source("Scripts/Figure_Themes.R")

# load protein tree
p_tree <- treeio::read.raxml("Data/phylo/Protein/RAxML_bipartitionsBranchLabels.elegans_partition")

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
  # xlim(0,2.5)+ 
  # geom_label2(aes(label=bootstrap), size = 2) +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) 
p_t

# load protein tree
d_tree <- treeio::read.raxml("Data/phylo/TREE/RAxML_bipartitionsBranchLabels.partition")

# original <- p_tree@phylo$tip.label
# clean_names <- c("GLCT-2", "GLCT-3", "GLCT-6", "GLCT-4","GLCT-5","GLCT-1")
# rename_df <- data.frame(label = original, newlab = clean_names)
# 
# p_tree_clean <- treeio::rename_taxa(p_tree, rename_df, label, newlab)

d_t <- ggtree(d_tree) + 
  scale_color_viridis_c(direction = -1) +
  theme(legend.position="right")+
  theme_classic(12) +
  geom_tiplab(size = 2) +
  # xlim(0,2.5)+ 
  # geom_label2(aes(label=bootstrap), size = 2) +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) 
d_t
