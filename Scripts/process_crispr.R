library(tidyverse)

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
setwd("..")

phen <- data.table::fread("Data/glct_causality/20190624_glct_del_swap.csv")

phen %>%
  ggplot()+
  aes(x = strain, y = survival, fill = strain)+
  geom_boxplot(alpha = 0.8)+
  scale_fill_manual(values=c("hotpink3","gray60","gray60","cadetblue3"))+
  scale_x_discrete(labels=c("BRC20067" = "BRC20067", 
                            "DL238" = "DL238", 
                            "BRC20067_del" = expression("BRC20067 "~ Delta~italic("glct-3")),
                            "BRC20067_swap" = "BRC20067 (Gly16*)"))+
  ggbeeswarm::geom_beeswarm(cex = 0.7, priority = "density")+
  theme_classic(20)+
  labs(y = "L1 Survival (%)") +
  # ylim(0,100)+
  theme(axis.title.x = element_blank(), legend.position = "none")

ggsave("Plots/Crispr_pheno.pdf", height = 4, width = 10)
ggsave("Plots/Crispr_pheno.png", height = 4, width = 10, dpi = 300)
