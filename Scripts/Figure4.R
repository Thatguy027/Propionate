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
                            "BRC20067_Del" = expression(atop("BRC20067",Delta~italic("glct-3"))),
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
             size = 1) +
  geom_point(data = strains_330 %>% dplyr::filter(!is.na(aa_change)), 
             aes(x=as.numeric(long), y=as.numeric(lat)), 
             fill = "cadetblue3",
             shape = 21, 
             size = 2) +
  scale_fill_manual(values = c("cadetblue3","hotpink3"))+
  theme_map()+
  scale_y_continuous(expand=c(0.01,0.01))+
  scale_x_continuous(expand=c(0.01,0.01))

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



#####################################################################################################################################################################
# load tree
tree <- ape::read.tree(glue::glue("Data/whole_genome_tree/330_genome.raxml.bestTree"))

# highlight branches for strains of interest
branch_strains <- list(CONNECT = alt_strains$strain)

tree_pt_h <- ggtree::groupOTU(tree, branch_strains)

glct_tree <- ggtree(tree_pt_h,
       branch.length="rate", 
       aes(color=group), size = 0.25) + 
  scale_color_manual(values=c("hotpink3", "cadetblue3"), 
                     name = "GLCT-3\nallele", 
                     labels=c("REF", "Gly16*")) + 
  theme(legend.position="right")+
  theme_classic(12) + 
  coord_flip() + 
  scale_x_reverse() +
  scale_y_continuous(expand=c(0.01,0.01))+
  theme(axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        axis.text.y = element_text(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())



# top_panel <- cowplot::plot_grid(lemon::reposition_legend(gwas_split + theme(legend.box.background = element_rect(),
#                                                                             legend.box.margin = margin(.1, .1, .1, .1)),'top right'),
#                                 pheno_plot + theme(axis.title.y = element_blank()), 
#                                 nrow = 1, 
#                                 rel_widths = c(2,3),
#                                 label_size = 14,
#                                 labels = c("A","B"), align = "hv")

top_panel <- cowplot::plot_grid(gwas_split + theme(legend.position = "none"),
                                pheno_plot + theme(axis.title.y = element_blank()), 
                                nrow = 1, 
                                rel_widths = c(2,4),
                                label_size = 14,
                                labels = c("A","B"), align = "hv")

fig4 <- cowplot::plot_grid(top_panel,
                   map, 
                   lemon::reposition_legend(glct_tree, "bottom left"), 
                   nrow = 3, 
                   label_size = 14,
                   labels = c(NA,"C","D"), align = "h", scale = c(1,1,0.95))



ggsave(plot = fig4, "Plots/SVG_PLOTS/Figure4.pdf", height = 10, width = 6.5, units = "in")
ggsave(plot = fig4, "Plots/SVG_PLOTS/Figure4.png", height = 10, width = 6.5, units = "in", dpi = 300)
ggsave(plot = fig4, "Plots/SVG_PLOTS/Figure4.svg", height = 10, width = 6.5, units = "in")
