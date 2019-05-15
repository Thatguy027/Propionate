library(data.table)
library(dplyr)
library(tidyr)
library(cegwas)
library(ggplot2)
library(broom)
library(cowplot)

p <- data.table::fread("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170103_set1-3.csv")%>%
  dplyr::select(strain, ratio = average, replicate, set)

# correct phenotypes of all GWAS strains based on the parameters identified above
pr_p <- p %>%
  dplyr::group_by(strain)%>%
  dplyr::mutate(m_rat = median(ratio, na.rm = T))%>%
  dplyr::select( replicate, set, strain, m_rat)%>%
  dplyr::ungroup()%>%
  dplyr::group_by(strain)%>%
  dplyr::summarise(phenotype = mean(m_rat))

pr_pheno <- process_pheno(pr_p)
maps <- gwas_mappings(pr_pheno)

pr_maps <- process_mappings(maps,pr_pheno, BF = 4)


save(pr_maps, file ="~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170402_most_updated/noRegressionMaps.Rda")



manplot(pr_maps)
ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170402_mean_analysis_output/manplot.pdf",
       width = 10,
       height = 5)
pxg_plot(pr_maps,color_strains = c("CX11307", "DL238", "EG4725", "JU775"))

ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170402_mean_analysis_output/pxg.pdf",
       width = 6,
       height = 5)




