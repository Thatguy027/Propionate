library(data.table)
library(dplyr)
library(tidyr)
library(cegwas)
library(ggplot2)
library(broom)
library(cowplot)

p <- data.table::fread("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170103_set1-3.csv")%>%
  dplyr::select(strain, ratio = average, replicate, set)

# pull out strains that were repeated multiple days, take median of replicates (mean & median produces same QTL)
rep_strains <- p %>%
  dplyr::filter(strain %in% c("CX11307", "DL238", "EG4725", "JU775", "JU830", "MY23"))%>%
  dplyr::group_by(strain,set,replicate)%>%
  dplyr::mutate(m_rat = median(ratio, na.rm = T))%>%
  dplyr::select( replicate, set, strain, m_rat)%>%
  dplyr::ungroup()%>%
  dplyr::group_by(strain)

# fit linear model to see how set number affects survival
set_fit <- lm(m_rat ~ set, data = rep_strains)
coeffs <- coef(set_fit)
# visualize model
plot(rep_strains$set, rep_strains$m_rat)
abline(8.663427, 3.579320, col = "red")

# correct phenotypes of all GWAS strains based on the parameters identified above
pr_p <- p %>%
  dplyr::group_by(strain,set,replicate)%>%
  # dplyr::mutate(m_rat = mean(ratio, na.rm = T))%>%
  dplyr::select( replicate, set, strain, ratio)%>%
  dplyr::ungroup()%>%
  dplyr::group_by(strain)%>%
  dplyr::filter(strain != "")

# adjust
adjusted_survival <- pr_p %>%
  dplyr::mutate(set_adj = coeffs[1] + coeffs[2]*ratio)

# visualize to make sure that the uncorrected vs. corrected values correspond to the fit line identified in the model above. 
plot(adjusted_survival$ratio, adjusted_survival$set_adj)
abline(8.663427, 3.579320, col = "red")

# fit linear model to see how set number affects survival
replicate_fit <- lm(set_adj ~ replicate, data = adjusted_survival)
coeffs <- coef(replicate_fit)
# visualize model
plot(adjusted_survival$replicate, adjusted_survival$ratio)
abline(84.36347, -10.66466, col = "red")

# correct phenotypes of all GWAS strains based on the parameters identified above
pr_p <- adjusted_survival %>%
  dplyr::select( replicate, set, strain, set_adj)%>%
  dplyr::ungroup()%>%
  dplyr::group_by(strain)%>%
  dplyr::filter(strain != "")

# adjust
adjusted_survival <- pr_p %>%
  dplyr::mutate(replicate_adj = coeffs[1] + coeffs[2]*set_adj)

# visualize to make sure that the uncorrected vs. corrected values correspond to the fit line identified in the model above. 
plot(adjusted_survival$set_adj, adjusted_survival$replicate_adj)
abline(84.36347, -10.66466, col = "red")

adjusted_survival_pr <-  adjusted_survival%>%
  group_by(strain)%>%
  dplyr::summarise(phenotype = mean(replicate_adj, na.rm = T))
  
pr_pheno <- process_pheno(adjusted_survival_pr)
maps <- gwas_mappings(pr_pheno)

pr_maps <- process_mappings(maps,pr_pheno, BF = 4)

save(pr_maps, file ="~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170402_most_updated/set-rep_regression-maps.Rda")



manplot(pr_maps)
ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170402_regression_set-replicate/manplot.pdf",
       width = 10,
       height = 5)

pxg_plot(pr_maps,color_strains = c("CX11307", "DL238", "EG4725", "JU775"))

ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170402_regression_set-replicate/pxg.pdf",
       width = 6,
       height = 5)


