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
  dplyr::mutate(m_rat = median(ratio, na.rm = T))%>%
  dplyr::select( replicate, set, strain, m_rat)%>%
  dplyr::ungroup()%>%
  dplyr::group_by(strain)

# adjust
adjusted_survival <- pr_p %>%
  dplyr::mutate(set_adj = coeffs[1] + coeffs[2]*m_rat  )

# visualize to make sure that the uncorrected vs. corrected values correspond to the fit line identified in the model above. 
plot(adjusted_survival$m_rat, adjusted_survival$set_adj)
abline(8.663427, 3.579320, col = "red")


adjusted_survival_pr <-  adjusted_survival%>%
  dplyr::select(strain, replicate, set_adj)%>%
  dplyr::filter(strain != "")%>%
  dplyr::group_by(strain,replicate)%>%
  dplyr::mutate(m_rat = mean(set_adj, na.rm = T))%>%
  dplyr::select( replicate, strain, m_rat)%>%
  dplyr::ungroup()%>%
  dplyr::distinct(strain, m_rat,replicate,.keep_all = T)%>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(mph = median(m_rat, na.rm = T),
                sph = sd(m_rat, na.rm = T))%>%
  dplyr::filter(m_rat <= 1.5*sph+mph,
                m_rat >= mph-1.5*sph)


mean_adjusted <- adjusted_survival_pr%>%
  dplyr::summarise(phenotype = mean(m_rat, na.rm =T))

pr_pheno <- process_pheno(mean_adjusted)
maps <- gwas_mappings(pr_pheno)

pr_maps <- process_mappings(maps,pr_pheno, BF = 4)

save(pr_maps, file ="~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170402_most_updated/correct_with_six.Rda")




