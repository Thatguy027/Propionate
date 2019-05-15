
library(data.table)
library(dplyr)
library(tidyr)
library(cegwas)
library(ggplot2)
library(broom)
library(cowplot)
library(lme4)

p <- data.table::fread("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170103_set1-3.csv")%>%
  dplyr::select(strain, ratio = average, replicate, set)%>%
  dplyr::filter(strain != "", strain != "AB1")

pfit <- lmer(ratio ~ (1|replicate) + (1|set) + strain, data = p)

mixed_p <- data.frame(strain = p$strain, 
                      replicate= p$replicate,
                      set= p$set,
                      ratio = residuals(pfit))%>%
  group_by(strain)%>%
  summarise(ratio = mean(ratio))%>%
  arrange(ratio)

pr_pheno <- process_pheno(mixed_p)
maps_strain_fe <- gwas_mappings(pr_pheno)
pr_maps_strain_fe <- process_mappings(maps_strain_fe,pr_pheno, BF = 4, snp_grouping = 100)

save(pr_maps_strain_fe, file = "~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170403_strain_fixed_mappings_noAB1.Rda")

# # # # 

p <- data.table::fread("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170103_set1-3.csv")%>%
  dplyr::select(strain, ratio = average, replicate, set)%>%
  dplyr::filter(strain != "")

pfit <- lmer(ratio ~ (1|replicate) + (1|set) + strain, data = p)

mixed_p <- data.frame(strain = p$strain, 
                      replicate= p$replicate,
                      set= p$set,
                      ratio = residuals(pfit))%>%
  group_by(strain)%>%
  summarise(ratio = mean(ratio))%>%
  arrange(ratio)

pr_pheno <- process_pheno(mixed_p)
maps_strain_fe <- gwas_mappings(pr_pheno)
pr_maps_strain_fe <- process_mappings(maps_strain_fe,pr_pheno, BF = 4, snp_grouping = 100)

save(pr_maps_strain_fe, file = "~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170403_strain_fixed_mappings_withAB1.Rda")

# # # # 

p <- data.table::fread("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170103_set1-3.csv")%>%
  dplyr::select(strain, ratio = average, replicate, set)%>%
  dplyr::filter(strain != "", strain != "AB1")

pfit <- lmer(ratio ~ replicate + set + (1|strain), data = p)

mixed_p <- data.frame(strain = p$strain, 
                      replicate= p$replicate,
                      set= p$set,
                      ratio = residuals(pfit))%>%
  group_by(strain)%>%
  summarise(ratio = mean(ratio))%>%
  arrange(ratio)

pr_pheno <- process_pheno(mixed_p)
maps_strain_re <- gwas_mappings(pr_pheno)
pr_maps_strain_re <- process_mappings(maps_strain_re,pr_pheno, BF = 4, snp_grouping = 100)

save(pr_maps_strain_fe, file = "~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170403_strain_random_mappings_noAB1.Rda")



# # # # strain as a random variable, replicate and set as fixed 

p <- data.table::fread("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170103_set1-3.csv")%>%
  dplyr::select(strain, ratio = average, replicate, set)%>%
  dplyr::filter(strain != "")

pfit <- lmer(ratio ~ replicate + set + (1|strain), data = p)

mixed_p <- data.frame(strain = p$strain, 
                      replicate= p$replicate,
                      set= p$set,
                      ratio = residuals(pfit))%>%
  group_by(strain)%>%
  summarise(ratio = mean(ratio))%>%
  arrange(ratio)

pr_pheno <- process_pheno(mixed_p)
maps_strain_re <- gwas_mappings(pr_pheno)
pr_maps_strain_re <- process_mappings(maps_strain_re,pr_pheno, BF = 4, snp_grouping = 100)

save(pr_maps_strain_fe, file = "~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170403_strain_random_mappings_withAB1.Rda")

# # # # # # find how set affects replicated strains, adjust other strains by these coefficients, and map

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
  dplyr::filter(m_rat <= 2*sph+mph,
                m_rat >= mph-2*sph)


mean_adjusted <- adjusted_survival_pr%>%
  dplyr::summarise(phenotype = mean(m_rat, na.rm =T))

pr_pheno <- process_pheno(mean_adjusted)
maps <- gwas_mappings(pr_pheno)

pr_maps <- process_mappings(maps,pr_pheno, BF = 4)

save(pr_maps, file ="~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/phenotype_correct_with_sixCoeffs_maps.Rda")

# # # # # # original way

p <- data.table::fread("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170103_set1-3.csv")%>%
  dplyr::select(strain, ratio = average, replicate, set)

pr_resid <- p%>%
  ungroup()%>%
  dplyr::do(augment(lm(ratio ~ set, data =.)))%>%
  dplyr::rename(resid = .resid)%>%
  data.frame()%>%
  dplyr::left_join(p, ., by = c("ratio","set"))%>%
  dplyr::distinct(ratio,set,replicate,strain,.keep_all = T) %>%
  dplyr::ungroup()%>%
  dplyr::select(strain, resid)%>%
  dplyr::group_by(strain)%>%
  dplyr::mutate(mph = median(resid, na.rm = T),
                sph = sd(resid, na.rm = T))%>%
  dplyr::filter(resid <= 2*sph+mph,
                resid >= mph-2*sph)%>%
  dplyr::mutate(m.p = mean(resid))%>%
  dplyr::arrange(desc(m.p))

pr_resid$strain <- factor(pr_resid$strain,levels = pr_resid$strain,labels = pr_resid$strain, ordered = T)

ggplot(pr_resid) + 
  aes(y = resid, x =strain , fill = strain)+
  geom_boxplot()

p_reg <- pr_resid%>%
  dplyr::summarise(phenotype = mean(resid, na.rm = T))

pr_pheno <- process_pheno(p_reg)
maps <- gwas_mappings(pr_pheno)

pr_maps <- process_mappings(maps,pr_pheno)

save(pr_maps, file ="~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170403_simple_setRegression_maps.Rda")

# # # # 

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
  dplyr::arrange(desc(m.p))

pr_resid$strain <- factor(pr_resid$strain,levels = pr_resid$strain,labels = pr_resid$strain, ordered = T)

ggplot(pr_resid) + 
  aes(y = resid, x =strain , fill = strain)+
  geom_boxplot()

p_reg <- pr_resid%>%
  dplyr::summarise(phenotype = mean(resid, na.rm = T))

pr_pheno <- process_pheno(p_reg)
maps <- gwas_mappings(pr_pheno)

pr_maps <- process_mappings(maps,pr_pheno)

save(pr_maps, file ="~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170403_simple_set-and-ratioRegression_maps.Rda")

# # # # take means and map

p <- data.table::fread("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170103_set1-3.csv")%>%
  dplyr::select(strain, ratio = average)

pr_resid <- p%>%
  group_by(strain)%>%
  dplyr::summarise(phenotype = mean(ratio, na.rm = T))%>%
  dplyr::filter(strain != "")

pr_pheno <- process_pheno(pr_resid)
maps <- gwas_mappings(pr_pheno)

pr_maps <- process_mappings(maps,pr_pheno, BF = 4)

save(pr_maps, file ="~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170403_mean_maps.Rda")

# # # # 


p <- data.table::fread("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170103_set1-3.csv")%>%
  dplyr::select(strain, ratio = average, replicate, set)%>%
  dplyr::mutate(fraction = ratio/100)%>%
  dplyr::mutate(asin_ratio = asin(sqrt(fraction)))%>%
  dplyr::select(strain, ratio = asin_ratio, replicate, set)

pr_resid <- p%>%
  ungroup()%>%
  dplyr::do(augment(lm(ratio ~ set, data =.)))%>%
  dplyr::rename(resid = .resid)%>%
  data.frame()%>%
  dplyr::left_join(p, ., by = c("ratio","set"))%>%
  dplyr::distinct(ratio,set,replicate,strain,.keep_all = T) %>%
  dplyr::ungroup()%>%
  dplyr::select(strain, resid)%>%
  dplyr::group_by(strain)%>%
  dplyr::mutate(mph = median(resid, na.rm = T),
                sph = sd(resid, na.rm = T))%>%
  dplyr::filter(resid <= 2*sph+mph,
                resid >= mph-2*sph)%>%
  dplyr::mutate(m.p = mean(resid))%>%
  dplyr::arrange(desc(m.p))

p_reg <- pr_resid%>%
  dplyr::summarise(phenotype = mean(resid, na.rm = T))

pr_pheno <- process_pheno(p_reg)
maps <- gwas_mappings(pr_pheno)

pr_maps <- process_mappings(maps,pr_pheno, BF = 4)

save(pr_maps, file ="~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170403_arcsin_SetRegress_maps.Rda")

# # # # # #

p <- data.table::fread("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170103_set1-3.csv")%>%
  dplyr::select(strain, ratio = average, replicate, set)%>%
  dplyr::mutate(fraction = ratio/100)%>%
  dplyr::mutate(asin_ratio = asin(sqrt(fraction)))%>%
  dplyr::select(strain, ratio = asin_ratio, replicate, set)%>%
  dplyr::filter(strain != "")

pfit <- lmer(ratio ~ (1|replicate) + (1|set) + strain, data = p)

mixed_p <- data.frame(strain = p$strain, 
                      replicate= p$replicate,
                      set= p$set,
                      ratio = residuals(pfit))%>%
  group_by(strain)%>%
  summarise(ratio = mean(ratio))%>%
  arrange(ratio)

pr_pheno <- process_pheno(mixed_p)
maps_strain_fe <- gwas_mappings(pr_pheno)
pr_maps_strain_fe <- process_mappings(maps_strain_fe,pr_pheno, BF = 3, snp_grouping = 100)

save(pr_maps_strain_fe, file ="~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170403_arcsin_mixed_strainFE_maps.Rda")

# # # # # # 


p <- data.table::fread("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170103_set1-3.csv")%>%
  dplyr::select(strain, ratio = average, replicate, set)%>%
  dplyr::mutate(fraction = ratio/100)%>%
  dplyr::mutate(asin_ratio = asin(sqrt(fraction)))%>%
  dplyr::select(strain, ratio = asin_ratio, replicate, set)%>%
  dplyr::filter(strain != "")

pr_resid <- p%>%
  ungroup()%>%
  dplyr::do(augment(lm(ratio ~ set+ replicate, data =.)))%>%
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
  dplyr::arrange(desc(m.p))

p_reg <- pr_resid%>%
  dplyr::summarise(phenotype = mean(resid, na.rm = T))

pr_pheno <- process_pheno(p_reg)
maps <- gwas_mappings(pr_pheno)

pr_maps <- process_mappings(maps,pr_pheno, BF = 4)

save(pr_maps, file ="~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170403_arcsin_SetRepRegress_maps.Rda")
