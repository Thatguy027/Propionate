kinship_correction <- function(y, kin = kinship) {
  
  K <- kin[colnames(kin) %in% y$strain, rownames(kin) %in% y$strain]
  y <- y %>% dplyr::filter(strain %in% colnames(K))
  
  model = regress::regress(as.vector(y$value)~1,~K, pos = c(TRUE, TRUE))	
  
  summary(model)
  
  #This err.cov is the same as err.cov in Dan's code using estVC
  err.cov = regress::summary.regress(model)$sigma[1]*K+regress::summary.regress(model)$sigma[2]*diag(nrow(K))
  
  eW = eigen(err.cov, symmetric = TRUE)
  
  if(min(eW$values) < 0 && abs(min(eW$values)) > sqrt(.Machine$double.eps)){
  }else{
    eW$values[eW$values <= 0] = Inf
  } 
  
  err.cov = eW$vector %*% diag(eW$values^-0.5) %*% t(eW$vector)
  new.pheno <- err.cov %*% y$value
  
  new.pheno <- data.frame(strain = y$strain, corrected_pheno = new.pheno)
  
  plot(new.pheno$corrected_pheno, y$value)
  
  return(new.pheno)
}

kin <- cegwas::kinship

# # # # # 


load("/Users/stefan/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170403_strain_fixed_mappings_withAB1.Rda")

mixed_mappings <- pr_maps_strain_fe
manplot(mixed_mappings)
ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/figures/20170403_strain_fixed_mappings_withAB1.pdf",
       width=11,
       height=5)
pxg_plot(mixed_mappings)

original <- na.omit(mixed_mappings)%>%
  dplyr::select(strain, trait, value)

out_strain_random <- kinship_correction(original,kin)

# # # # # # 

load("/Users/stefan/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170403_strain_fixed_mappings_noAB1.Rda")

noAB1_mixed_mappings <- pr_maps_strain_fe
manplot(noAB1_mixed_mappings)
ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/figures/20170403_strain_fixed_mappings_noAB1.pdf",
       width=11,
       height=5)
pxg_plot(noAB1_mixed_mappings)

original <- na.omit(noAB1_mixed_mappings)%>%
  dplyr::select(strain, trait, value)%>%
  dplyr::distinct(strain, trait, value)

out_strain_random <- kinship_correction(original,kin)

# # # # # # 

load("/Users/stefan/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/phenotype_correct_with_sixCoeffs_maps.Rda")

six_strain_correction_mappings <- pr_maps
manplot(six_strain_correction_mappings)
ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/figures/phenotype_correct_with_sixCoeffs.pdf",
       width=11,
       height=5)
pxg_plot(six_strain_correction_mappings)

original <- na.omit(six_strain_correction_mappings)%>%
  dplyr::select(strain, trait, value)%>%
  dplyr::distinct(strain, trait, value)

out_strain_random <- kinship_correction(original,kin)

# # # # # ## 

load("/Users/stefan/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170403_simple_setRegression_maps.Rda")

set_regression_mappings <- pr_maps
manplot(set_regression_mappings)
ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/figures/20170403_simple_setRegression_maps.pdf",
       width=11,
       height=5)
pxg_plot(set_regression_mappings)

original <- na.omit(set_regression_mappings)%>%
  dplyr::select(strain, trait, value)%>%
  dplyr::distinct(strain, trait, value)

out_strain_random <- kinship_correction(original,kin)

# # # # # ## 

load("/Users/stefan/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170403_simple_set-and-ratioRegression_maps.Rda")

setReplicate_regression_mappings <- pr_maps
manplot(setReplicate_regression_mappings)
ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/figures/20170403_simple_setReplicateRegression_maps.pdf",
       width=11,
       height=5)
pxg_plot(setReplicate_regression_mappings)

original <- na.omit(setReplicate_regression_mappings)%>%
  dplyr::select(strain, trait, value)%>%
  dplyr::distinct(strain, trait, value)

out_strain_random <- kinship_correction(original,kin)

# # # # # # 


load("/Users/stefan/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170403_mean_maps.Rda")

mean_mappings <- pr_maps
manplot(mean_mappings)
ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/figures/20170403_mean_maps.pdf",
       width=11,
       height=5)
pxg_plot(mean_mappings)

original <- na.omit(mean_mappings)%>%
  dplyr::select(strain, trait, value)%>%
  dplyr::distinct(strain, trait, value)

out_strain_random <- kinship_correction(original,kin)

# # ##### ### 

load("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170403_arcsin_mixed_strainFE_maps.Rda")

arcsin_maps_strain_fe <- pr_maps_strain_fe
manplot(arcsin_maps_strain_fe)
ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/figures/20170403_mean_maps.pdf",
       width=11,
       height=5)
pxg_plot(arcsin_maps_strain_fe)

original <- na.omit(arcsin_maps_strain_fe)%>%
  dplyr::select(strain, trait, value)%>%
  dplyr::distinct(strain, trait, value)

out_strain_random <- kinship_correction(original,kin)

# # # # # # ## 

load("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170403_arcsin_mixed_strainFE_maps.Rda")

arcsin_maps_strain_fe <- pr_maps_strain_fe
manplot(arcsin_maps_strain_fe)
ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/figures/20170403_mean_maps.pdf",
       width=11,
       height=5)
pxg_plot(arcsin_maps_strain_fe)

original <- na.omit(arcsin_maps_strain_fe)%>%
  dplyr::select(strain, trait, value)%>%
  dplyr::distinct(strain, trait, value)

out_strain_random <- kinship_correction(original,kin)

# # # # # # ## 

load("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170403_arcsin_SetRegress_maps.Rda")

arcsin_setRegress <- pr_maps
manplot(arcsin_setRegress)
ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/figures/20170403_arcsin_SetRegress_maps.pdf",
       width=11,
       height=5)
pxg_plot(arcsin_setRegress)

original <- na.omit(arcsin_setRegress)%>%
  dplyr::select(strain, trait, value)%>%
  dplyr::distinct(strain, trait, value)

out_strain_random <- kinship_correction(original,kin)

# # # # # # ## 

load("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/20170403_arcsin_SetRepRegress_maps.Rda")

arcsin_setRepRegress <- pr_maps
manplot(arcsin_setRepRegress)
ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/figures/20170403_arcsin_SetRepRegress_maps.pdf",
       width=11,
       height=5)
pxg_plot(arcsin_setRepRegress)

original <- na.omit(arcsin_setRepRegress)%>%
  dplyr::select(strain, trait, value)%>%
  dplyr::distinct(strain, trait, value)

out_strain_random <- kinship_correction(original,kin)

