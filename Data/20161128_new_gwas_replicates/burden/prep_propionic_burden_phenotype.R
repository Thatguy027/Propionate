library(ggplot2)
library(dplyr)
library(tidyr)
library(broom)
library(cegwas)



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
  dplyr::arrange(desc(m.p))%>%
  dplyr::distinct(strain,m.p)

pr_resid$strain <- factor(pr_resid$strain,levels = pr_resid$strain,labels = pr_resid$strain, ordered = T)

ggplot(pr_resid) + 
  aes(y = m.p, x =strain , fill = strain)+
  geom_boxplot()

# p_reg <- pr_resid%>%
#   dplyr::summarise(phenotype = mean(resid, na.rm = T))

pr_pheno <- process_pheno(pr_resid)
maps <- gwas_mappings(pr_pheno,mapping_snp_set = F)

pr_maps <- process_mappings(maps,pr_pheno)

l1_manplot <- pr_maps%>%
  dplyr::filter(CHROM!="MtDNA")%>%
  manplot()
l1_manplot[[1]] + labs(title = "L1 Survival")

ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/burden/singlemarker_propionic_l1.pdf",
       width = 12,
       height =4)

pr_maps%>%
  na.omit()%>%
  dplyr::mutate(col_st = ifelse(strain == "BRC20067", "hotpink3",
                                ifelse(strain == "DL238", "cadetblue3","red")))%>%
  ggplot(.)+
  aes(x=factor(as.character(allele), labels = c("REF","ALT")), y = value)+
  geom_boxplot(fill = "gray90",outlier.colour = NA)+
  geom_jitter(width = .25, aes(fill = col_st, size = col_st,alpha = col_st), shape = 21, color = "black")+
  scale_fill_manual(values = c("hotpink3","cadetblue3","gray50"))+
  scale_size_manual(values = c(6,6,2))+
  scale_alpha_manual(values = c(1,1,.5))+
  facet_grid(~marker)+
  theme_bw()+
  theme(legend.position = 'none')+
  labs(x="",y = "L1 Survival")+
  theme(axis.text.x = ggplot2::element_text(size = 14),
        axis.text.y = ggplot2::element_text(size = 14),
        legend.position = "none",
        axis.title.x = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3),
        axis.title.y = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3),
        strip.text.x = element_text(size = 14, face = "bold"))

ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/burden/singlemarker_propionic_l1_pxg.pdf", 
       height = 6, 
       width = 6)

# extend CI in light of the NIL results
pr_maps$endPOS <- gsub(4129900, 4929900, pr_maps$endPOS)

genes<- process_correlations(variant_correlation(pr_maps, variant_severity = "ALL", condition_trait = F))

save(genes,file = "~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/burden/propionic_l1_arrest_genes.Rda")

genes%>%
  dplyr::filter(strain=="DL238" | strain == "BRC20067")%>%
  dplyr::group_by(CHROM, POS, gene_id, aa_change, nt_change)%>%
  dplyr::mutate(variant_count = ifelse(GT == "ALT", 1, 0))%>%
  dplyr::mutate(sum_alt = sum(variant_count))%>%
  dplyr::ungroup()%>%
  dplyr::filter(strain == "DL238")%>%
  ggplot(.)+
  aes(x=POS/1e6, y = -log10(corrected_spearman_cor_p))+
  geom_point(aes(fill = factor(variant_count)),alpha = 0.7, shape = 21, color = "black", size =2)+
  scale_fill_manual(values=c("cadetblue3","hotpink3"))+
  theme_bw()+
  theme(legend.position = 'none')+
  labs(x="Genomic Position (Mb)",y = expression(-log[10](italic(p))))+
  theme(axis.text.x = ggplot2::element_text(size = 14),
        axis.text.y = ggplot2::element_text(size = 14),
        legend.position = "none",
        axis.title.x = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3),
        axis.title.y = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3),
        strip.text.x = element_text(size = 14, face = "bold"))

ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/burden/singlemarker_propionic_l1_finemapping.pdf", 
       height = 6, 
       width = 10)

genes%>%
  dplyr::filter(strain=="DL238" | strain == "BRC20067")%>%
  dplyr::group_by(CHROM, POS, gene_id, aa_change, nt_change)%>%
  dplyr::mutate(variant_count = ifelse(GT == "ALT", 1, 0))%>%
  dplyr::mutate(sum_alt = sum(variant_count),
                log10p = -log10(corrected_spearman_cor_p))%>%
  dplyr::ungroup()%>%
  dplyr::filter(strain == "DL238")%>%
  dplyr::select(CHROM,POS,aa_change, gene_id, gene_name,log10p)%>%
  View()

# prepare phenotype ped file
prop_burden <- pr_resid%>%
  # dplyr::mutate(trait = "prop_L1_arrest")%>%
  dplyr::mutate(Fam = "elegans", Sample = strain, Paternal = 0, Maternal = 0, Sex = 2)%>%
  # tidyr::spread(trait,value)%>%
  dplyr::ungroup()%>%
  dplyr::select(-strain)%>%
  dplyr::select(Fam,Sample,Paternal,Maternal,Sex, phenotype=m.p)

prop_burden[is.na(prop_burden)] <- -9


write.table(prop_burden,"~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/burden/propionic_l1_arrest.ped",
            quote = F, col.names = F, row.names = F)

# docker run -i -t -v $(pwd):/home -w /home zhanxw/rvtests /rvtests/executable/rvtest --pheno propionic_l1_arrest.ped --inVcf snp_indel_snpEff.vcf.gz --freqUpper .05 --freqLower 0.003 --out propionic_burden --geneFile refFlat.ws245.txt --burden cmc --vt price --kernel skat

# process burden mappings
skat <- data.table::fread("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/burden/propionic_burden.Skat.assoc")%>%
  dplyr::rowwise()%>%
  dplyr::mutate(CHROM = strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][1],
                POS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][1]),
                endPOS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][2]))%>%
  dplyr::mutate(size = abs(POS-endPOS))

skat_pr  <- dplyr::filter(skat,CHROM!="MtDNA")%>%
  dplyr::ungroup()%>%
  # dplyr::filter(NumVar > 1,size >500)%>%
  dplyr::filter(NumVar > 1)%>%
  dplyr::mutate(significant = ifelse(Pvalue < .05/n(), TRUE,FALSE ))%>%
  dplyr::arrange(Pvalue)

write.table(skat_pr, 
            file ="~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/burden/processed_SKAT.tsv", 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

skat_pr%>%
  # dplyr::filter(significant == TRUE)%>%
  ggplot()+
  aes(x = POS/1e6, y = -log10(Pvalue), alpha = 0.5, color = significant)+
  geom_point()+
  scale_color_manual(values=c("black","red"))+
  facet_grid(.~CHROM, scales = "free", space = "free")+
  theme_bw()+
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16),
                 axis.text.y = ggplot2::element_text(size = 16),
                 axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                 axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.x = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 plot.title = ggplot2::element_text(size = 24, face = "bold", vjust = 1), 
                 panel.background = ggplot2::element_rect(color = "black",size = 1.2),
                 legend.position = "none")+
  labs(x = "Genomic Position (Mb)", y = expression(-log[10](italic(p))))

ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/burden/skat_propionic_l1.pdf",
       width = 12,
       height =4)

vtprice <- data.table::fread("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/burden/propionic_burden.VariableThresholdPrice.assoc")%>%
  dplyr::rowwise()%>%
  dplyr::mutate(CHROM = strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][1],
                POS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][1]),
                endPOS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][2]))%>%
  dplyr::mutate(size = abs(POS-endPOS))

vtprice_pr<- vtprice%>%
  ungroup()%>%
  dplyr::filter(CHROM!="MtDNA",NumVar>1)%>%
  dplyr::mutate(significant = ifelse(PermPvalue < .05/n(), TRUE,FALSE ))%>%
  dplyr::arrange(desc(Stat))
  
write.table(vtprice_pr, 
            file ="~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/burden/processed_VTprice.tsv", 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

# %>%
#   dplyr::filter(CHROM!="MtDNA",NumVar>1,size >500)

vtprice_pr%>%
  ggplot()+
  aes(x = POS/1e6, y = Stat, size = NumVar, alpha = 0.5, color = significant)+
  geom_point()+
  scale_color_manual(values=c("black","red"))+
  facet_grid(.~CHROM, scales = "free", space = "free")+
  theme_bw()+
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16),
                 axis.text.y = ggplot2::element_text(size = 16),
                 axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                 axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.x = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 plot.title = ggplot2::element_text(size = 24, face = "bold", vjust = 1), 
                 panel.background = ggplot2::element_rect(color = "black",size = 1.2),
                 legend.position = "none")+
  labs(x = "Genomic Position (Mb)", y = "Test Statisitic")



ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/burden/VTprice_propionic_l1.pdf",
       width = 12,
       height =4)

# process burden mappings
cmc <- data.table::fread("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/burden/propionic_burden.CMC.assoc")%>%
  dplyr::rowwise()%>%
  dplyr::mutate(CHROM = strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][1],
                POS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][1]),
                endPOS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][2]))%>%
  dplyr::mutate(size = abs(POS-endPOS))

cmc_pr  <- dplyr::filter(cmc,CHROM!="MtDNA")%>%
  dplyr::ungroup()%>%
  dplyr::filter(NumVar > 1,size >500)%>%
  dplyr::mutate(significant = ifelse(Pvalue < .05/n(), TRUE,FALSE ))

cmc_pr%>%
  # dplyr::filter(significant == TRUE)%>%
  ggplot()+
  aes(x = POS/1e6, y = -log10(Pvalue), alpha = 0.5, color = significant)+
  geom_point()+
  scale_color_manual(values=c("black","red"))+
  facet_grid(.~CHROM, scales = "free", space = "free")+
  theme_bw()+
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16),
                 axis.text.y = ggplot2::element_text(size = 16),
                 axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                 axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.x = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 plot.title = ggplot2::element_text(size = 24, face = "bold", vjust = 1), 
                 panel.background = ggplot2::element_rect(color = "black",size = 1.2),
                 legend.position = "none")+
  labs(x = "Genomic Position (Mb)", y = expression(-log[10](italic(p))))

ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/burden/cmc_propionic_l1.pdf",
       width = 12,
       height =4)


glct3_stop <- c("MY23","JU2879","ED3049","CX11276","ECA363","QX1793","QX1791","ECA189","ED3046","ECA396","ECA191","DL238","ECA36","QX1792","NIC252","XZ1513","ECA372")

pheno_by_glct3 <- pr_resid%>%
  dplyr::mutate(glct3stop = ifelse(strain %in% glct3_stop, "STOP","WT"))%>%
  ggplot()+
  aes(x = glct3stop, y = m.p)+
  geom_jitter(width=0.25)+
  geom_boxplot(alpha = 0.5, outlier.colour = NA)+
  theme_bw()+
  labs(x = "glct-3 stop", y = "L1 survival")+
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16),
                 axis.text.y = ggplot2::element_text(size = 16),
                 axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                 axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.x = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 plot.title = ggplot2::element_text(size = 24, face = "bold", vjust = 1), 
                 panel.background = ggplot2::element_rect(color = "black",size = 1.2),
                 legend.position = "none")
pheno_by_glct3
ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20161128_new_gwas_replicates/burden/glct3_stop_pxg.pdf",
       width = 6,
       height =4)
