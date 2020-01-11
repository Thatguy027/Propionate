
glct3_gene <- cegwas2::query_vcf("WBGene00011781", vcf = "~/UCLA/Genomics_Data/VCFs/328_Soft_csq_annotated.vcf.gz")


glct_gene_matrix <- glct3_gene %>%
  dplyr::filter(FILTER=="PASS") %>%
  dplyr::select(POS, SAMPLE, ALT, a1,a2) %>%
  dplyr::mutate(nGT = ifelse(a1 == ALT | a2 == ALT, 1, 0)) %>%
  dplyr::select(POS, SAMPLE, nGT) %>%
  dplyr::distinct(POS, SAMPLE, .keep_all=T) %>%
  tidyr::spread(POS, nGT)

glct_gene_matrix$altct <- rowSums(glct_gene_matrix[,2:ncol(glct_gene_matrix)], na.rm = T)

glct_gene_matrix <- glct_gene_matrix %>% dplyr::filter(altct > 0 | SAMPLE %in% c("BRC20067", "N2"))

row.names(glct_gene_matrix) <- glct_gene_matrix$SAMPLE

glct_gene_tree = ape::nj(ape::dist.gene(glct_gene_matrix, pairwise.deletion = T,method = "percentage"))


plot(glct_gene_tree)
ggtree(glct_gene_tree) +
  geom_tiplab(size = 2) +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) 

min(finemap$POS)
max(finemap$POS)


region <- cegwas2::query_vcf("I:12000000-12500000", vcf = "~/UCLA/Genomics_Data/VCFs/328_Soft_csq_annotated.vcf.gz")


glct_region <- region %>%
  dplyr::filter(FILTER=="PASS") %>%
  dplyr::select(POS, SAMPLE, ALT, a1,a2) %>%
  dplyr::mutate(nGT = ifelse(a1 == ALT | a2 == ALT, 1, 0)) %>%
  dplyr::select(POS, SAMPLE, nGT) %>%
  dplyr::distinct(POS, SAMPLE, .keep_all=T) %>%
  tidyr::spread(POS, nGT)

row.names(glct_region) <- glct_region$SAMPLE

glct_region_tree = ape::nj(ape::dist.gene(glct_region))
ggtree(glct_region_tree) +
  geom_tiplab(size = 2) +
  # geom_label2(aes(label=bootstrap), size = 2) +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) 



