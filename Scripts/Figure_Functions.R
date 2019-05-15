background_color <- "white"

Plot_Peak_LD <- function(processed_mapping, genotype_matrix){
  snp_df <- processed_mapping %>% na.omit()
  ld_snps <- dplyr::filter(genotype_matrix, CHROM %in% snp_df$CHROM, POS %in% snp_df$POS)
  ld_snps <- data.frame(snp_id = paste(ld_snps$CHROM, ld_snps$POS, sep = "_"), data.frame(ld_snps)[, 5:ncol(ld_snps)])
  sn <- list()
  for (i in 1:nrow(ld_snps)) {
    sn[[i]] <- genetics::genotype(as.character(gsub(1, "T/T", gsub(-1, "A/A", ld_snps[i, 4:ncol(ld_snps)]))))
  }
  
  test <- data.frame(sn)
  colnames(test) <- (ld_snps$snp_id)
  ldcalc <- t(genetics::LD(test)[[4]])^2
  diag(ldcalc) <- 1
  LDs <- tbl_df(data.frame(ldcalc) %>% dplyr::add_rownames(var = "SNP1")) %>% 
    tidyr::gather(SNP2, r2, -SNP1) %>% dplyr::arrange(SNP1) %>% 
    tidyr::separate(SNP1, sep = "_", into = c("CHROM1", "POS1"), remove = F) %>% 
    dplyr::arrange(CHROM1,  as.numeric(POS1))
  
  ldplot <- ggplot2::ggplot(LDs) + 
    ggplot2::aes(x = factor(SNP1, levels = unique(SNP1), ordered = T), 
                 y = factor(SNP2, levels = unique(SNP1), ordered = T)) + 
    ggplot2::geom_tile(ggplot2::aes(fill = r2), color = background_color) + 
    ggplot2::geom_text(ggplot2::aes(label = signif(r2, 3)), fontface = "bold", size = 4.5) + 
    scale_x_discrete(labels = function(x) {
      gsub("_", ":", x)
    }, expand = c(0, 0)) + 
    scale_y_discrete(position = "right",  
                     labels = function(x) {gsub("_", ":", x)}, 
                     expand = c(0, 0)) + 
    scale_fill_continuous(high = "#FF0000", low = "white", na.value = "white")
  
  return(list(ldplot, LDs))
}
