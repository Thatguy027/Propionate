library(data.table)
library(tidyverse)
library(cegwas)
library(cegwas2)
library(broom)
library(cowplot)

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
setwd("..")

source("Scripts/Figure_Functions.R")


c2_prmaps <- data.table::fread("Processed_Data/GWAS_processed_mapping_cegwas2.tsv")
genes <- cegwas2::query_vcf(c("I:12374204-12388791"))
gm <- readr::read_tsv(glue::glue("Processed_Data/Genotype_Matrix.tsv"))


skat_maps <- data.table::fread("Processed_Data/GWAS_Skat.assoc") %>%
  dplyr::rowwise()%>%
  dplyr::mutate(CHROM = strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][1],
                POS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][1]),
                endPOS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][2]))%>%
  dplyr::mutate(size = abs(POS-endPOS)) %>%
  dplyr::filter(CHROM!="MtDNA")%>%
  dplyr::ungroup()%>%
  dplyr::filter(NumVar > 1)%>%
  dplyr::mutate(significant = ifelse(Pvalue < .05/n(), TRUE,FALSE ))%>%
  dplyr::arrange(Pvalue)

########

norm_pheno <- na.omit(c2_prmaps) %>%
  dplyr::distinct(strain, value) %>%
  dplyr::ungroup()%>%
  dplyr::arrange(value) %>%
  dplyr::mutate(norm_pheno_temp = ifelse(value == min(value), 0, 1))%>%
  dplyr::mutate(delta_pheno = ifelse(norm_pheno_temp == 0, 0, abs(dplyr::lag(value) - value)))%>%
  dplyr::mutate(norm_pheno = cumsum(delta_pheno)) %>%
  dplyr::mutate(final_pheno = norm_pheno/max(norm_pheno)) %>%
  dplyr::select(strain, final_pheno)

gene_pr <- dplyr::filter(genes, SAMPLE %in% unique(c2_prmaps$strain)) %>%
  dplyr::select(CHROM, POS, strain = SAMPLE, REF, ALT, gene_id, gene_name, allele = a1, impact, effect, aa_change) %>%
  dplyr::mutate(TGT = ifelse(allele == ALT, "ALT", "REF")) %>%
  dplyr::left_join(., norm_pheno, by = "strain") %>%
  dplyr::group_by(strain, gene_id) %>%
  dplyr::mutate(vt_ct = sum(grepl("ALT",TGT)))

# my10 does not have glct and is divergent in the area
# JU1395 has a frame shift
gene_ref <- dplyr::filter(gene_pr,  vt_ct == 0, strain != "MY10" ) %>%
  dplyr::mutate(rmju = ifelse(gene_name == "glct-3" & strain == "JU1395", "remove", "keep")) %>%
  dplyr::filter(rmju != "remove") %>%
  dplyr::distinct(strain, final_pheno, gene_id, .keep_all = T) %>%
  dplyr::mutate(plot_gt = "REF")

gene_high <- dplyr::filter(gene_pr, impact == "HIGH", TGT == "ALT" | (strain == "JU1395" & gene_name == "glct-3")) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::arrange(impact) %>%
  dplyr::filter(impact == "HIGH") %>%
  dplyr::mutate(plot_gt = effect) 

gene_alt <- dplyr::filter(gene_pr, TGT == "ALT") %>%
  # dplyr::filter(!strain %in% gene_high$strain | strain %in% c("BRC20067", "DL238")) %>%
  dplyr::group_by(strain, gene_id) %>%
  dplyr::mutate(plot_gt = paste(unique(aa_change), collapse = " ")) %>%
  dplyr::mutate(plot_gt = gsub("p\\.","", plot_gt)) %>%
  dplyr::distinct(strain, final_pheno, gene_id,plot_gt, .keep_all =T)

sig_genes <- skat_maps%>%
  dplyr::filter(significant == T, grepl("I:",RANGE))%>% dplyr::pull(Gene)

gene_all <- dplyr::bind_rows(list(gene_ref,gene_high,gene_alt)) %>%
  dplyr::filter(gene_id %in% sig_genes)



countSpaces <- function(s) { sapply(gregexpr(" ", s), function(p) { sum(p>=0) } ) }

gene_plots <- list()
for(g in 1:length(unique(gene_all$gene_id))){
  
  temp_all <- gene_all %>% 
    dplyr::filter(gene_id == unique(gene_all$gene_id)[g]) %>%
    dplyr::ungroup()
  
  temp_all$ids <- temp_all %>%
    dplyr::group_by(plot_gt) %>%
    group_indices(., plot_gt) 
  
  temp_all <- temp_all %>%
    dplyr::rowwise() %>%
    dplyr::mutate(alt_ct = countSpaces(plot_gt)) %>%
    dplyr::mutate(alt_ct = ifelse(plot_gt!="REF" & !plot_gt %in% gene_high$plot_gt, paste0("Hap:", ids), plot_gt)) %>%
    dplyr::mutate(alt_ct = gsub("\\&","\n",alt_ct)) %>%
    dplyr::group_by(alt_ct, gene_id) %>%
    dplyr::mutate(mean_gt = median(final_pheno))
  
  temp_df <- temp_all %>%
    dplyr::ungroup()%>%
    dplyr::arrange(desc(mean_gt)) %>%
    dplyr::mutate(pt_gt = factor(alt_ct, levels = unique(c(alt_ct[!alt_ct %in% "REF"], "REF")))) %>%
    dplyr::mutate(pt_col = factor(strain, 
                                  levels = c("BRC20067", "DL238", unique(strain)[!unique(strain) %in% c("BRC20067", "DL238")]) ,
                                  labels = c("BRC20067", "DL238", rep("other", length(unique(strain)[!unique(strain) %in% c("BRC20067", "DL238")]))))) %>%
    dplyr::distinct(strain, final_pheno, alt_ct,pt_col,.keep_all = T)
  
  gene_plots[[g]] <- ggplot(temp_df)+
    ggbeeswarm::geom_beeswarm(aes(x =pt_gt, y = final_pheno, 
                                  fill = pt_col, shape = pt_col, size = pt_col))+
    scale_fill_manual(values = c("hotpink3", "cadetblue3", "gray50"))+
    scale_size_manual(values = c(5, 5, 2))+
    scale_shape_manual(values = c(23, 23, 21))+
    geom_point(aes(x = pt_gt, y = mean_gt), 
               data = dplyr::distinct(temp_df, alt_ct, gene_id, mean_gt,.keep_all=T),
               shape = 23, 
               size = 3,
               fill = "red")+
    theme_bw(12) +
    facet_wrap(~gene_name, scale = "free_y") +
    labs(x = "Genotype",
         y = "Normalized L1 Survival") +
    theme(legend.position = "none",strip.text = element_text(face = "italic"))+
    coord_flip()
}

length(gene_plots)

all_gene_plot <- cowplot::plot_grid(gene_plots[[1]]+ theme(axis.title.x = element_blank()) ,
                                    gene_plots[[2]]+ theme(axis.title.y = element_blank(),axis.title.x = element_blank()),
                                    gene_plots[[3]],
                                    gene_plots[[4]]+theme(axis.title.y = element_blank()), 
                                    ncol = 2, 
                                    label_size = 14, align = "vh", labels = "AUTO")
all_gene_plot


load("/Users/Stefan/github_repos/Propionate/vcf/gene_ref_flat.Rda")

gene_df <- gene_ref_flat %>%
  dplyr::filter(wbgene %in% c(gene_pr$gene_id, "WBGene00011781", "WBGene00011650","WBGene00008479","WBGene00008160","WBGene00008293")) %>%
  dplyr::select(gene_id = wbgene, strand, txstart, txend, feature_id = gene) %>%
  dplyr::arrange(txstart, feature_id)%>%
  dplyr::distinct(gene_id, feature_id, .keep_all = TRUE) %>%
  dplyr::distinct(gene_id, peak_marker, CHROM, strand, txstart, txend, start_pos, end_pos) 



roi_ld <- data.table::fread("vcf/propionate_chrI_d.ld")

glct_stop <- "I:12385811"

pr_roi_ld <- roi_ld %>%
  dplyr::mutate(peak_marker = gsub("_", ":", glct_stop)) %>%
  dplyr::mutate(marker = gsub(":", "_", SNP_B)) %>%
  dplyr::select(peak_marker, peak_maf = MAF_A, marker, maf_marker_b = MAF_B, ld_r2 = R2)  %>%
  tidyr::separate(marker, c("CHROM","POS"),sep ="_", remove =F, convert = T)

peak_roi_marker <- dplyr::filter(pr_roi_ld, marker == "I_12385811")

peak_roi_marker

ggplot(pr_roi_ld) +
  aes(x = POS/1e6) +
  geom_point(aes(fill = ld_r2, y = ld_r2), shape = 23, size = 3) +
  geom_point(aes(y = ld_r2), shape = 23, size = 3, fill = "red",
             data = peak_roi_marker) +
  scale_fill_viridis_c(name = "R2") +
  theme_bw(15)+
  labs(x = "Genomic Position (Mb)",
       y = expression(-log[10](italic(p)))) + xlim(12.3,12.46)

ggplot(pr_roi_ld) +
  geom_rect(aes(xmin = ifelse(strand == "+", txstart/1e6, txend/1e6),
                xmax = ifelse(strand == "+", txend/1e6, txstart/1e6),
                ymin = 1.05,
                ymax = 0), size = 1, data = gene_df, alpha = 0.25) +
  geom_point(aes(x = POS/1e6, fill = ld_r2, y = ld_r2), shape = 23, size = 3) +
  geom_point(aes(x = POS/1e6, y = ld_r2), shape = 23, size = 3, fill = "red",
             data = peak_roi_marker) +
  scale_fill_viridis_c(name = "R2") +
  theme_bw(15)+
  geom_segment(aes(x = ifelse(strand == "+", txstart/1e6, txend/1e6),
                   xend = ifelse(strand == "+", txend/1e6, txstart/1e6),
                   y = 1.05,
                   yend = 1.05),
               arrow = arrow(length = unit(5, "points")), size = 1, data = gene_df) +

  labs(x = "Genomic Position (Mb)",
       y = expression(-log[10](italic(p)))) + xlim(12.3,12.46)




genes_to_ld <- dplyr::filter(genes, SAMPLE %in% unique(c2_prmaps$strain)) %>%
  dplyr::select(CHROM, POS, strain = SAMPLE, REF, ALT, gene_id, gene_name, allele = a1, impact, effect, aa_change) %>%
  dplyr::mutate(TGT = ifelse(allele == ALT, "ALT", "REF")) %>%
  dplyr::select(CHROM, POS, strain, TGT) %>%
  tidyr::unite(marker, CHROM, POS) %>%
  dplyr::mutate(TGT = ifelse(TGT == "REF", -1, 1)) %>%
  dplyr::distinct() %>%
  tidyr::spread(strain, TGT)
  

genes_to_ld[5596:5598,]

gene_pr %>%
  dplyr::ungroup() %>%
  dplyr::select(CHROM, POS, strain, TGT) %>%
  tidyr::unite(marker, CHROM, POS) %>%
  dplyr::mutate(TGT)



sn <- list()
for (i in 1:nrow(genes_to_ld)) {
  sn[[i]] <- genetics::genotype(as.character(gsub(1, "T/T", gsub(-1, "A/A", genes_to_ld[i, 4:ncol(genes_to_ld)]))))
}

test <- data.frame(sn)
colnames(test) <- (genes_to_ld$marker)
ldcalc <- t(genetics::LD(test)[[4]])^2
diag(ldcalc) <- 1
LDs <- tbl_df(data.frame(ldcalc) %>% dplyr::add_rownames(var = "SNP1")) %>% 
  tidyr::gather(SNP2, r2, -SNP1) %>% dplyr::arrange(SNP1) %>% 
  tidyr::separate(SNP1, sep = "_", into = c("CHROM1", "POS1"), remove = F) %>% 
  dplyr::arrange(CHROM1,  as.numeric(POS1))

ldplot <- ggplot2::ggplot(LDs) + 
  ggplot2::aes(x = factor(SNP1, levels = unique(SNP1), ordered = T), 
               y = factor(SNP2, levels = unique(SNP1), ordered = T)) + 
  ggplot2::geom_tile(ggplot2::aes(fill = r2), color = NA) + 
  scale_x_discrete(labels = function(x) {
    gsub("_", ":", x)
  }, expand = c(0, 0)) + 
  scale_y_discrete(position = "right",  
                   labels = function(x) {gsub("_", ":", x)}, 
                   expand = c(0, 0)) + 
  scale_fill_continuous(high = "#FF0000", low = "white", na.value = "white")+
  labs(fill = expression("R"^2))





# 
# ggsave(plot = all_gene_plot, filename = "Plots/chr1_gene_geno_pheno.pdf", height = 12, width = 20)
# ggsave(plot = all_gene_plot, filename = "Plots/chr1_gene_geno_pheno.png", height = 12, width = 20, dpi = 300)
# ggsave(plot = all_gene_plot, filename = "Plots/chr1_gene_geno_pheno.svg", height = 12, width = 20)