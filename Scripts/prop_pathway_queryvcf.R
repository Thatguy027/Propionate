prop_pathway_genes <- c("WBGene00017864","WBGene00018701","WBGene00008415","WBGene00014202",
                        "WBGene00016943","WBGene00001155","WBGene00017301","WBGene00012608",
                        "WBGene00000114","WBGene00018138","WBGene00020169","WBGene00016144",
                        "WBGene00006510","WBGene00010988","WBGene00015512","WBGene00013577",
                        "WBGene00003411","WBGene00004062","WBGene00013855","WBGene00003609","WBGene00003658")
  
pp_variant <- list()
ct <- 1
for(g in prop_pathway_genes){
  pp_variant[[ct]] <- cegwas2::query_vcf(g, impact = "ALL")
  ct <- 1+ct
}  

pp_variant <- dplyr::bind_rows(pp_variant)
pp_variant_high <- dplyr::filter(pp_variant, SAMPLE == "DL238", a1 == ALT | a2 == ALT)

pp_variant_high <- dplyr::filter(pp_variant, impact == "HIGH" | impact == "MODERATE", a1 == ALT | a2 == ALT)
