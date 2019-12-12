library(ggplot2)  # FYI you need v2.0
library(dplyr)    # yes, i could have not done this and just used 'subset' instead of 'filter'
library(ggalt)    # devtools::install_github("hrbrmstr/ggalt")
library(ggthemes) # theme_map and tableau colors
library(cegwas2)
library(tidyverse)
library(ggtree)
library(ape)
library(ggrepel)

get_vcf2 <- function() {
  # Function for fetching the path the VCF
  path <- glue::glue("~/Dropbox/Andersenlab/Reagents/WormReagents/_SEQ/WI/WI-{cendr_dataset_release}/vcf/WI.{cendr_dataset_release}.snpeff.vcf.gz") # nolint
  if (file.exists(path)) {
    message("Using local VCF")
  } else {
    message("Using remote VCF")
    path <- glue::glue("http://storage.googleapis.com/elegansvariation.org/releases/{cendr_dataset_release}/variation/WI.{cendr_dataset_release}.soft-filter.vcf.gz") # nolint
  }
  path
}

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
setwd("..")



# modify settings
background_color <- "white"
axes_text_size <- "black"
axis_color <- "black"
axes_text_size <- 10

# load world map and remove antartica 
world <- map_data("world")
world <- world[world$region != "Antarctica",] # intercourse antarctica

# download data from CeNDR
isolation_info <- readr::read_tsv("https://elegansvariation.org/strain/strain_data.tsv")

# filter specific strains
strains_to_plot <- c("N2", "CB4856")

cendr_dataset_release <- "20180527"
allele_of_interest <- cegwas2::query_vcf("glct-3", vcf = get_vcf2()) 

allele_of_interest <- allele_of_interest %>%
  dplyr::filter(aa_change == "p.Gly16*") %>%
  dplyr::select(strain = SAMPLE, REF, ALT, a1, a2) %>%
  dplyr::mutate(a_col = ifelse(a1 == "T" | a2 == "T", "ALT", "REF"))

g3 <- cegwas2::query_vcf(c("glct-3"), vcf = get_vcf2()) 


alt_strains <-g3%>%
  dplyr::rowwise() %>%
  dplyr::filter( impact == "HIGH" & (a1==ALT || a2 == ALT || grepl(a1, ALT) )) %>% 
  dplyr::filter(SAMPLE!="ECA701" & query == "glct-3") %>%
  dplyr::select(strain = SAMPLE, gene_name, aa_change, REF, ALT, a1, a2) %>%
  dplyr::mutate(a_col = ifelse((a1 == "T" | a2 == "T") & gene_name == "glct-3", "ALT", 
                               ifelse(a1 == ALT | a2 == ALT,"ALT","REF")) )
  


strains_330 <- isolation_info%>%
  dplyr::filter(reference_strain == "True")%>%
  dplyr::select(strain = isotype, long = longitude, lat = latitude) %>%
  dplyr::filter(strain%in% alt_strains$strain) %>%
  dplyr::left_join(., alt_strains, by = "strain")  


ggplot()+ 
  geom_map(data=world, map=world,
           aes(x=long, y=lat, map_id=region),
           color=axis_color, fill=background_color, size=0.5)+
  geom_point(data = strains_330 %>% dplyr::filter(a_col == "ALT"), 
             aes(x=as.numeric(long), y=as.numeric(lat), fill = a_col), 
             color = "black",
             shape = 21, 
             size = 4) +
  geom_label_repel(
    data          = strains_330 %>% dplyr::filter(a_col == "ALT"),
    segment.size  = 0.2,
    segment.color = "grey50",
    aes(label = strain, x=as.numeric(long), y=as.numeric(lat))
    # direction = "x"
  ) +
  scale_fill_manual(values = c("cadetblue3","hotpink3"))+
  theme_map()+
  labs(fill = "Gly16*") +
  theme(panel.background = element_rect(fill = background_color, colour = NA),
        text = element_text(size = axes_text_size),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.key.size = unit(1, "cm")) 

ggsave("Plots/world_distribution.pdf",height = 10, width = 20)
ggsave("Plots/world_distribution.svg",height = 10, width = 20)




strains_330 <- isolation_info%>%
  dplyr::filter(reference_strain == "True")%>%
  dplyr::select(strain = isotype, long = longitude, lat = latitude, elevation:substrate) %>%
  # dplyr::filter(strain%in% alt_strains$strain) %>%
  dplyr::left_join(., alt_strains, by = "strain")  


strains_330$gene_name[is.na(strains_330$gene_name)] <- "ref"

strains_330 <- strains_330 %>%
  dplyr::select(strain, substrate, gene_name) %>%
  dplyr::mutate(gene_name = ifelse(strain == "ECA369" |strain == "ECA347" | strain == "MY10" | gene_name == "glct-3", "glct-3", NA)) %>%
  dplyr::mutate(val = 1) %>%
  dplyr::filter(substrate!="None") 

strains_330%>%
  ggplot() + 
  geom_bar(aes(y = val, x = gene_name, fill = substrate), position="fill", stat="identity", color = "black", size = 0.1)+
  scale_fill_viridis_d()



# load tree
tree <- ape::read.tree(glue::glue("Data/whole_genome_tree/330_genome.raxml.bestTree"))

# highlight branches for strains of interest
branch_strains <- list(CONNECT = alt_strains$strain)

tree_pt_h <- ggtree::groupOTU(tree, branch_strains)

ggtree(tree_pt_h,
       branch.length="rate", 
       aes(color=group)) + 
  scale_color_manual(values=c("hotpink3", "cadetblue3"), 
                     name = "GLCT-3\nAllele", 
                     labels=c("REF", "Gly16*")) + 
  theme(legend.position="right")+
  theme_tree2() 

ggsave("Plots/genomewide_tree.pdf",height = 32, width = 12)
ggsave("Plots/genomewide_tree.svg",height = 32, width = 12)

# basic tree - circle
tree_pt<-ggtree(tree,
                branch.length="rate", layout = "circular")

clean_strains <- strains_in_tree %>%
  dplyr::filter(highlight_strains == "high")

tree_pt %<+% clean_strains + 
  geom_tiplab2(aes(color=highlight_strains))+
  scale_color_manual(values=c("high" =highlight_color , "not" = background_color))



# glct3 <- c(12385766, 12388791)
# 
# mcolor_grp <- plot_df %>% dplyr::select(haplotype, color) %>% dplyr::distinct()
# mcolor <- mcolor_grp$color
# names(mcolor) <- mcolor_grp$haplotype
# 
# strain_labels <- plot_df %>%
#   dplyr::mutate(alt_glct = ifelse(isotype %in% alt_strains, "ALT", "REF")) %>%
#   dplyr::distinct(isotype, plotpoint, alt_glct) %>%
#   dplyr::arrange(alt_glct, (plotpoint)) %>%
#   dplyr::mutate(plotpoint = row_number())
# 
# plotdf <- plot_df %>%
#   dplyr::select(-plotpoint) %>%
#   dplyr::left_join(.,strain_labels, by = "isotype")
# 
# plotdf %>%
#   dplyr::filter(chromosome == "I") %>%
#   dplyr::mutate(alt_glct = ifelse(isotype %in% alt_strains, "ALT", "REF")) %>%
#   ggplot(.,
#          aes(xmin = start/1E6, xmax = stop/1E6,
#              ymin = plotpoint - 0.5, ymax = plotpoint + 0.5,
#              fill = haplotype)) +
#   geom_rect() +
#   scale_fill_manual(values = mcolor) +
#   scale_y_continuous(breaks = strain_labels$plotpoint,
#                      labels = strain_labels$isotype,
#                      expand = c(0, 0)) +
#   xlab("Position (Mb)") +
#   theme_bw() +
#   coord_cartesian(xlim=c(12, 13)) +
#   geom_vline(aes(xintercept = glct3[1]/1e6)) +
#   geom_vline(aes(xintercept = glct3[2]/1e6)) +
#   facet_grid(alt_glct~chromosome, scales="free_y", space="free_y") +
#   theme(legend.position="none",
#         axis.text.y = element_blank(), 
#         axis.ticks.y = element_blank())
# 
# ggsave("Plots/haplotype.png",height = 12, width = 8)


allele_of_interest <- cegwas2::query_vcf(c("glct-3","glct-1","glct-2","glct-4","glct-5","glct-6"), vcf = get_vcf2()) 

allele_of_interest%>%
  dplyr::rowwise() %>%
  dplyr::filter( impact == "HIGH", grepl(a1, ALT), a1!=REF)%>%
  dplyr::group_by(query) %>%
  dplyr::summarise(ct = n())

allele_of_interest%>%
  dplyr::rowwise() %>%
  dplyr::filter( impact == "HIGH", grepl(a1, ALT), a1!=REF)%>%
  dplyr::group_by(SAMPLE) %>%
  dplyr::summarise(ct = n()) %>%
  View()

glct_alts <-allele_of_interest%>%
  dplyr::rowwise() %>%
  dplyr::filter( impact == "HIGH" & (a1==ALT || a2 == ALT || grepl(a1, ALT) ), a1!=REF) 

write.table(glct_alts, "Processed_Data/glct_gene_high_variant.txt", col.names = T, row.names = F, quote = F, sep = " ")
