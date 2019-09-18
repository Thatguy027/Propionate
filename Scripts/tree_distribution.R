library(ggplot2)  # FYI you need v2.0
library(dplyr)    # yes, i could have not done this and just used 'subset' instead of 'filter'
library(ggalt)    # devtools::install_github("hrbrmstr/ggalt")
library(ggthemes) # theme_map and tableau colors
library(cegwas2)
library(tidyverse)
library(ggtree)
library(ape)

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

strains_330 <- isolation_info%>%
  dplyr::filter(reference_strain == "True")%>%
  dplyr::select(strain = isotype, long = longitude, lat = latitude) %>%
  dplyr::filter(strain%in% allele_of_interest$strain) %>%
  dplyr::left_join(., allele_of_interest, by = "strain")  

ggplot()+ geom_map(data=world, map=world,
                   aes(x=long, y=lat, map_id=region),
                   color=axis_color, fill=background_color, size=0.5)+
  geom_point(data = strains_330 %>% dplyr::filter(a_col == "REF"), 
             aes(x=as.numeric(long), y=as.numeric(lat), fill = a_col), 
             color = "black",
             shape = 21, 
             size = 3) +
  geom_point(data = strains_330 %>% dplyr::filter(a_col == "ALT"), 
             aes(x=as.numeric(long), y=as.numeric(lat), fill = a_col), 
             color = "black",
             shape = 21, 
             size = 3) +
  scale_fill_manual(values = c("cadetblue3","hotpink3"))+
  theme_map()+
  labs(fill = "Gly16*") +
  theme(panel.background = element_rect(fill = background_color, colour = NA),
        text = element_text(size = axes_text_size),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.key.size = unit(1, "cm")) 


ggsave("Plots/world_distribution.pdf",height = 5, width = 10)

alt_strains <- dplyr::filter(strains_330, a_col == "ALT") %>%
  dplyr::pull(strain)

# load tree
tree <- ape::read.tree(glue::glue("Data/whole_genome_tree/330_genome.raxml.bestTree"))

# highlight branches for strains of interest
branch_strains <- list(CONNECT = alt_strains)

tree_pt_h <- ggtree::groupOTU(tree, branch_strains)

ggtree(tree_pt_h,
       branch.length="rate", 
       aes(color=group)) + 
  geom_tiplab(align = T) +
  scale_color_manual(values=c("hotpink3", "cadetblue3"), 
                     name = "Presence of TALT", 
                     breaks=c("0", "TALT"),
                     labels=c("FALSE", "TRUE")) + 
  theme(legend.position="right")+
  theme_tree2() 

ggsave("Plots/genomewide_tree.pdf",height = 32, width = 12)


# basic tree - circle
tree_pt<-ggtree(tree,
                branch.length="rate", layout = "circular")

clean_strains <- strains_in_tree %>%
  dplyr::filter(highlight_strains == "high")

tree_pt %<+% clean_strains + 
  geom_tiplab2(aes(color=highlight_strains))+
  scale_color_manual(values=c("high" =highlight_color , "not" = background_color))