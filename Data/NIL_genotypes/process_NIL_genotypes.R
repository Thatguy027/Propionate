setwd("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/NIL_genotypes/")
nana_nils <- data.table::fread("gt_hmm.tsv")

length(unique(nana_nils$sample))
# df - nil genotype df in segment format
# chr - NIL chromosome
# left.cb - if looking at NIL breakups, this number corresponds to the left flank for CB4856 NILs
# left.n2 - if looking at NIL breakups, this number corresponds to the left flank for N2 NILs
# left.bound - left boundary for region of interest
# right.bound - right boundary for region of interest
# scan.range - cutoff for small genomic regions when identifying left nils
# qtl.left - where to draw red line #1
# qtl.right - where to draw red line #2
# qtl. peak - where to draw cyan line

nil_plot <- function(df, chr, left.cb, left.n2 , left.bound,right.bound, scan.range, qtl.left, qtl.right, qtl.peak,all.chr=F){
  # # # determine if NILs are CB or N2
  nilsII_sort_type <- dplyr::filter(df)%>%    
    dplyr::group_by(sample)%>%
    dplyr::mutate(size = end - start )%>%
    dplyr::group_by(sample, start)%>%
    dplyr::mutate(gt_ct = sum(size))%>%
    dplyr::ungroup()%>%
    dplyr::group_by(sample, gt)%>%
    dplyr::mutate(major_gt = sum(gt_ct))%>%
    dplyr::ungroup()%>%
    dplyr::arrange(desc(major_gt))%>%
    dplyr::distinct(sample, .keep_all = T)%>%
    dplyr::mutate(nil_type = ifelse(gt == 2, "DL", "BRC"))%>%
    dplyr::select(sample, nil_type)
  
  # # # keep NILs that lost NIL genotype on right side
  nilsII_left <- dplyr::filter(df)%>%    
    dplyr::filter(chrom == chr)%>%
    dplyr::left_join(.,nilsII_sort_type, by = "sample")%>%
    dplyr::filter((start > left.cb - .3e6 & start < left.cb + .3e6 & nil_type == "DL") |
                    (start > left.n2 - .3e6 & start < left.n2 + .3e6 & nil_type == "BRC") )%>%
    dplyr::filter(gt == 2 & nil_type == "DL" | gt == 1 & nil_type == "BRC")%>%
    dplyr::mutate(side = "LEFT")%>%
    dplyr::mutate(size = end - start )%>%
    dplyr::ungroup()%>%
    dplyr::arrange(desc(size))%>%
    dplyr::distinct(sample, .keep_all = T)%>%
    dplyr::filter(size > scan.range)%>% # # # remove small (likely wrong calls) around interval site
    dplyr::select(sample, side)
  
  # # # keep NILs that lost NIL genotype on left side
  nilsII_right <- dplyr::filter(df)%>%    
    dplyr::filter(chrom == chr)%>%
    dplyr::left_join(.,nilsII_sort_type, by = "sample")%>%
    dplyr::filter(!sample %in% nilsII_left$sample)%>%
    dplyr::mutate(side = "RIGHT")%>%
    dplyr::distinct(sample,.keep_all = T)%>%
    dplyr::select(sample, side)
  
  nil_sides <- bind_rows(nilsII_left,nilsII_right)
  
  
  nilsII_sort_left <- dplyr::filter(df)%>%    
    dplyr::filter(chrom == chr)%>%
    dplyr::left_join(.,nilsII_sort_type, by = "sample")%>%
    dplyr::left_join(.,nil_sides, by = "sample")%>%
    dplyr::filter(gt == 2 & nil_type == "DL" | gt == 1 & nil_type == "BRC")%>%
    dplyr::group_by(sample)%>%
    dplyr::filter(start > left.bound & start < right.bound)%>%
    dplyr::filter(end > left.bound | end < right.bound)%>%
    # dplyr::filter( ((start > left.n2 -.1e6 | start < left.n2 +.1e6)  & start < right.bound & nil_type == "DL")|
    #                  ((start > left.cb -.1e6 | start < left.cb +.1e6) & start < right.bound & nil_type == "BRC"))%>%
    dplyr::filter(end > left.bound | end < right.bound)%>%
    dplyr::mutate(size = end - start )%>%
    dplyr::group_by(sample, start)%>%
    dplyr::mutate(gt_ct = sum(size))%>%
    dplyr::ungroup()%>%
    dplyr::filter(side == "LEFT")%>%
    dplyr::arrange( desc(gt_ct))%>%
    dplyr::distinct(sample, .keep_all = T)%>%
    dplyr::arrange(nil_type, desc(gt_ct))
  
  nilsII_sort_right <- dplyr::filter(df)%>%    
    dplyr::filter(chrom == chr)%>%
    dplyr::left_join(.,nilsII_sort_type, by = "sample")%>%
    dplyr::left_join(.,nil_sides, by = "sample")%>%
    dplyr::filter(side == "RIGHT")%>%
    dplyr::filter(gt == 2 & nil_type == "DL" | gt == 1 & nil_type == "BRC")%>%
    dplyr::group_by(sample)%>%
    dplyr::filter(start > left.bound & start < right.bound)%>%
    dplyr::filter(end > left.bound | end < right.bound)%>%
    dplyr::mutate(size = end - start )%>%
    dplyr::group_by(sample, start)%>%
    dplyr::mutate(gt_ct = sum(size))%>%
    dplyr::ungroup()%>%
    dplyr::arrange( desc(gt_ct))%>%
    dplyr::distinct(sample, .keep_all = T)%>%
    dplyr::ungroup()%>%
    dplyr::arrange( nil_type, gt_ct)
  
  nilsII_sort <- bind_rows(nilsII_sort_right, nilsII_sort_left)
  
  if ( all.chr == T){
    
    nilsII <- dplyr::filter(df)%>%    
      dplyr::filter(chrom != "MtDNA")%>%
      dplyr::left_join(.,nilsII_sort_type, by = "sample")
    
    nilsII$sample <- factor(nilsII$sample, levels = unique(nilsII_sort$sample), labels = unique(nilsII_sort$sample), ordered = T)
    nilsII$gt <- as.character(nilsII$gt)
    
    nl.pl <- ggplot(nilsII)+
      geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt, size = 2))+
      facet_grid(nil_type~chrom, scales = "free",  space = "free_y")+
      scale_color_manual(values=c("1"="cadetblue3","2"="hotpink3"))+
      theme(axis.text.x = element_text(size=12, face="bold", color="black"),
            axis.text.y = element_text(size=12, face="bold", color="black"),
            axis.title.x = element_text(size=14, face="bold", color="black"),
            axis.title.y = element_text(size=14, face="bold", color="black"),
            strip.text.x = element_text(size=20,face="bold", color="black"),
            strip.text.y = element_text(size=12, angle =0, face="bold", color="black"),
            plot.title = element_text(size=24, face="bold"),
            legend.position = "none")+
      labs(x = "Genomic Position (Mb)", y = "NIL")
    
  }else{
    nilsII <- dplyr::filter(df)%>%    
      dplyr::filter(chrom == chr, chrom != "MtDNA")%>%
      dplyr::left_join(.,nilsII_sort_type, by = "sample")
    
    nilsII$sample <- factor(nilsII$sample, levels = unique(nilsII_sort$sample), labels = unique(nilsII_sort$sample), ordered = T)
    nilsII$gt <- as.character(nilsII$gt)
    
    nl.pl <- ggplot(nilsII)+
      geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt, size = 2))+
      facet_grid(nil_type~chrom, scales = "free",  space = "free_y")+
      scale_color_manual(values=c("1"="cadetblue3","2"="hotpink3"))+
      theme(axis.text.x = element_text(size=12, face="bold", color="black"),
            axis.text.y = element_text(size=12, face="bold", color="black"),
            axis.title.x = element_text(size=14, face="bold", color="black"),
            axis.title.y = element_text(size=14, face="bold", color="black"),
            strip.text.x = element_text(size=20,face="bold", color="black"),
            strip.text.y = element_text(size=12, angle =0, face="bold", color="black"),
            plot.title = element_text(size=24, face="bold"),
            legend.position = "none")+
      labs(x = "Genomic Position (Mb)", y = "NIL")
  }
  
  return(list(nl.pl,nilsII_sort, nilsII))
}

pt_nils <- dplyr::filter(nana_nils, sample %in% c("Nana_B9_E09", "Nana_B6_E04", "Nana_C6_C10", "Nana_B2_G01", "Nana_A6_E08", "Nana_B8_E05"))

output <- nil_plot(pt_nils, "V", 3.6e6, 4.8e6, 3e6, 5e6, 3e4, 3.1, 4.3, 4.03,T )

output[[1]] 

ggsave("whole_genome_genotypes.pdf",height = 2, width = 12)

output <- nil_plot(pt_nils, "V", 3.6e6, 4.8e6, 3e6, 5e6, 3e4, 3.1, 4.3, 4.03,F )

output[[1]] 

ggsave("chrV_genotypes.pdf",height = 2, width = 12)

output[[1]] + 
  xlim(c(3,5))+
  theme_bw()+
  theme(legend.position = "none")

ggsave("chrV_zoom_genotypes.pdf",height = 2, width = 12)
output[[2]] 
output[[3]] 
