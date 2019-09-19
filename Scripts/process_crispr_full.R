library(tidyverse)

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
setwd("..")

phen <- data.table::fread("Data/glct_causality/glct_deletion_expt_full.csv")

phen %>%
  ggplot()+
  aes(x = Strain, y = Survival, fill = Strain)+
  geom_boxplot(alpha = 0.8)+
  scale_fill_manual(values=c("hotpink3","gray60","gray60","cadetblue3"))+
  scale_x_discrete(labels=c("BRC20067" = "BRC20067", 
                            "DL238" = "DL238", 
                            "BRC20067_Del" = expression("BRC20067 "~ Delta~italic("glct-3")),
                            "BRC20067_Swap" = "BRC20067 (Gly16*)"))+
  ggbeeswarm::geom_beeswarm(cex = 0.7, priority = "density", aes(color = Day))+
  theme_classic(20)+
  labs(y = "L1 Survival (%)") +
  # ylim(0,100)+
  theme(axis.title.x = element_blank(), legend.position = "none")


reg_pheno <- phen %>%
  dplyr::do(broom::augment(lm(Survival ~ Day, data =.)))%>%
  dplyr::select(-Survival) %>%
  dplyr::rename(Survival = .resid) %>%
  dplyr::mutate(Strain = phen$Strain) %>%
  dplyr::select(Strain, Survival) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Survival) %>%
  dplyr::mutate(norm_pheno_temp = ifelse(Survival == min(Survival), 0, 1))%>%
  dplyr::mutate(delta_pheno = ifelse(norm_pheno_temp == 0, 0, abs(dplyr::lag(Survival) - Survival)))%>%
  dplyr::mutate(norm_pheno = cumsum(delta_pheno)) %>%
  dplyr::mutate(final_pheno = norm_pheno/max(norm_pheno)) %>%
  dplyr::select(-Survival) %>%
  dplyr::rename(Survival = final_pheno)
  
reg_pheno%>%
  ggplot()+
  aes(x = Strain, y = Survival, fill = Strain)+
  geom_boxplot(alpha = 0.8)+
  scale_fill_manual(values=c("hotpink3","gray60","gray60","cadetblue3"))+
  scale_x_discrete(labels=c("BRC20067" = "BRC20067", 
                            "DL238" = "DL238", 
                            "BRC20067_Del" = expression("BRC20067 "~ Delta~italic("glct-3")),
                            "BRC20067_Swap" = "BRC20067 (Gly16*)"))+
  ggbeeswarm::geom_beeswarm(cex = 0.7, priority = "density", alpha = 0.5)+
  theme_classic(20)+
  labs(y = "Normalized\nL1 Survival") +
  # ylim(0,100)+
  theme(axis.title.x = element_blank(), legend.position = "none")


ggsave("Plots/Crispr_pheno.pdf", height = 4, width = 10)
ggsave("Plots/Crispr_pheno.png", height = 4, width = 10, dpi = 300)


run_TukeyHSD <- function(df){
  stat_df <- df %>%
    dplyr::select(Strain, phenotype=Survival)
  
  aov_res <- aov(stat_df$phenotype ~ stat_df$Strain)
  summary(aov_res)
  tuk <- TukeyHSD(aov_res)
  
  psig=as.numeric(apply(tuk$`stat_df$Strain`[,2:3],1,prod)>=0)+1
  op=par(mar=c(4.2,9,3.8,2))
  plot(tuk,col=psig,yaxt="n")
  for (j in 1:length(psig)){
    axis(2,at=j,labels=rownames(tuk$`stat_df$Strain`)[length(psig)-j+1],
         las=1,cex.axis=.8,col.axis=psig[length(psig)-j+1])
  }
  par(op)
  
  pwtuk <- TukeyHSD(aov_res)
  
  return(pwtuk)
}
run_TukeyHSD(phen)

run_TukeyHSD(reg_pheno)
