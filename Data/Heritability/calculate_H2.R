library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(data.table)

# x should be a data frame in long format that contains a column "pheno" and a column "strain"

H2.fun <- function(x){
  pdata <- x
  pdata = split(pdata$pheno, pdata$strain)
  pdata.notNAcnt = sapply(pdata, function(x){sum(!is.na(x))})
  pdata[pdata.notNAcnt<2]=NULL
  pdata.melted = melt(pdata)
  names(pdata.melted)=c('pheno', 'strain')
  pdata.melted$strain=as.factor(pdata.melted$strain)
  reffMod = lmer(pheno ~ 1 + (1|strain), data=pdata.melted)
  Var_Random_effect <- as.numeric(VarCorr(reffMod))
  Var_Residual <- attr(VarCorr(reffMod), "sc")^2
  H2 <- Var_Random_effect/(Var_Random_effect+Var_Residual)
  print(H2)
}


raw_df <- fread("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/Heritability/20151022_PA_doseresponse_herit.csv")
raw_df <- fread("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/Combined_Heritability/Data/20151022_PA_doseresponse_herit.csv")

long_df <- raw_df %>%
  gather(strain_n, ratio, -DATE, -PA_conc) 

clean_df <- data.frame(date = long_df$DATE, pa_conc = long_df$PA_conc, strain = long_df$strain_n, pheno = long_df$ratio)%>%
  separate(strain, into = c("strain", "replicate"), sep = "_")

herits <- list()
for(i in 1:length(unique(clean_df$pa_conc))){
  for_calc <- filter(clean_df, pa_conc == unique(clean_df$pa_conc)[i]) %>%
    dplyr::select(strain, pheno)
  
  herits[[i]] <- data.frame(H2 = signif(H2.fun(for_calc),2), pa_conc = unique(clean_df$pa_conc)[i])
  
}

herits_df <- rbindlist(herits)


ggplot(herits_df)+
  aes(x = pa_conc, y = H2)+
  geom_point(size = 4)+
  theme_bw()+
  annotate("text", x = herits_df$pa_conc, y = herits_df$H2+.05, label = herits_df$H2)+
  theme(axis.text.x = element_text(size=24, face="bold", color="black"),
        axis.text.y = element_text(size=24, face="bold", color="black"),
        axis.title.x = element_text(size=24, face="bold", color="black", vjust=-.3),
        axis.title.y = element_text(size=24, face="bold", color="black"),
        strip.text.x = element_text(size=24, face="bold", color="black"),
        strip.text.y = element_text(size=16, face="bold", color="black"),
        plot.title = element_text(size=24, face="bold", vjust = 1),
        legend.position="none",
        panel.background = element_rect( color="black",size=1.2),
        strip.background = element_rect(color = "black", size = 1.2))+
  labs(x = "PA mM", y = "Broad Sense Heritability")

ggsave(filename = "~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/9replicate_H2.pdf",
       height=8,
       width=8)


# sampling 

clean_df$strain <- as.character(clean_df$strain)

sampled_h2 <- list()

for(j in 1:12){
  sampled_df <- clean_df%>%
    dplyr::select(-date)%>%
    group_by(strain, pa_conc) %>%
    sample_n(3)
  
  s_herits <- list()
  for(i in 1:length(unique(sampled_df$pa_conc))){
    for_calc <- filter(sampled_df, pa_conc == unique(sampled_df$pa_conc)[i]) %>%
      dplyr::select(strain, pheno)
    
    s_herits[[i]] <- data.frame(H2 = signif(H2.fun(for_calc),2), pa_conc = unique(sampled_df$pa_conc)[i], sample_n = j)
    
  }
  
  s_herits_df <- rbindlist(s_herits)
  
  sampled_h2[[j]] <- s_herits_df
}

sampled_h2_df <- rbindlist(sampled_h2)

ggplot(sampled_h2_df)+
  aes(x = pa_conc, y = H2, group = pa_conc, fill = factor(pa_conc))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Set1")+
  geom_jitter()+
  theme_bw()+
  # annotate("text", x = sampled_h2_df$pa_conc, y = sampled_h2_df$H2+herits_df$H2*.075, label = sampled_h2_df$H2)+
  theme(axis.text.x = element_text(size=16, face="bold", color="black"),
        axis.text.y = element_text(size=16, face="bold", color="black"),
        axis.title.x = element_text(size=24, face="bold", color="black", vjust=-.3),
        axis.title.y = element_text(size=24, face="bold", color="black"),
        strip.text.x = element_text(size=24, face="bold", color="black"),
        strip.text.y = element_text(size=16, face="bold", color="black"),
        plot.title = element_text(size=24, face="bold", vjust = 1),
        legend.position="none",
        panel.background = element_rect( color="black",size=1.2),
        strip.background = element_rect(color = "black", size = 1.2))+
  labs(x = "PA mM", y = "Broad Sense Heritability", title = "Sampled 3")

ggsave(filename = "~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/sampled_H2.pdf",
       height=6,
       width=8)