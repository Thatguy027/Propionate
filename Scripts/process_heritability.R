library(data.table)
library(tidyverse)
library(lme4)

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
setwd("..")

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


##################################### ##################################### Complete Dose response
raw_df <- fread("Data/Heritability/PA_doseresponse_herit.csv")

long_df <- raw_df %>%
  gather(strain_n, ratio, -DATE, -PA_conc) 

complete_dr <- data.frame(date = long_df$DATE, pa_conc = long_df$PA_conc, strain = long_df$strain_n, pheno = long_df$ratio)%>%
  separate(strain, into = c("strain", "replicate"), sep = "_")

# # # # # DR

# complete_drplt <- complete_dr %>%
#   dplyr::group_by(strain, pa_conc) %>%
#   dplyr::mutate(dose_mean = mean(pheno)) %>%
#   ggplot()+
#   aes(x = pa_conc, y = dose_mean, color = strain)+
#   geom_smooth(span=.8) +
#   theme_classic(20) +
#   labs(x = "Propionate Concentration (mM)", y = "L1 Survival", color = "Strain")

complete_dr_means <- complete_dr %>%
  dplyr::group_by(strain, pa_conc) %>%
  dplyr::mutate(dose_mean = mean(pheno))  

complete_drplt <- ggplot()+
  geom_smooth(aes(x = pa_conc, y = dose_mean, group = strain), color = "black", fill = "gray90",
              data = complete_dr_means %>% dplyr::filter(!(strain %in% c("N2", "DL238"))),
              span=.8) +
  geom_smooth(aes(x = pa_conc, y = dose_mean, color = strain), 
              data = complete_dr_means %>% dplyr::filter(strain %in% c("N2", "DL238", "BRC20067")),
              span=.8) +
  scale_x_continuous(limits = c(0, 150), breaks = unique(complete_dr_means$pa_conc)) +
  scale_y_continuous(limits = c(-10, 120), breaks = c(0,25,50,75,100)) +
  scale_color_manual(values = c("cadetblue3", "orange")) +
  geom_hline(aes(yintercept = 50), linetype = 2, color = "red", alpha = 0.7) +
  theme_classic(20) +
  labs(x = "Propionate Concentration (mM)", y = "L1 Survival (%)", color = "Strain")


complete_drplt

ggsave("Plots/Complete_Propionate_DR.pdf", 
       plot = complete_drplt,
       height = 4, width = 8)
ggsave("Plots/Complete_Propionate_DR.png", 
       plot = complete_drplt,
       height = 4, width = 8, dpi = 300)
ggsave("Plots/Complete_Propionate_DR.svg", 
       plot = complete_drplt,
       height = 4, width = 8)

# # # # # H2

herits <- list()
for(i in 1:length(unique(complete_dr$pa_conc))){
  for_calc <- dplyr::filter(complete_dr, pa_conc == unique(complete_dr$pa_conc)[i]) %>%
    dplyr::select(strain, pheno)
  
  herits[[i]] <- data.frame(H2 = signif(H2.fun(for_calc),2), pa_conc = unique(complete_dr$pa_conc)[i])
  
}

herits_df <- rbindlist(herits)

complete_h2 <- ggplot(herits_df)+
  aes(x = pa_conc, y = H2)+
  geom_point(size = 4)+
  theme_classic(20)+
  annotate("text", x = herits_df$pa_conc, y = herits_df$H2+.05, label = herits_df$H2)+
  scale_x_continuous(limits = c(0, 150), breaks = unique(complete_dr_means$pa_conc)) +
  labs(x = "Propionate Concentration (mM)", y = "Broad-Sense Heritability")

ggsave("Plots/Complete_Propionate_DR_H2.pdf", 
       plot = complete_h2,
       height = 4, width = 6)
ggsave("Plots/Complete_Propionate_DR_H2.png", 
       plot = complete_h2,
       height = 4, width = 6, dpi = 300)
ggsave("Plots/Complete_Propionate_DR_H2.svg", 
       plot = complete_h2,
       height = 4, width = 6)

##################################### ##################################### Fine-scale Dose response

pheno <- fread("Data/20150219_heritability_v2/cleaned_phenotypes.csv") %>%
  data.frame()%>%
  gather(strain, pheno, -Dose, -day)

pheno$strain <- gsub(pattern = "strain.", replacement = "", pheno$strain)
pheno$strain <- gsub(pattern = "\\.1$|\\.2$", replacement = "", pheno$strain)

# # # # # Fix strain names

pheno2 <- fread("Data/20150219_heritability_v2/strain_lookup.csv")
colnames(pheno2) <- c("strain", "strain_name")

pheno$strain <- as.numeric(pheno$strain)

pheno3 <- pheno %>%
  separate(Dose, into = c("Concentration", "value"), sep = " ",remove = FALSE)%>%
  mutate(n_conc = as.numeric(Concentration))%>%
  arrange(n_conc) %>%
  left_join(.,pheno2, by="strain")

# # # # # DR
fine_scale_dr_mean <- pheno3 %>%
  dplyr::group_by(strain_name, n_conc) %>%
  dplyr::mutate(dose_mean = mean(pheno, na.rm = T)) 

fine_scale_dr <- ggplot(fine_scale_dr_mean)+
  geom_smooth(aes(x = n_conc, y = dose_mean, group = strain_name), color = "black", fill = "gray80",
              data = fine_scale_dr_mean %>% dplyr::filter(!(strain %in% c("N2", "DL238"))),
              span=1) +
  geom_smooth(aes(x = n_conc, y = dose_mean, color = strain_name), 
              data = fine_scale_dr_mean %>% dplyr::filter(strain_name %in% c("N2", "DL238", "BRC20067")),
              span=1) +
  scale_color_manual(values = c("cadetblue3", "orange")) +
  scale_x_continuous(limits = c(70, 120), breaks = as.numeric(unique(fine_scale_dr_mean$Concentration))) +
  theme_classic(20) +
  geom_hline(aes(yintercept = 50), linetype = 2, color = "red", alpha = 0.7) +
  labs(x = "Propionate Concentration (mM)", y = "L1 Survival", color = "Strain")

ggsave("Plots/Finescale_Propionate_DR.pdf", 
       plot = fine_scale_dr, 
       height = 4, width = 8)
ggsave("Plots/Finescale_Propionate_DR.png", 
       plot = fine_scale_dr, 
       height = 4, width = 8, dpi = 300)
ggsave("Plots/Finescale_Propionate_DR.svg", 
       plot = fine_scale_dr, 
       height = 4, width = 8)


# # # # # H2 
H2.fun1 <- function(x){
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

clean_df <- pheno

clean_df$strain <- as.character(clean_df$strain)

sampled_h2 <- list()

for(j in 1:12){
  sampled_df <- clean_df%>%
    # dplyr::select(-date)%>%
    group_by(strain, Dose) %>%
    sample_n(3)
  
  s_herits <- list()
  for(i in 1:length(unique(sampled_df$Dose))){
    for_calc <- filter(sampled_df, Dose == unique(sampled_df$Dose)[i]) %>%
      dplyr::select(strain, pheno)
    
    s_herits[[i]] <- data.frame(H2 = signif(H2.fun1(for_calc),2), Dose = unique(sampled_df$Dose)[i], sample_n = j)
    
  }
  
  s_herits_df <- rbindlist(s_herits)
  
  sampled_h2[[j]] <- s_herits_df
}

sampled_h2_df <- rbindlist(sampled_h2)

sampled_h2_df <- sampled_h2_df %>% 
  separate(Dose, into = c("Concentration", "value"), sep = " ",remove = FALSE)%>%
  mutate(n_conc = as.numeric(Concentration))%>%
  arrange(n_conc)

finescale_h2 <- ggplot(sampled_h2_df)+
  aes(x = factor(n_conc,ordered = TRUE), y = H2, group = Dose)+
  geom_boxplot()+
  ggbeeswarm::geom_beeswarm()+
  theme_classic(20)+
  labs(x = "Propionate Concentration (mM)", y = "Broad-Sense Heritability")

ggsave("Plots/Finescale_Propionate_DR_H2.pdf",
       plot = finescale_h2,
       width = 8,
       height = 6)
ggsave("Plots/Finescale_Propionate_DR_H2.png",
       plot = finescale_h2,
       width = 8,
       height = 6,
       dpi = 300)

# # # # # Cowplots

# complete DR



cowplot::plot_grid(lemon::reposition_legend(complete_drplt, 'top right'),
                   complete_h2,
                   ncol =2, label_size = 20,
                   labels = "AUTO")

ggsave("Plots/Complete_DR_H2_Figure.pdf",
       width = 12,
       height = 6)
ggsave("Plots/Complete_DR_H2_Figure.png",
       width = 12,
       height = 6,
       dpi = 300)
ggsave("Plots/Complete_DR_H2_Figure.svg",
       width = 12,
       height = 6)

# finescale DR

cowplot::plot_grid(lemon::reposition_legend(fine_scale_dr, 'top right'),
                   finescale_h2,
                   ncol =2, label_size = 20,
                   labels = "AUTO")

ggsave("Plots/Finescale_DR_H2_Figure.pdf",
       width = 12,
       height = 6)
ggsave("Plots/Finescale_DR_H2_Figure.png",
       width = 12,
       height = 6,
       dpi = 300)
ggsave("Plots/Finescale_DR_H2_Figure.svg",
       width = 12,
       height = 6)
