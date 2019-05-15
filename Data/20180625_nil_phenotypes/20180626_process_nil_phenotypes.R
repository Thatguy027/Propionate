library(dplyr)
library(ggplot2)

setwd("~/Dropbox/AndersenLab/LabFolders/Stefan/Collaborations/propionic_acid/20180625_nil_phenotypes/")

prop_pheno <- data.table::fread("20180625_nil_phenos.csv")%>%
  na.omit()
prop_pheno$strain <- gsub("-","_",prop_pheno$strain)
prop_pheno$day <- as.factor(prop_pheno$day)

prop_zero <- prop_pheno %>%
  dplyr::group_by(strain,day)%>%
  dplyr::mutate(mean_surv = mean(survival))%>%
  dplyr::ungroup()%>%
  dplyr::distinct(strain, day, mean_surv)%>%
  dplyr::group_by(day)%>%
  dplyr::filter(mean_surv == max(mean_surv))%>%
  dplyr::select(max_surv = mean_surv, day)

prop_sub_max <- prop_pheno %>%
  dplyr::left_join(.,prop_zero,by="day")%>%
  dplyr::mutate(dif_surv = survival-max_surv)
  
ggplot(prop_sub_max)+
  aes(x = strain, y = dif_surv, color = day)+
  geom_point()

prop_sub_max <- prop_pheno %>%
  dplyr::left_join(.,prop_zero,by="day")%>%
  dplyr::mutate(dif_surv = survival/max_surv)%>%
  # dplyr::filter(strain != "DL238")%>%
  dplyr::mutate(resids = residuals(lm(dif_surv ~ day)))

ggplot(prop_sub_max)+
  aes(x = strain, y = resids, color = day)+
  geom_boxplot(outlier.colour = NA, color = "black")+
  geom_jitter(width = 0.2)+
  theme_bw(15)+
  labs(x = "Strain", y = "Normalized Survival")

ggsave("Normalzied_survival.pdf",height = 6, width = 10)

aov_res <- aov(prop_sub_max$resids ~ prop_sub_max$strain)
summary(aov_res)
tuk <- TukeyHSD(aov_res)

psig=as.numeric(apply(tuk$`prop_sub_max$strain`[,2:3],1,prod)>=0)+1
op=par(mar=c(4.2,9,3.8,2))
plot(tuk,col=psig,yaxt="n")
for (j in 1:length(psig)){
  axis(2,at=j,labels=rownames(tuk$`prop_sub_max$strain`)[length(psig)-j+1],
       las=1,cex.axis=.8,col.axis=psig[length(psig)-j+1])
}
par(op)

