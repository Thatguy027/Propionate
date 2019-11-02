library(pwr)
library(data.table)
library(tidyverse)
library(lme4)
library(ggbeeswarm)

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
setwd("..")

prop <- data.table::fread("Data/Power/20160920_cleaned_survival.csv")

ggplot(prop)+
  aes(x=strain,y= day2/day0)+
  geom_boxplot()+
  theme_classic(18)+
  geom_beeswarm()


dd = seq(0.01,1,.01)

samp.all <- list()

samps <- list()
for(i in 2:8){
  print(i)
  for(j in 1:10){
    print(j)

    for(d in 1:100){
      st.sd <- dplyr::filter(prop, strain =="DL238")%>%
        group_by(plate)%>%
        sample_n(i)%>%
        ungroup()%>%
        summarise(sds = sd(day2/day0))
      
      if(!exists("sds")){
        sds <- as.numeric(st.sd[1,1])
      } else {
        sds <- base::append(sds,as.numeric(st.sd[1,1]))
      }
      
    }
    
    pwrs <- list()
    for(e in 1:length(dd)){
      pwrs[[e]] = data.frame(p = power.t.test(n=nrow(dplyr::filter(prop, strain =="DL238")%>%
                                                       group_by(plate)%>%
                                                       sample_n(i)), 
                                              delta=dd[e],
                                              sd = mean(sds),
                                              sig.level=.00001)$power,
                             n = nrow(dplyr::filter(prop, strain =="DL238")%>%
                                        group_by(plate)%>%
                                        sample_n(i)),
                             d = dd[e],
                             sd = mean(sds))
    }
    samps[[j]] <- bind_rows(pwrs)%>%
      mutate(iteration = j)
    
    rm(sds)
  }
  samp.all[[i]] <- bind_rows(samps)
}

samp.all.df <- bind_rows(samp.all) %>%
  group_by(n,d)%>%
  mutate(m.sd = mean(p),
         sd.sd = sd(p))

write_tsv(samp.all.df, path = "Processed_Data/Power_samples.tsv",col_names = T)

ggplot(samp.all.df)+
  aes(x = d, y = m.sd, 
      color = factor(n),
      fill = factor(n))+
  geom_line() +
  theme_classic(18) +
  geom_hline(aes(yintercept = 0.8), color = "red")+
  geom_ribbon(aes(ymin=m.sd-sd.sd,ymax=m.sd+sd.sd), alpha= 0.05, linetype=4)+
  scale_color_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                     name = "Sample Size")+
  scale_fill_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                    name = "Sample Size")+
  xlim(0,.75) +
  labs(x = "Difference Between Means", y = "Power")

ggsave(filename = "Plots/Power_Calculation.pdf", height = 6, width = 10)
ggsave(filename = "Plots/Power_Calculation.png", height = 6, width = 10, dpi = 300)
ggsave(filename = "Plots/Power_Calculation.svg", height = 6, width = 10)
