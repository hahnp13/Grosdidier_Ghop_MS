library(tidyverse)
library(vegan) #Community Ecology
library(car) #Anova
library(DHARMa) #Model Testing: QQPlots, Predicted vs. Residuals, etc.
library(emmeans) #Least-Squared Means
library(glmmTMB) #Generalized linear mixed models
library(viridis) #Viridis Color Package
library(GGally) #GGpairs
library(performance) #Check for collinearity using vif for glmmTMB

#COMMUNITY SAMPLING ANALYSIS####
comm<-read.csv('Data/Grasshopper ID_nymphs ided_05202023.csv')%>% filter(HABITAT!='F')
fire<-read.csv('Data/Burn_by_Site Data_with frequency_05232023.csv') %>% filter(HABITAT!='F')
temp<-read.csv('Data/OSBS_Site Temp Data.csv') %>% filter(HABITAT!='F')


## acharum and apteno abundance ####
comm0check <-comm %>% 
  mutate(DATE=gsub('6-Jul-22','07/06/2022',DATE),
    DATE=as.Date(DATE, format = '%m/%d/%Y'),
         YEAR=substr(DATE, 1, 4),
    YEAR=paste('Y',YEAR,sep="")) %>% 
  #filter(DATE > as.Date('2022-01-01'), IDed=='Y') %>% 
  mutate(SUBFAMILY=ifelse(GENUS=='MELANOPLUS','MEL',SUBFAMILY),
         SPECIES=paste(substr(GENUS,1,2),substr(SPECIES,1,2),sep = "")) %>% 
  select(-c(HABITAT,GENUS,SEX,AGE,IDed,MORPH.CODE,SPP.NOTES)) %>% 
  group_by(YEAR,TIME.OF.YEAR,SITE,SPECIES) %>% 
  summarise(X..INDV=sum(X..INDV)) %>% 
  pivot_wider(names_from = c(YEAR,SPECIES), values_from = X..INDV) %>% 
  replace(is.na(.),0)

# write.csv(comm0check, "Outputs/Ordway_AcacarAptsph.csv")

ggplot(comm0check, aes(x=Y2021_ACCA, y=Y2022_ACCA)) +
  geom_point()+
  geom_smooth(method="glm", method.args=list(family=poisson))+
  theme_bw(base_size = 20)+
  facet_wrap(~TIME.OF.YEAR, scales="free")

ggplot(comm0check, aes(x=Y2021_ACCA, y=Y2022_ACCA)) +
  geom_point()+
  geom_smooth(method="glm", method.args=list(family=poisson))+
  theme_bw(base_size = 20)+
  facet_wrap(~TIME.OF.YEAR, scales="free")

ggplot(comm0check, aes(x=Y2021_APSP, y=Y2022_APSP)) +
  geom_point()+
  geom_smooth(method="glm", method.args=list(family=poisson))+
  theme_bw(base_size = 20)+
  facet_wrap(~TIME.OF.YEAR, scales="free")
