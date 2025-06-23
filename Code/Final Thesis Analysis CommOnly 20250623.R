library(tidyverse)
library(vegan) #Community Ecology
library(car) #Anova
library(DHARMa) #Model Testing: QQPlots, Predicted vs. Residuals, etc.
library(factoextra) #PCA
library(RVAideMemoire) #PermANOVA
library(emmeans) #Least-Squared Means
library(FD) #Functional Diversity, Use for trait matrix
library(glmmTMB) #Generalized linear mixed models
library(lme4) #Linear Mixed-Effects Models
library(lmerTest) #Tests in Linear Mixed-Effects Models
library(viridis) #Viridis Color Package
library(viridisLite) #Additional color-blind palettes
library(stringr)
library(Hmisc)
library(GGally) #GGpairs
library(ggrepel)
library(performance) #Check for collinearity using vif for glmmTMB
library(hillR) #Hill Numbers
library(cowplot) #Just used to get the legend from the feeding niche figures
library(MuMIn) #Extract R^2 values from models generated using glmmTMB, (r.squaredGLMM())
library(patchwork) #Piece plots together
library(picante)
library(ggeffects)

#GRASSHOPPER TRAIT PCA####

#Call in trait data
t1<-read.csv("Data/GHopp Trait Measurements_2022_withCN2.csv")

#Calculate FTs
#Here, locomotory traits are scaled by body size.
t2<-t1 %>% mutate(
  BV=(pi/4)*(BODY_LENGTH)*((DIAMETER_HEAD_WIDTH+DIAMETER_THORAX+DIAMETER_ABDOMEN_BASE+DIAMETER_ABDOMEN_TIP)/4)^2,
    IS=((HEAD_LENGTH*HEAD_HEIGHT*DIAMETER_HEAD_WIDTH)*(LA/LI)*(1/RI)),
  EYE_AREA=(EYE_WIDTH_DIAMETER*EYE_HEIGHT_DIAMETER),
  FEMUR_AREA=(FEMUR_LENGTH*FEMUR_WIDTH)/BV,
  TIBIA_LENGTH=TIBIA_LENGTH/BV,
  WINGSPAN=WINGSPAN/BV,
  CN_RATIO=(C/N))

#Remove extraneous columns and make species codes
t3<-t2 %>% select(-c(HABITAT,SITE,SAMPLE,LA,LI,RI,HEAD_LENGTH,HEAD_HEIGHT,BODY_LENGTH,DIAMETER_HEAD_WIDTH,
                     DIAMETER_THORAX,DIAMETER_ABDOMEN_BASE,DIAMETER_ABDOMEN_TIP,EYE_HEIGHT_DIAMETER,
                     EYE_WIDTH_DIAMETER,FEMUR_WIDTH,FEMUR_LENGTH,Comments,ID,SEX,C,N)) %>% 
  mutate(SPECIES=paste(substr(GENUS,1,2),substr(SPECIES,1,2),sep = "")) %>% 
  select(-c(GENUS)) %>% 
  mutate(SPECIES=case_when(SPECIES=='MEIM'~'MEKE',SPECIES!='MEIM'~SPECIES)) %>% 
  filter(SPECIES!='CHAU'&SPECIES!='PSFE')

#Look for outliers, covariates, and non normal distributions
ggpairs(t3)

#Transform non-normal data
t2<-t1 %>% mutate(
  BV=log((pi/4)*(BODY_LENGTH)*((DIAMETER_HEAD_WIDTH+DIAMETER_THORAX+DIAMETER_ABDOMEN_BASE+DIAMETER_ABDOMEN_TIP)/4)^2),
  IS=log(((HEAD_LENGTH*HEAD_HEIGHT*DIAMETER_HEAD_WIDTH)*(LA/LI)*(1/RI))),
  EYE_AREA=(EYE_WIDTH_DIAMETER*EYE_HEIGHT_DIAMETER),
  ANTENNAL_LENGTH=ANTENNAL_LENGTH,
  FEMUR_AREA=(FEMUR_LENGTH*FEMUR_WIDTH)/BV,
  TIBIA_LENGTH=TIBIA_LENGTH/BV,
  WINGSPAN=WINGSPAN/BV,
  CN_RATIO=(C/N))

t3<-t2 %>% mutate(SUBFAMILY=case_when(SPECIES=='ACCA' |
                                        SPECIES=='DIVI' |
                                        SPECIES=='ORPE' |
                                        SPECIES=='SYAD' |
                                        SPECIES=='EROB' |
                                        SPECIES=='AMMY' |
                                        SPECIES=='MEPI'~'Gomphocerinae',
                                      SPECIES=='APSP' |
                                        SPECIES=='SCAM' |
                                        SPECIES=='SCDA'~'Cyrtacanthacrididae',
                                      SPECIES=='MEKE' |
                                        SPECIES=='MERO'~'Melanoplinae',
                                      SPECIES=='PAPH' |
                                        SPECIES=='SPMA'~'Oedipodinae',
                                      SPECIES=='ODAP'~'Tettigoniidae')) 


t3<-t2 %>% select(-c(HABITAT,SITE,SAMPLE,LA,LI,RI,HEAD_LENGTH,HEAD_HEIGHT,BODY_LENGTH,DIAMETER_HEAD_WIDTH,
                     DIAMETER_THORAX,DIAMETER_ABDOMEN_BASE,DIAMETER_ABDOMEN_TIP,EYE_HEIGHT_DIAMETER,
                     EYE_WIDTH_DIAMETER,FEMUR_WIDTH,FEMUR_LENGTH,Comments,ID,SEX,C,N)) %>% 
  mutate(SPECIES=paste(substr(GENUS,1,2),substr(SPECIES,1,2),sep = "")) %>% 
  select(-c(GENUS)) %>% 
  mutate(SPECIES=case_when(SPECIES=='MEIM'~'MEKE',SPECIES!='MEIM'~SPECIES)) %>% 
  filter(SPECIES!='CHAU'&SPECIES!='PSFE',SPECIES!='ODAP',CN_RATIO<10)
ggpairs(t3)

t4<-t3 %>% mutate(SUBFAMILY=case_when(SPECIES=='ACCA' |
                                        SPECIES=='DIVI' |
                                        SPECIES=='ORPE' |
                                        SPECIES=='SYAD' |
                                        SPECIES=='EROB' |
                                        SPECIES=='AMMY' |
                                        SPECIES=='MEPI'~'Gomphocerinae',
                                      SPECIES=='APSP' |
                                        SPECIES=='SCAM' |
                                        SPECIES=='SCDA'~'Cyrtacanthacrididae',
                                      SPECIES=='MEKE' |
                                        SPECIES=='MERO'~'Melanoplinae',
                                      SPECIES=='PAPH' |
                                        SPECIES=='SPMA'~'Oedipodinae')) 
                                     #SPECIES=='ODAP'~'Tettigoniidae')) 

ggplot(t4 %>% arrange(SUBFAMILY), aes(x=SUBFAMILY, y=BV, color=SPECIES))+
  geom_boxplot(outlier.shape = NA)+
  #geom_jitter(width=.2, height=0)+
  scale_color_viridis(discrete = T, direction = -1, end=.8)+
  theme_bw(base_size = 16)

ggplot(t4 %>% arrange(SUBFAMILY), aes(x=SUBFAMILY, y=IS, color=SPECIES))+
  geom_boxplot(outlier.shape = NA)+
  #geom_jitter(width=.2, height=0)+
  scale_color_viridis(discrete = T, direction = -1, end=.8)+
  theme_bw(base_size = 16)

ggplot(t4 %>% arrange(SUBFAMILY), aes(x=SUBFAMILY, y=CN_RATIO, color=SPECIES))+
  geom_boxplot()+
  #geom_jitter(width=.2, height=0)+
  scale_color_viridis(discrete = T, direction = -1, end=.8)+
  theme_bw(base_size = 16)

ggpairs(t4[c(5,6,9)])+theme_bw(base_size = 16)

#Right skewed data, specifically IS and BV were log transformed. Distribution now follows Gaussian trend.

## Average traits across species ####
t5<-t4 %>% group_by(SPECIES,SUBFAMILY) %>%
  summarise(BV=mean(BV,na.rm=T),IS=mean(IS,na.rm=T),CN_RATIO=mean(CN_RATIO,na.rm=T),
            WINGSPAN=mean(WINGSPAN,na.rm=T),FEMUR_AREA=mean(FEMUR_AREA,na.rm=T),
            TIBIA_LENGTH=mean(TIBIA_LENGTH,na.rm=T),ANTENNAL_LENGTH=mean(ANTENNAL_LENGTH,na.rm=T),
            EYE_AREA=mean(EYE_AREA,na.rm=T))

## Principle Component Analysis ####
# t5<-t4 %>% mutate(SUBFAMILY=case_when(SPECIES=='ACCA' |
#                                SPECIES=='DIVI' |
#                                SPECIES=='ORPE' |
#                                SPECIES=='SYAD' |
#                                SPECIES=='EROB' |
#                                SPECIES=='AMMY' |
#                                  SPECIES=='MEPI'~'Gomphocerinae',
#                              SPECIES=='APSP' |
#                                SPECIES=='SCAM' |
#                                SPECIES=='SCDA'~'Cyrtacanthacrididae',
#                              SPECIES=='MEKE' |
#                                SPECIES=='MERO'~'Melanoplinae',
#                              SPECIES=='PAPH' |
#                                SPECIES=='SPMA'~'Oedipodinae',
#                              SPECIES=='ODAP'~'Tettigoniidae')) 

t4a <- t4[c(1,5,6,9,10)] %>% na.omit()

t.pca <- prcomp(t4a[2:4] %>% na.omit(), scale = TRUE)
t.pca$rotation
summary(t.pca)
t4a$pc1 <- t.pca$x[,1]
t4a$pc2 <- t.pca$x[,2]

t5a<-t4a %>% group_by(SPECIES,SUBFAMILY) %>%
  summarise(BV=mean(BV,na.rm=T),IS=mean(IS,na.rm=T),CN_RATIO=mean(CN_RATIO,na.rm=T),
            pc1=mean(pc1,na.rm=T),pc2=mean(pc2,na.rm=T))

t.pca1 <- prcomp(t5[3:5], scale = TRUE)
t.pca1$rotation
summary(t.pca1)
biplot(t.pca1)

ggpairs(t5a)


# t.pca_alltraits<-fviz_pca_biplot(t.pca, label = "var", pointsize = 4) +
#   geom_point(aes(x = t.pca$x[, 1], y = t.pca$x[, 2], color = t4a$SUBFAMILY), size=5) +
#   #geom_text_repel(aes(x = t.pca$x[,1], y = t.pca$x[,2], label = t4a$SPECIES),size = 4, force = 10) +
#   scale_color_viridis(discrete = T, option = 'C', end = .875) +
#   xlab('PCA1 (85.6%)') + ylab('PCA2 (8.9%)') +  ggtitle('') + labs(color='Subfamily') +
#   theme_bw(base_size = 24) +
#   theme(panel.grid = element_blank()) ;t.pca_alltraits

t.pca_alltraits<-fviz_pca_biplot(t.pca1, label = "var", pointsize = 4) +
  geom_point(aes(x = t.pca1$x[, 1], y = t.pca1$x[, 2], color = t5$SUBFAMILY), size=5) +
  geom_text_repel(aes(x = t.pca1$x[,1], y = t.pca1$x[,2], label = t5$SPECIES),size = 4, force = 10) +
  scale_color_viridis(discrete = T, option = 'C', end = .875) +
  xlab('PCA1 (85.6%)') + ylab('PCA2 (8.9%)') +  ggtitle('') + labs(color='Subfamily') +
  theme_bw(base_size = 24) +
  theme(panel.grid = element_blank()) ;t.pca_alltraits

t.pca_alltraits2<-ggplot() +
  geom_point(data=t4a, aes(x=pc1, y=pc2, color=SUBFAMILY), pch=16, size=2, alpha=.5) +
  stat_ellipse(data=t4a, aes(x=pc1, y=pc2, color=SUBFAMILY), level=.66, lwd=2)+
  geom_point(data=t5a, aes(x=pc1, y=pc2, color=SUBFAMILY), size=5) +
  #geom_text_repel(aes(x = t.pca$x[,1], y = t.pca$x[,2], label = t4a$SPECIES),size = 4, force = 10) +
  scale_color_viridis(discrete = T, option = 'C', end = .875) +
  scale_fill_viridis(discrete = T, option = 'C', end = .875) +
  xlab('PCA1 (85.6%)') + ylab('PCA2 (8.9%)') +  
  ggtitle('') + labs(color='Subfamily') +
  theme_bw(base_size = 24) +
  theme(panel.grid = element_blank());t.pca_alltraits2


ggsave("t.pca_alltraits2.tiff", t.pca_alltraits, width=10, height=6, units="in", dpi=600, compression = "lzw", path="Outputs")
ggsave("t.pca_reduced2.tiff", t.pca_reduced, width=13, height=10, units="in", dpi=600, compression = "lzw", path="Outputs")

## MANOVA on functional traits ####
tniche <- manova(cbind(BV, IS,CN_RATIO) ~ SUBFAMILY+SPECIES, data=t4a)
summary(tniche)
summary.aov(tniche)

t_bv0 <- glmmTMB(BV ~ 1 , data=t4)
t_bv1 <- glmmTMB(BV ~ 1 + (1|SPECIES), data=t4)
t_bv2 <- glmmTMB(BV ~ 1 + (1|SUBFAMILY/SPECIES), data=t4)
anova(t_bv0,t_bv1,t_bv2)

t_is0 <- glmmTMB(IS ~ 1 , data=t4)
t_is1 <- glmmTMB(IS ~ 1 + (1|SPECIES), data=t4)
t_is2 <- glmmTMB(IS ~ 1 + (1|SUBFAMILY), data=t4)
t_is3 <- glmmTMB(IS ~ 1 + (1|SUBFAMILY/SPECIES), data=t4)
anova(t_is0,t_is1,t_is2,t_is3)
summary(t_is3)

t_bv0 <- glmmTMB(CN_RATIO ~ 1 , data=t4)
t_bv1 <- glmmTMB(CN_RATIO ~ 1 + (1|SPECIES), data=t4)
t_bv2 <- glmmTMB(CN_RATIO ~ 1 + (1|SUBFAMILY/SPECIES), data=t4)
anova(t_bv0,t_bv1,t_bv2)


t_bv <- glmmTMB(CN_RATIO ~ SUBFAMILY+SPECIES, data=t4)
Anova(t_bv)
summary(t_bv)


#COMMUNITY SAMPLING ANALYSIS####
comm<-read.csv('Data/Grasshopper ID_nymphs ided_05202023.csv')
fire<-read.csv('Data/Burn_by_Site Data_with frequency_05232023.csv')
temp<-read.csv('Data/OSBS_Site Temp Data.csv')

## acharum and apteno abundance ####
comm0check <-comm %>% 
  mutate(DATE=as.Date(DATE, format = '%m/%d/%Y')) %>% 
  filter(DATE > as.Date('2022-01-01'), IDed=='Y') %>% 
  mutate(SUBFAMILY=ifelse(GENUS=='MELANOPLUS','MEL',SUBFAMILY),
         SPECIES=paste(substr(GENUS,1,2),substr(SPECIES,1,2),sep = "")) %>% 
  dplyr::select(-c(DATE,HABITAT,GENUS,SEX,AGE,IDed,MORPH.CODE,SPP.NOTES)) %>% 
  group_by(TIME.OF.YEAR,SITE,SPECIES) %>% 
  summarise(X..INDV=sum(X..INDV)) %>% 
  pivot_wider(names_from = SPECIES, values_from = X..INDV) %>% 
  replace(is.na(.),0) %>%
  filter(TIME.OF.YEAR=='EARLY') %>% 
  select(TIME.OF.YEAR, SITE, ACCA, APSP)
write.csv(comm0check, "Outputs/Ordway_AcacarAptsph.csv")

#FUNCTIONAL TRAIT ANALYSIS ON COMMUNITY DATA
comm1<-comm %>% 
  mutate(DATE=as.Date(DATE, format = '%m/%d/%Y')) %>% 
  filter(DATE > as.Date('2022-01-01'), IDed=='Y') %>% 
  mutate(SUBFAMILY=ifelse(GENUS=='MELANOPLUS','MEL',SUBFAMILY),
         SPECIES=paste(substr(GENUS,1,2),substr(SPECIES,1,2),sep = "")) %>% 
  select(-c(DATE,HABITAT,GENUS,SEX,AGE,IDed,MORPH.CODE,SPP.NOTES)) %>% 
  group_by(TIME.OF.YEAR,SITE,SPECIES) %>% 
  summarise(X..INDV=sum(X..INDV)) %>% 
  pivot_wider(names_from = SPECIES, values_from = X..INDV) %>% 
  replace(is.na(.),0)

t7<-t2 %>% select(-c(HABITAT,SITE,SAMPLE,LA,LI,RI,HEAD_LENGTH,HEAD_HEIGHT,BODY_LENGTH,DIAMETER_HEAD_WIDTH,
                     DIAMETER_THORAX,DIAMETER_ABDOMEN_BASE,DIAMETER_ABDOMEN_TIP,EYE_HEIGHT_DIAMETER,
                     EYE_WIDTH_DIAMETER,FEMUR_WIDTH,FEMUR_LENGTH,Comments,ID,SEX,C,N)) %>% 
  mutate(SPECIES=paste(substr(GENUS,1,2),substr(SPECIES,1,2),sep = "")) %>% 
  select(-c(GENUS)) %>% filter(SPECIES!='CHAU'&SPECIES!='PSFE') %>% 
  group_by(SPECIES) %>% 
  summarise(IS=mean(IS,na.rm=T),
            CN_RATIO=mean(CN_RATIO,na.rm=T),
            WINGSPAN=mean(WINGSPAN,na.rm=T),
            ANTENNAL_LENGTH=mean(ANTENNAL_LENGTH,na.rm=T)) %>% 
  mutate(SUBFAMILY=case_when(SPECIES=='ACCA' |
                               SPECIES=='DIVI' |
                               SPECIES=='ORPE' |
                               SPECIES=='SYAD' |
                               SPECIES=='EROB' |
                               SPECIES=='AMMY' |
                               SPECIES=='MEPI'~'Gomphocerinae',
                             SPECIES=='APSP' |
                               SPECIES=='SCAM' |
                               SPECIES=='SCDA'~'Cyrtacanthacrididae',
                             SPECIES=='MEKE' |
                               SPECIES=='MERO' |
                               SPECIES=='MEIM'~'Melanoplinae',
                             SPECIES=='PAPH' |
                               SPECIES=='SPMA'~'Oedipodinae',
                             SPECIES=='ODAP'~'Tettigoniidae'),
         FEEDING_GUILD=case_when(SUBFAMILY=='Gomphocerinae'~0,
                                 SUBFAMILY=='Cyrtacanthacrididae' |
                                   SUBFAMILY=='Oedipodinae' |
                                   SUBFAMILY=='Tettigoniidae'~1,
                                 SUBFAMILY=='Melanoplinae'~2)) %>% 
  select(-SUBFAMILY) %>% 
  remove_rownames %>% column_to_rownames(var="SPECIES") %>% as.matrix()

comm2<-comm1 %>% select(-c(BESU,ARXA,ARNA,NANA,APAP,HIOC,SCFU)) %>% 
  .[3:18] %>% select(order(colnames(.))) %>% as.matrix()

functraits2 <- functcomp(t7, comm2)
comm3<-cbind(comm1,functraits2) %>% rename(TIME_OF_YEAR = TIME.OF.YEAR)

#ADD IN FIRE DATA
fire1<-fire %>% filter(YEAR==2022) %>% select(-c(YEAR,HABITAT))
plot(fire1$TIME_SINCE_FIRE_DAYS,fire1$FIRE_FREQUENCY) #Log transform TIME_SINCE_FIRE_DAYS in the joined dataset
comm4<-full_join(comm3,fire1,by = c('TIME_OF_YEAR','SITE')) %>% 
  mutate(TIME_SINCE_FIRE_DAYS=log(TIME_SINCE_FIRE_DAYS))

head(comm4)

ggplot(comm4, aes(x=TIME_SINCE_FIRE_DAYS, y=SCAM)) +
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~TIME_OF_YEAR)


#ADD IN TEMP DATA
temp[320,1]="9/17/2022"
temp[320,2]="14:00"

temp1<-temp %>% 
  mutate(DATE=as.Date(DATE, format = '%m/%d/%Y')) %>% filter(DATE > as.Date('2022-01-01')) %>% 
  mutate(TIME_OF_YEAR=case_when(DATE > '2022-08-01'~'LATE',DATE < '2022-08-01'~'EARLY')) %>% 
  mutate(DATE)%>% 
  group_by(TIME_OF_YEAR,SITE,DATE,TIME,SAMPLE) %>% 
  pivot_longer(TEMP.1:TEMP.4,names_to = 'QUADRAT', values_to = 'TEMP') %>% 
  group_by(TIME_OF_YEAR,SITE,DATE,TIME) %>% 
  summarise(TEMP_M=mean(TEMP, na.rm=T), TEMP_SD=sd(TEMP, na.rm=T)) %>% 
  mutate(TEMP_CV=(TEMP_SD/TEMP_M)*100) %>% 
  select(-c(TEMP_SD))
plot(temp1$TEMP_M,temp1$TEMP_CV) #No need to log transform variables
comm5<-full_join(comm4,temp1,by = c('TIME_OF_YEAR','SITE'))

#CALCULATE DIVERSITY METRICES USING HILL NUMBERS
comm5$richness<-hill_taxa(comm5[3:25], q=0)
comm5$shannon<-hill_taxa(comm5[3:25], q=1)
comm5$simpsonsinv<-hill_taxa(comm5[3:25], q=2)

## COMMUNITY SES ####
t7a<-t2 %>% select(-c(HABITAT,SITE,SAMPLE,LA,LI,RI,HEAD_LENGTH,HEAD_HEIGHT,BODY_LENGTH,DIAMETER_HEAD_WIDTH,
                      DIAMETER_THORAX,DIAMETER_ABDOMEN_BASE,DIAMETER_ABDOMEN_TIP,EYE_HEIGHT_DIAMETER,
                      EYE_WIDTH_DIAMETER,FEMUR_WIDTH,FEMUR_LENGTH,Comments,ID,SEX,C,N)) %>% 
  mutate(SPECIES=paste(substr(GENUS,1,2),substr(SPECIES,1,2),sep = "")) %>% 
  select(-c(GENUS)) %>% filter(SPECIES!='CHAU'&SPECIES!='PSFE') %>% 
  group_by(SPECIES) %>% 
  summarise(BV=mean(BV, na.rm=T),
            IS=mean(IS,na.rm=T),
            CN_RATIO=mean(CN_RATIO,na.rm=T),
            WINGSPAN=mean(WINGSPAN,na.rm=T),
            ANTENNAL_LENGTH=mean(ANTENNAL_LENGTH,na.rm=T)) %>% 
  mutate(SUBFAMILY=case_when(SPECIES=='ACCA' |
                               SPECIES=='DIVI' |
                               SPECIES=='ORPE' |
                               SPECIES=='SYAD' |
                               SPECIES=='EROB' |
                               SPECIES=='AMMY' |
                               SPECIES=='MEPI'~'Gomphocerinae',
                             SPECIES=='APSP' |
                               SPECIES=='SCAM' |
                               SPECIES=='SCDA'~'Cyrtacanthacrididae',
                             SPECIES=='MEKE' |
                               SPECIES=='MERO' |
                               SPECIES=='MEIM'~'Melanoplinae',
                             SPECIES=='PAPH' |
                               SPECIES=='SPMA'~'Oedipodinae',
                             SPECIES=='ODAP'~'Tettigoniidae'),
         FEEDING_GUILD=case_when(SUBFAMILY=='Gomphocerinae'~0,
                                 SUBFAMILY=='Cyrtacanthacrididae' |
                                   SUBFAMILY=='Oedipodinae' |
                                   SUBFAMILY=='Tettigoniidae'~1,
                                 SUBFAMILY=='Melanoplinae'~2)) %>% 
  remove_rownames %>% column_to_rownames(var="SPECIES")
t8 <- t5[c(1:5)] %>% remove_rownames %>% column_to_rownames(var="SPECIES");

comm3 <- comm2 %>% as.data.frame() %>% 
  select(ACCA,AMMY,APSP,DIVI,EROB,MEKE,MEPI,MERO,ORPE,PAPH,SCAM,SCDA,SPMA,SYAD)

comm_roa <- comm4 %>% select(TIME_OF_YEAR,SITE,TIME_SINCE_FIRE_DAYS,FIRE_FREQUENCY)

## Functional dispersion ####
### Rao's randomization -- All traits ####

frao_All <- hill_func(comm3, t8)[1,]
frao_Allr <- data.frame(data.frame(matrix(ncol = 9, nrow = 60)))

for (i in 1:999) {#For each row in the matrix (for each site)
  #select randomly, as many species as the species richness of the site:
  comm3rand <- comm3
  #colnames(comm3rand) <- sample(colnames(comm3), replace=F)
  comm3rand <- comm3rand %>% select(order(colnames(.)))
   
  t8rand <- t8
  t8rownames <- rownames(t8)
  t8rand <- t8rand[sample(nrow(t8rand)),]
  rownames(t8rand) <- t8rownames
  #colnames(t8rand) <- sample(colnames(t8), replace=F) 
  #t8rand <- t8rand %>% select(order(colnames(.)))
  
  frao_Allr[,i] <- hill_func(comm3rand, t8rand)[1,]
}

All_SEStab <- frao_Allr %>% 
  cbind(.,frao_All) %>%
  rowwise() %>% 
  mutate(rand_mean = mean(c_across(2:1000)),rand_sd = sd(c_across(2:1000))) %>% 
  select(frao_All,rand_mean,rand_sd) %>% 
  mutate(SES_All=(frao_All-rand_mean)/rand_sd, pval=pnorm(q=SES_All, lower.tail = T)) 

summary(All_SEStab)  
hist(All_SEStab$SES_All)
comm_roa$SES_All <- All_SEStab$SES_All

ggplot(comm_roa, aes(x=FIRE_FREQUENCY, y=SES_All))+
  geom_point()+
  geom_smooth(method="lm", formula=y~x+I(x^2))

#mutate(rand_mean = mean(c_across(2:1000)),rand_sd = sd(c_across(2:1000)), rank=(ecdf(c_across(2:1000))(.$frao_IS[1]))) %>% 

### Rao's randomization -- BV trait ####

frao_BV <- hill_func(comm3, t8[2])[1,]
frao_BVr <- data.frame(data.frame(matrix(ncol = 9, nrow = 60)))

for (i in 1:999) {#For each row in the matrix (for each site)
  #select randomly, as many species as the species richness of the site:
  comm3rand <- comm3
  #colnames(comm3rand) <- sample(colnames(comm3), replace=F) 
  comm3rand <- comm3rand %>% select(order(colnames(.)))
  
  t8rand <- t8
  t8rownames <- rownames(t8)
  t8rand <- t8rand[sample(nrow(t8rand)),]
  rownames(t8rand) <- t8rownames
  
  frao_BVr[,i] <- hill_func(comm3rand, t8rand[2])[1,]
}

BV_SEStab <- frao_BVr %>% 
  cbind(.,frao_BV) %>%
  rowwise() %>% 
  mutate(rand_mean = mean(c_across(2:1000)),rand_sd = sd(c_across(2:1000))) %>% 
  select(frao_BV,rand_mean,rand_sd) %>% 
  mutate(SES_BV=(frao_BV-rand_mean)/rand_sd, pval=pnorm(q=SES_BV, lower.tail = T)) 

summary(BV_SEStab)  
hist(BV_SEStab$SES_BV)
comm_roa$SES_BV <- BV_SEStab$SES_BV

### Rao's randomization -- IS trait ####

frao_IS <- hill_func(comm3, t8[3])[1,]
frao_ISr <- data.frame(data.frame(matrix(ncol = 9, nrow = 60)))

for (i in 1:999) {#For each row in the matrix (for each site)
  #select randomly, as many species as the species richness of the site:
  comm3rand <- comm3
  #colnames(comm3rand) <- sample(colnames(comm3), replace=F) 
  comm3rand <- comm3rand %>% select(order(colnames(.)))
  
  t8rand <- t8
  t8rownames <- rownames(t8)
  t8rand <- t8rand[sample(nrow(t8rand)),]
  rownames(t8rand) <- t8rownames
  
  frao_ISr[,i] <- hill_func(comm3rand, t8rand[3])[1,]
}

IS_SEStab <- frao_ISr %>% 
  cbind(.,frao_IS) %>%
  rowwise() %>% 
  mutate(rand_mean = mean(c_across(2:1000)),rand_sd = sd(c_across(2:1000)), rank=(ecdf(c_across(2:1000))(.$frao_IS[1]))) %>% 
  select(frao_IS,rand_mean,rand_sd) %>% 
  mutate(SES_IS=(frao_IS-rand_mean)/rand_sd, pval=pnorm(q=SES_IS, lower.tail = T)) 

summary(IS_SEStab)  
hist(IS_SEStab$SES_IS)
comm_roa$SES_IS <- IS_SEStab$SES_IS

### Rao's randomization -- CN trait ####

frao_CN <- hill_func(comm3, t8[4])[1,]
frao_CNr <- data.frame(data.frame(matrix(ncol = 9, nrow = 60)))

for (i in 1:999) {#For each row in the matrix (for each site)
  #select randomly, as many species as the species richness of the site:
  comm3rand <- comm3
  #colnames(comm3rand) <- sample(colnames(comm3), replace=F) 
  comm3rand <- comm3rand %>% select(order(colnames(.)))
  
  t8rand <- t8
  t8rownames <- rownames(t8)
  t8rand <- t8rand[sample(nrow(t8rand)),]
  rownames(t8rand) <- t8rownames
  
  frao_CNr[,i] <- hill_func(comm3rand, t8rand[4])[1,]
}

CN_SEStab <- frao_CNr %>% 
  cbind(.,frao_CN) %>%
  rowwise() %>% 
  mutate(rand_mean = mean(c_across(2:1000)),rand_sd = sd(c_across(2:1000))) %>% 
  select(frao_CN,rand_mean,rand_sd) %>% 
  mutate(SES_CN=(frao_CN-rand_mean)/rand_sd, pval=pnorm(q=SES_CN, lower.tail = T)) 

summary(CN_SEStab)  
hist(CN_SEStab$SES_CN)
comm_roa$SES_CN <- CN_SEStab$SES_CN

### Rao's randomization -- SF trait ####

frao_SF <- hill_func(comm3, t8[1])[1,]
frao_SFr <- data.frame(data.frame(matrix(ncol = 9, nrow = 60)))

for (i in 1:999) {#For each row in the matrix (for each site)
  #select randomly, as many species as the species richness of the site:
  comm3rand <- comm3
  #colnames(comm3rand) <- sample(colnames(comm3), replace=F) 
  comm3rand <- comm3rand %>% select(order(colnames(.)))
  
  t8rand <- t8
  t8rownames <- rownames(t8)
  t8rand <- t8rand[sample(nrow(t8rand)),]
  rownames(t8rand) <- t8rownames
  
  frao_SFr[,i] <- hill_func(comm3rand, t8rand[1])[1,]
}

SF_SEStab <- frao_SFr %>% 
  cbind(.,frao_SF) %>%
  rowwise() %>% 
  mutate(rand_mean = mean(c_across(2:1000)),rand_sd = sd(c_across(2:1000))) %>% 
  select(frao_SF,rand_mean,rand_sd) %>% 
  mutate(SES_SF=(frao_SF-rand_mean)/rand_sd, pval=pnorm(q=SES_SF, lower.tail = T)) 

summary(SF_SEStab)  
hist(SF_SEStab$SES_SF)
comm_roa$SES_SF <- SF_SEStab$SES_SF

## Rao figures ####
comm_roa2 <- comm_roa %>% group_by(TIME_OF_YEAR,SITE,TIME_SINCE_FIRE_DAYS,FIRE_FREQUENCY) %>% 
  pivot_longer(cols=5:9, values_to = "SES", names_to = "Trait")


trait_names <- c(`SES_All` = "All traits", `SES_BV`="Body Volume", `SES_CN`="C:N", `SES_IS`="Incisor strength",
                 `SES_SF`="Subfamily")

raoSESfig <- sesTrait <- ggplot(comm_roa2, aes(x=TIME_OF_YEAR, y=SES, color=TIME_OF_YEAR)) + 
  geom_boxplot(outlier.shape=NA, lwd=1)+
  geom_jitter(height=0, width=.2) + 
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_color_viridis(option="G", discrete=T, direction=-1, end=.6, name = "Time of year", labels=c("July","Sept."))+
  scale_y_continuous(name="Standardized Effect Size \n Rao's Func. Disp.")+
  #scale_x_discrete(name = "Time of year", labels=c("July","Sept."))+
  facet_wrap(~Trait, labeller=as_labeller(trait_names))+
  theme_bw(base_size = 16)+ 
  theme(legend.position = c(0.85, 0.2), # c(0,0) bottom left, c(1,1) top-right.
        legend.background = element_rect(fill = "white", colour = NA))

raoSESFirefig <- sesTrait <- ggplot(comm_roa2, aes(x=TIME_SINCE_FIRE_DAYS, y=SES, color=TIME_OF_YEAR)) + 
  geom_point() +
  geom_smooth(method="lm")+
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_color_viridis(option="G", discrete=T, direction=-1, end=.6, name = "Time of year", labels=c("July","Sept."))+
  scale_y_continuous(name="Standardized Effect Size \n Rao's Func. Disp.")+
  #scale_x_discrete(name = "Time of year", labels=c("July","Sept."))+
  facet_wrap(~Trait, labeller=as_labeller(trait_names))+
  theme_bw(base_size = 16)+ 
  theme(legend.position = c(0.85, 0.2), # c(0,0) bottom left, c(1,1) top-right.
        legend.background = element_rect(fill = "white", colour = NA))

ggsave("raoSESfig_revised.tiff", raoSESfig, width=8, height=6, units="in", dpi=600, compression = "lzw", path="Outputs")


glmrao0 <- glmmTMB(SES ~ Trait*TIME_OF_YEAR + (1|SITE/Trait), data=comm_roa2, dispformula = ~Trait)
#glmrao <- glmmTMB(SES ~ Trait*TIME_OF_YEAR + (1|SITE/TIME_OF_YEAR), data=comm_roa2)
#anova(glmrao0,glmrao)
Anova(glmrao0)
summary(glmrao0)

emmeans(glmrao0, ~Trait|TIME_OF_YEAR, infer=T)
emmeans(glmrao0, ~Trait, infer=T)

#FEEDING NICHE AND TRAIT MATCHING########################################################################
#LOAD IN PLANT DATA AND PROCESS PLANT DATA
sla<-read.csv('Data/Leaf SLA_2022.csv') %>% select(-c(DRY_WT,AREA))
ldmc<-read.csv('Data/Plant CNP.csv') %>% select(-c(WET_WT,DRY_WT))
cndata<-read.csv('Data/Hahn Plant TC TN results march 2023.csv')

cndata2<-cndata %>% mutate(TUBE.ID=as.numeric(gsub("^(\\d+).*", "\\1", cndata$Sample.ID)),
                           CN_RATIO=wt..C/wt..N) %>% 
  select(-c(Date,X,X.1,wt..C,wt..N))

cnldmc<-full_join(cndata2,ldmc,by = "TUBE.ID")
ptrait<-full_join(cnldmc,sla,by = c('PLANT_ID','SPECIES')) %>% 
  select(-Sample.ID) %>% select(TUBE.ID,PLANT_ID,SPECIES,everything()) %>% rename(TUBE_ID=TUBE.ID)
ptrait2<-ptrait %>% group_by(SPECIES) %>% 
  summarise(CN_RATIO=mean(CN_RATIO, na.rm=T), LDMC=mean(LDMC,na.rm=T), SLA=mean(SLA,na.rm=T))

#LOOK AT PLANT TRAIT GRADIENTS AND CORRELATIONS
ggplot(data = ptrait)+geom_boxplot(aes(x=reorder(SPECIES,CN_RATIO), y=CN_RATIO))
ggplot(data = ptrait)+geom_boxplot(aes(x=reorder(SPECIES,LDMC), y=LDMC))+ylim(200,600)
ggplot(data = ptrait)+geom_boxplot(aes(x=reorder(SPECIES,SLA), y=SLA))
plot(SLA~CN_RATIO,data = ptrait2)
abline(lm(SLA~CN_RATIO,data = ptrait2))
plot(SLA~LDMC,data = ptrait2)
abline(lm(SLA~LDMC,data = ptrait2))
plot(LDMC~CN_RATIO,data = ptrait2)
abline(lm(LDMC~CN_RATIO,data = ptrait2))

ggpairs(ptrait[4:6] %>% na.omit())
ggpairs(ptrait2[2:4])

# LOAD IN FEEDING TRIAL DATA ####
exp<-read.csv('Data/GHopp Feeding Trials_2022.csv')

## check ACCA and APSP data ####
aadat <- exp %>% filter(GRASSHOPPER_SPECIES == 'ACCA' | GRASSHOPPER_SPECIES == 'APSP'  ) %>% 
  pivot_longer(cols=ARBE:VAMY, names_to = "PlantSpp", values_to = "Amt_Eaten")

aadatsum <- aadat %>% group_by(GRASSHOPPER_SPECIES, PlantSpp) %>% summarise(Amt_Eaten_m = mean(Amt_Eaten, na.rm=T))
## Microniche Breadth ####
nmat1 <- exp %>% filter(GRASSHOPPER_SPECIES == 'ACCA' | GRASSHOPPER_SPECIES == 'APSP'| GRASSHOPPER_SPECIES == 'DIVI'  ) 

nmat <- nmat1 %>% 
  dplyr::select(GRASSHOPPER_SPECIES, ARBE:VAMY) %>% 
  rowwise() %>% 
  mutate(Consumed=sum(c_across(2:17)))

nmat2 <- (nmat[2:17]/nmat$Consumed)
nmat3 <- nmat2 %>% rowwise() %>% 
  mutate(B=1/(sum(c_across(everything()))), Ba=(B-1)/(16-1))
nmat$B <- nmat3$B
nmat$Ba <- nmat3$Ba

library(janitor)
nmat4 <- t(nmat[1:16]) %>% row_to_names(row_number = 1)



#CALCULATE CWM PLANT FUNCTIONAL TRAIT IN GHOP DIET
ptrait3<-ptrait2 %>% remove_rownames %>% column_to_rownames(var="SPECIES") %>% as.matrix()
exp2 <- exp[6:21] %>% select(order(colnames(.))) %>% as.matrix()
functraits <- functcomp(ptrait3, exp2)

exp3 <- cbind(exp,functraits) %>% 
  mutate(SUBFAMILY=case_when(GRASSHOPPER_SPECIES=='ACCA' |
                               GRASSHOPPER_SPECIES=='DIVI' |
                               GRASSHOPPER_SPECIES=='ORPE' |
                               GRASSHOPPER_SPECIES=='SYAD'~'Gomphocerinae',
                             GRASSHOPPER_SPECIES=='APSP' |
                               GRASSHOPPER_SPECIES=='SCAM' |
                               GRASSHOPPER_SPECIES=='SCDA'~'Cyrtacanthacrididae',
                             GRASSHOPPER_SPECIES=='MEKE' |
                               GRASSHOPPER_SPECIES=='MERO'~'Melanoplinae',
                             GRASSHOPPER_SPECIES=='PAPH' |
                               GRASSHOPPER_SPECIES=='SPMA'~'Oedipodinae'))

## Microniche Breadth ####
nmat <- exp3 %>% filter(GRASSHOPPER_SPECIES == 'ACCA' | GRASSHOPPER_SPECIES == 'APSP'| GRASSHOPPER_SPECIES == 'DIVI'  ) %>% 
  dplyr::select(GRASSHOPPER_SPECIES, ARBE:VAMY) %>% 
  rowwise() %>% 
  mutate(Consumed=sum(c_across(2:17)))

nmat2 <- (nmat[2:17]/nmat$Consumed)^2
nmat3 <- nmat2 %>% rowwise() %>% 
  mutate(B=1/(sum(c_across(everything()))), Ba=(B-1)/(16-1))
exp3$B <- nmat3$B
exp3$Ba <- nmat3$Ba

library(janitor)
nmat4 <- t(nmat[1:16]) %>% row_to_names(row_number = 1)


#RUN A MANOVA TO SEE IF GHOPP SPECIES DIFFERENTIATE THEIR FEEDING NICHES BASED ON EACH PLANT TRAIT
fniche <- manova(cbind(SLA, LDMC,CN_RATIO) ~ SUBFAMILY+GRASSHOPPER_SPECIES, data=exp3)
summary(fniche)
summary.aov(fniche)

fn_sla <- glmmTMB(CN_RATIO ~ GRASSHOPPER_SPECIES + (1|SUBFAMILY), data=exp3)
Anova(fn_sla)
summary(fn_sla)

#THEY DO!
fig1 <- ggplot(chem2plot , aes(fill=compound, x=paste(Region,Chemo,sep="-"), y=log(conc+1))) + 
  geom_bar(position="stack",stat="identity")+
  #facet_wrap(~Season) + 
  #scale_fill_viridis(discrete=T,labels=c('Other','α-Thujene','Myrcene','Octen-3-ol','α-Terpinene', 'γ-Terpinene',
  #                                      'p-Cymene','Thymoquinone','Carvacrol','Thymol'), name="Compound")+
  scale_fill_jco(labels=c('Other','α-Thujene','Myrcene','Octen-3-ol','α-Terpinene', 'γ-Terpinene',
                          'p-Cymene','Thymoquinone','Carvacrol','Thymol'), name="Compound")+
  #ggtitle("B")+
  xlab('Origin - Chemo')+ ylab('Terpene concentration (log(mg/g))')+
  theme_bw(base_size = 24)
fig1

## reshape long
library(ggsci)
exp3long <- exp3 %>% select(GRASSHOPPER_SPECIES, MESOCOSM,SUBFAMILY, ARBE:VAMY) %>% 
  group_by(GRASSHOPPER_SPECIES, MESOCOSM,SUBFAMILY) %>% 
  pivot_longer(cols=ARBE:VAMY, names_to = 'Plant_spp', values_to = 'Amount_cons_mm2') %>% 
  group_by(Plant_spp,GRASSHOPPER_SPECIES,SUBFAMILY) %>% summarise(Amount_cons_mm2=mean(Amount_cons_mm2)) %>% 
  filter(Plant_spp=='ARBE'| Plant_spp=='SCSC'|Plant_spp=='SOSE'|Plant_spp=='PIGR'|Plant_spp=='QULA'|
           Plant_spp=='CHNI'|Plant_spp=='AMAR'|Plant_spp=='LEHI'|Plant_spp=='ERTO')

ggplot(exp3long, aes(fill=Plant_spp, x=GRASSHOPPER_SPECIES, y=Amount_cons_mm2))+
  geom_bar(position="fill",stat="identity")+
  scale_fill_jco()+
  theme_bw(base_size = 24)+
  facet_wrap(~SUBFAMILY, scales = 'free_x')


#NOW SUMMARIZE THE GRASSHOPPER DIETS WITHIN SPECIES AND SEPARATE GHOPPS BY SUBFAMILY
exp4 <- exp3 %>% group_by(GRASSHOPPER_SPECIES, SUBFAMILY) %>% 
  dplyr::summarize(SLA_M=mean(SLA), SLA_SD=sd(SLA),
            LDMC_M=mean(LDMC), LDMC_SD=sd(LDMC),
            CN_RATIO_M=mean(CN_RATIO), CN_RATIO_SD=sd(CN_RATIO))

## PLOT THE DIFFERENT TRAIT PAIRINGS ####
#LDMC AND SLA
fn1 <- exp4 %>%
  ggplot(aes(x =CN_RATIO_M ,y = LDMC_M)) +
  geom_pointrange(data=exp4, aes(ymin=(LDMC_M-LDMC_SD),ymax=(LDMC_M+LDMC_SD))) +
  geom_pointrange(data=exp4, aes(xmin=(CN_RATIO_M-CN_RATIO_SD),xmax=(CN_RATIO_M+CN_RATIO_SD))) +
  geom_point(aes(fill=SUBFAMILY), pch=21, color = "black", stroke = 1, size=4) +
  scale_fill_viridis(discrete = T,  option = "plasma") +
  theme_bw(base_size = 16) +
  ylab("Weighted mean LDMC in diet (mg/g)") +
  xlab("Weighted mean C:N in diet ") +
  labs(fill = "Subfamily")+
  annotate(geom="text", label="A)", x=20, y=530, size=6)+
  theme(
    panel.border = element_rect(color="black", fill=NA, size=1.5),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.position = "none", # remove legend
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.text = element_text(size=15));fn1

fn2 <- exp4 %>%
  ggplot(aes(x =SLA_M ,y = LDMC_M)) +
  geom_pointrange(data=exp4, aes(ymin=(LDMC_M-LDMC_SD),ymax=(LDMC_M+LDMC_SD))) +
  geom_pointrange(data=exp4, aes(xmin=(SLA_M-SLA_SD),xmax=(SLA_M+SLA_SD))) +
  geom_point(aes(fill=SUBFAMILY), pch=21, color = "black", stroke = 1, size=4) +
  scale_fill_viridis(discrete = T,  option = "plasma") +
  theme_bw(base_size = 16) +
  ylab("Weighted mean LDMC in diet (mg/g)") +
  xlab("Weighted mean SLA in diet (mm^2/mg)") +
  labs(fill = "Subfamily")+
  annotate(geom="text", label="B)", x=6, y=530, size=6)+
  theme(
    panel.border = element_rect(color="black", fill=NA, size=1.5),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.position = "none", # remove legend
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.text = element_text(size=15));fn2


## PLOT AVAILABLE PLANT TRAIT SPACE
ptraitave<-ptrait %>% group_by(SPECIES) %>% 
  #summarise(CN_RATIO_M=mean(CN_RATIO, na.rm=T), LDMC_M=mean(LDMC,na.rm=T), SLA_M=mean(SLA,na.rm=T),
   #         CN_RATIO_SD=sd(CN_RATIO, na.rm=T), LDMC_SD=sd(LDMC,na.rm=T), SLA_SD=sd(SLA,na.rm=T))%>% 
  mutate(FuncGroup=case_when(SPECIES=='AMAR' |
                               SPECIES=='CRAR'  |
                               SPECIES=='ERAR'  |
                               SPECIES=='ERTO'  |
                               SPECIES=='MOPU'  |
                               SPECIES=='PIGR'  |
                               SPECIES=='SOOD'  |
                               SPECIES=='STSY'~'Forb',
                             SPECIES=='CHNI' |
                               SPECIES=='LEHI'  |
                               SPECIES=='TEVI'~'Legume',
                             SPECIES=='ARBE' |
                               SPECIES=='SCSC'  |
                               SPECIES=='SOSE'~'Grass',
                             SPECIES=='QULA' |
                               SPECIES=='VAMY'~'Woody'))

fn1a <- ggplot() +
  stat_ellipse(data=ptraitave%>% filter(LDMC<800), aes(x =CN_RATIO ,y = LDMC, color=FuncGroup), level=.95, linewidth=1.25) +
  geom_pointrange(data=exp4, aes(x =CN_RATIO_M ,y = LDMC_M, ymin=(LDMC_M-LDMC_SD),ymax=(LDMC_M+LDMC_SD))) +
  geom_pointrange(data=exp4, aes(x =CN_RATIO_M ,y = LDMC_M, xmin=(CN_RATIO_M-CN_RATIO_SD),xmax=(CN_RATIO_M+CN_RATIO_SD))) +
  geom_point(data=exp4, aes(x =CN_RATIO_M ,y = LDMC_M,fill=SUBFAMILY), pch=21, color = "black", stroke = 1, size=4) +
  scale_fill_viridis(discrete = T,  option = "plasma") +
  scale_color_viridis(discrete = T,  option = "D") +
  theme_bw(base_size = 16) +
  ylab("Weighted mean LDMC in diet (mg/g)") +
  xlab("Weighted mean C:N in diet ") +
  labs(fill = "Subfamily", color="Functional group")+
  annotate(geom="text", label="A)", x=13, y=610, size=6)+
  theme(
    panel.border = element_rect(color="black", fill=NA, size=1.5),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.position = "none", # remove legend
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.text = element_text(size=15));fn1a

fn2a <- ggplot() +
  stat_ellipse(data=ptraitave%>% filter(LDMC<800), aes(x =SLA ,y = LDMC, color=FuncGroup), level=.95, linewidth=1.25) +
  geom_pointrange(data=exp4, aes(x=SLA_M ,y=LDMC_M,ymin=(LDMC_M-LDMC_SD),ymax=(LDMC_M+LDMC_SD))) +
  geom_pointrange(data=exp4, aes(x=SLA_M ,y=LDMC_M,xmin=(SLA_M-SLA_SD),xmax=(SLA_M+SLA_SD))) +
  geom_point(data=exp4, aes(x=SLA_M ,y=LDMC_M, fill=SUBFAMILY), pch=21, color = "black", stroke = 1, size=4) +
  scale_fill_viridis(discrete = T,  option = "plasma") +
  #geom_point(data=ptraitave %>% filter(LDMC<800), aes(x =SLA ,y = LDMC, color=FuncGroup), pch=16, alpha=.5) +
  scale_color_viridis(discrete = T,  option = "D") +
  theme_bw(base_size = 16) +
  ylab("Weighted mean LDMC in diet (mg/g)") +
  xlab("Weighted mean SLA in diet (mm^2/mg)") +
  labs(fill = "Subfamily", color = "Functional group")+
  annotate(geom="text", label="B)", x=2, y=610, size=6)+
  theme(
    panel.border = element_rect(color="black", fill=NA, size=1.5),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.position = "none", # remove legend
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.text = element_text(size=15));fn2a

#LDMC AND SLA
fn3 <- ptrait2ave %>%
  ggplot(aes(x =CN_RATIO_M ,y = LDMC_M)) +
  geom_pointrange(data=ptrait2ave, aes(ymin=(LDMC_M-LDMC_SD),ymax=(LDMC_M+LDMC_SD))) +
  geom_pointrange(data=ptrait2ave, aes(xmin=(CN_RATIO_M-CN_RATIO_SD),xmax=(CN_RATIO_M+CN_RATIO_SD))) +
  geom_point(aes(fill=FuncGroup), pch=21, color = "black", stroke = 1, size=4) +
  scale_fill_viridis(discrete = T,  option = "D") +
  theme_bw(base_size = 16) +
  ylab("Mean LDMC offered in trials") +
  xlab("Mean C:N offered in trials") +
  labs(fill = "Functional group")+
  annotate(geom="text", label="C)", x=15, y=700, size=6)+
  theme(
    panel.border = element_rect(color="black", fill=NA, size=1.5),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.position = "none", # remove legend
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.text = element_text(size=15));fn3

fn4 <- ptrait2ave %>%
  ggplot(aes(x =SLA_M ,y = LDMC_M)) +
  geom_pointrange(data=ptrait2ave, aes(ymin=(LDMC_M-LDMC_SD),ymax=(LDMC_M+LDMC_SD))) +
  geom_pointrange(data=ptrait2ave, aes(xmin=(SLA_M-SLA_SD),xmax=(SLA_M+SLA_SD))) +
  geom_point(aes(fill=FuncGroup), pch=21, color = "black", stroke = 1, size=4) +
  scale_fill_viridis(discrete = T,  option = "D") +
  theme_bw(base_size = 16) +
  ylab("Mean LDMC offered in trials") +
  xlab("Mean SLA offered in trials") +
  labs(fill = "Functional group")+
  annotate(geom="text", label="D)", x=4, y=700, size=6)+
  theme(
    panel.border = element_rect(color="black", fill=NA, size=1.5),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.position = "none", # remove legend
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.text = element_text(size=15));fn4

feednicheplot <- fn1+fn2 +
  plot_layout(guides = "collect") & theme(legend.position = "top") 
feedplantplot <- fn3+fn4 +
  plot_layout(guides = "collect") & theme(legend.position = "top") 

feedplotnew <- feednicheplot/feedplantplot
ggsave("feedplotnew.tiff", feedplotnew, width=10, height=10, units="in", dpi=600, compression = "lzw", path="Outputs")


feednicheplota <- fn1a+fn2a +
  plot_layout(guides = "collect") & theme(legend.position = "top",
                                          legend.title = element_text(size=11), 
                                          legend.text = element_text(size=10)) 

ggsave("Figure4.tiff", feednicheplot, width=12, height=6, units="in", dpi=600, compression = "lzw", path="Outputs")
ggsave("feednicheplota.tiff", feednicheplota, width=12, height=6, units="in", dpi=600, compression = "lzw", path="Outputs")



#PLANT AND HERBIVORE TRAIT LINKAGES ####
## JOIN PLANT AND GRASSHOPPER TRAIT DATA ####
t6<-t5a %>% filter(SPECIES!='EROB',SPECIES!='AMMY',SPECIES!='MEPI') %>% rename(CN_RATIO_GRASSHOPPER=CN_RATIO)
exp5<-exp4 %>% rename(SPECIES=GRASSHOPPER_SPECIES,CN_RATIO_PLANT=CN_RATIO_M)
link<-full_join(t5a[-c(2,5,7),],exp5,by=c('SPECIES','SUBFAMILY')) %>% 
  select(SPECIES,SUBFAMILY,everything()) %>% 
  mutate(GENUS=substr(SPECIES, start=1, stop=2))

exp3a <- exp3 %>% mutate(SPECIES=GRASSHOPPER_SPECIES) %>% 
  select(SPECIES,SUBFAMILY,SEX,MESOCOSM,CN_RATIO,LDMC,SLA) 
link2<-full_join(t5a,exp3a,by=c('SPECIES','SUBFAMILY')) %>% 
  select(SPECIES,SUBFAMILY,everything())%>% 
  mutate(GENUS=substr(SPECIES, start=1, stop=2))


## multivar feeding models ####

mv_LDMC <- glmmTMB(LDMC ~ pc1 + (1|SUBFAMILY), data=link2 %>% filter(SEX!='U'))

mv_mod <- manova(cbind(CN_RATIO_PLANT, LDMC_M, SLA_M) ~ BV, data=link)
anova(mv_mod, by="terms")
plot(mv_mod)
summary(mv_mod)

Anova(mv_mod)
check_collinearity(mv_mod)
summary(mv_mod)
summary.aov(mv_mod)

## trait link models ####

lm.IS_LDMC <- glmmTMB(LDMC ~ IS + (1|SUBFAMILY), data=link2 )#%>% filter(LDMC_M<420) )

summary(lm.IS_LDMC)
coef(lm.IS_LDMC)
Anova(lm.IS_LDMC)
#anova(lm.IS_LDMC, ddf="Kenward-Roger")
r.squaredGLMM(lm.IS_LDMC)
lm.IS_LDMC_simres<-simulateResiduals(lm.IS_LDMC);plot(lm.IS_LDMC_simres)
hist(resid(lm.IS_LDMC))

lm.BV_LDMC <- glmmTMB(SLA ~ IS + (1|SUBFAMILY), data=link2 )
summary(lm.BV_LDMC)
coef(lm.BV_LDMC)
Anova(lm.BV_LDMC)
#anova(lm.IS_LDMC, ddf="Kenward-Roger")
r.squaredGLMM(lm.BV_LDMC)
simulateResiduals(lm.BV_LDMC, plot=T)
hist(resid(lm.IS_LDMC))

lm.CN_CN <- glmmTMB(CN_RATIO.y ~ CN_RATIO.x + (1|SUBFAMILY), data=link2 )
summary(lm.CN_CN)
coef(lm.CN_CN)
Anova(lm.CN_CN)
#anova(lm.IS_LDMC, ddf="Kenward-Roger")
r.squaredGLMM(lm.CN_CN)
simulateResiduals(lm.CN_CN, plot=T)
hist(resid(lm.CN_CN))

lm.bv_CN <- glmmTMB(CN_RATIO.y ~ BV + (1|SUBFAMILY/GENUS), data=link2 )
summary(lm.bv_CN)
coef(lm.bv_CN)
Anova(lm.bv_CN)
#anova(lm.IS_LDMC, ddf="Kenward-Roger")
r.squaredGLMM(lm.bv_CN)
simulateResiduals(lm.bv_CN, plot=T)
hist(resid(lm.bv_CN))

### trait link plots ####
ISLD_tab1 <- ggpredict(lm.IS_LDMC, terms = c("IS [2.5:4.6, by=.01]"), 
                       ci.lvl = .95, type = "fixed") %>% as.data.frame() 

ghoppIS_LDMC <- link %>% ggplot() +
  geom_line(data=ISLD_tab1, aes(x=x,y=predicted), lwd=1.5, lty=1)+
  #geom_smooth(data=link, aes(x=IS, y=LDMC_M), se=T, method="lm")+
  geom_ribbon(data=ISLD_tab1, aes(x=x, ymin = conf.low, ymax = conf.high), alpha = .2) +
  geom_point(data=link %>% filter(LDMC_M<420), aes(x=IS, y=LDMC_M, color=SUBFAMILY), size=3) +
  geom_point(data=link %>% filter(LDMC_M>420), aes(x=IS, y=LDMC_M, color=SUBFAMILY), size=3, shape=17) +  
  scale_color_viridis(discrete=T, option = 'plasma') + 
  #guides(color = guide_legend(override.aes=list(shape = 16))) +
  theme_bw(base_size = 16) +
  xlab("log Grasshopper Incisor strength") +
  ylab("Weighted mean Plant LDMC in diet") +
  annotate(geom="text", label="A)", x=2.5, y=475, size=6)+
  annotate(geom="text", label="R2m=0.26", x=3, y=440, size=5)+
  annotate(geom="text", label="R2c=0.58", x=3, y=425, size=5)+
  annotate(geom="text", label="p = 0.04", x=3, y=410, size=5)+
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none", # remove legend
    legend.text = element_text(size=15));ghoppIS_LDMC

ghoppIS_LDMC <- link2 %>% ggplot() +
  geom_line(data=ISLD_tab1, aes(x=x,y=predicted), lwd=1.5, lty=1)+
  #geom_smooth(data=link, aes(x=IS, y=LDMC_M), se=T, method="lm")+
  geom_ribbon(data=ISLD_tab1, aes(x=x, ymin = conf.low, ymax = conf.high), alpha = .2) +
  geom_point(data=link2 , aes(x=IS, y=LDMC, color=SUBFAMILY), size=3) +
  #geom_point(data=link2 %>% filter(LDMC>500), aes(x=IS, y=LDMC, color=SUBFAMILY), size=3, shape=17) +  
  scale_color_viridis(discrete=T, option = 'plasma') + 
  #guides(color = guide_legend(override.aes=list(shape = 16))) +
  theme_bw(base_size = 16) +
  xlab("log Grasshopper Incisor strength") +
  ylab("Weighted mean Plant LDMC in diet") +
  annotate(geom="text", label="A)", x=2.5, y=575, size=6)+
  annotate(geom="text", label="R2m=0.05", x=3, y=520, size=5)+
  annotate(geom="text", label="R2c=0.20", x=3, y=500, size=5)+
  annotate(geom="text", label="p = 0.02", x=3, y=480, size=5)+
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none", # remove legend
    legend.text = element_text(size=15));ghoppIS_LDMC

### CN and LDMC ####
CN_tab1 <- ggpredict(lm.CN_CN, terms = c("CN_RATIO.x [4.09:4.6, by=.01]"), 
                     ci.lvl = .95, type = "fixed") %>% as.data.frame() 

ghoppCN_CN <- link2 %>% ggplot() +
  geom_line(data=CN_tab1, aes(x=x,y=predicted), lwd=1.5, lty=1)+
  #geom_smooth(data=link, aes(x=CN_RATIO, y=CN_RATIO_PLANT), se=T, method="lm")+
  geom_ribbon(data=CN_tab1, aes(x=x, ymin = conf.low, ymax = conf.high), alpha = .2) +
  geom_point(data=link2, aes(x=CN_RATIO.x, y=CN_RATIO.y, color=SUBFAMILY), size=3) +
  scale_color_viridis(discrete=T, option = 'plasma') + 
  #guides(color = guide_legend(override.aes=list(shape = 16))) +
  theme_bw(base_size = 16) +
  xlab("Grasshopper C:N Ratio") +
  ylab("Weighted mean Plant C:N Ratio in diet") +
  xlim(4.075,4.61) +
  annotate(geom="text", label="B)", x=4.075, y=44, size=6)+
  annotate(geom="text", label="R2m=0.08", x=4.14, y=25, size=5)+
  annotate(geom="text", label="R2c=0.47", x=4.14, y=22.5, size=5)+
  annotate(geom="text", label="p < 0.01", x=4.14, y=20, size=5)+
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none", # remove legend
    legend.text = element_text(size=15));ghoppCN_CN

## body size CN plant
bv_tab1 <- ggpredict(lm.bv_CN, terms = c("BV [4:7.1, by=.01]"), 
                     ci.lvl = .95, type = "fixed") %>% as.data.frame() 

ghoppBV_CN <- link2 %>% ggplot() +
  geom_line(data=bv_tab1, aes(x=x,y=predicted), lwd=1.5, lty=1)+
  #geom_smooth(data=link, aes(x=CN_RATIO, y=CN_RATIO_PLANT), se=T, method="lm")+
  geom_ribbon(data=bv_tab1, aes(x=x, ymin = conf.low, ymax = conf.high), alpha = .2) +
  geom_point(data=link2, aes(x=BV, y=CN_RATIO.y, color=SUBFAMILY), size=3) +
  scale_color_viridis(discrete=T, option = 'plasma') + 
  #guides(color = guide_legend(override.aes=list(shape = 16))) +
  theme_bw(base_size = 16) +
  xlab("Grasshopper body volume (mm3)") +
  ylab("Weighted mean Plant C:N Ratio in diet") +
  annotate(geom="text", label="C)", x=3.95, y=43, size=6)+
  annotate(geom="text", label="R2m=0.03", x=5.2, y=42.5, size=5)+
  annotate(geom="text", label="R2c=0.46", x=5.2, y=40.5, size=5)+
  annotate(geom="text", label="p = 0.05", x=5.2, y=38.5, size=5)+
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin=unit(c(2,1,1,1),"cm"),
    legend.position = c(-1,1.1), # remove legend
    legend.direction = "horizontal",
    legend.text = element_text(size=15));ghoppBV_CN

linkplot <- ghoppIS_LDMC + ghoppCN_CN+ghoppBV_CN
  plot_layout(guides = "collect") #& theme(legend.position = "top")

ggsave("linkplot.tiff", linkplot, width=15, height=6, units="in", dpi=600, compression = "lzw", path="Outputs")


ggsave("feednicheplot.tiff", feednicheplot, width=10, height=5, units="in", dpi=600, compression = "lzw")
## trait link plots
IS_LDMC2 <- link2 %>% ggplot() +
  geom_smooth(data=link2 %>% filter(LDMC<550),aes(x=IS, y=LDMC), method='glm') +
  geom_point(data=link2 , aes(x=IS, y=LDMC, color=SUBFAMILY), size=4) +
  scale_color_viridis(discrete=T, option = 'plasma') + guides(color = guide_legend(override.aes=list(shape = 16))) +
  theme_bw(base_size = 25) +
  xlab("Log Grasshopper Incisor Strength (N)") +
  ylab("Weighted mean LDMC in diet (mg/g)") +
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none", # remove legend
    legend.text = element_text(size=15));IS_LDMC2




#INCISOR STRENGTH
IS_LDMC <- link %>% ggplot() +
  geom_smooth(data=link %>% filter(LDMC_M<420),aes(x=IS, y=LDMC_M), color='black', method='glm') +
  geom_point(data=link %>% filter(LDMC_M<420), aes(x=IS, y=LDMC_M, color=SUBFAMILY), size=8) +
  geom_point(data=link %>% filter(LDMC_M>420), aes(x=IS, y=LDMC_M, color=SUBFAMILY), size=8, shape=17) +
  scale_color_viridis(discrete=T, option = 'plasma') + guides(color = guide_legend(override.aes=list(shape = 16))) +
  theme_bw(base_size = 25) +
  xlab("Log Grasshopper Incisor Strength (N)") +
  ylab("Weighted mean LDMC in diet (mg/g)") +
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none", # remove legend
    legend.text = element_text(size=15));IS_LDMC

PC_LDMC <- link %>% ggplot() +
  geom_smooth(data=link ,aes(x=pc2, y=LDMC_M), color='black', method='glm') +
  geom_point(data=link , aes(x=pc1, y=LDMC_M, color=SUBFAMILY), size=8) +
  #geom_point(data=link %>% filter(LDMC_M>420), aes(x=IS, y=LDMC_M, color=SUBFAMILY), size=8, shape=17) +
  scale_color_viridis(discrete=T, option = 'plasma') + guides(color = guide_legend(override.aes=list(shape = 16))) +
  theme_bw(base_size = 25) +
  xlab("Grasshopper PC Axis 1") +
  ylab("Weighted mean LDMC in diet (mg/g)") +
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none", # remove legend
    legend.text = element_text(size=15));PC_LDMC





IS_SLA <- link %>% ggplot() +
  geom_smooth(data=link,aes(x=IS, y=SLA_M), color='black', method='glm') +
  geom_point(data=link, aes(x=IS, y=SLA_M, color=SUBFAMILY), size=5) +
  scale_color_viridis(discrete=T, option = 'plasma') + guides(color = guide_legend(override.aes=list(shape = 16))) +
  theme_bw(base_size = 25) +
  xlab("Log Grasshopper Incisor Strength (N)") +
  ylab("Weighted mean SLA in diet (mm^2/mg)") +
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size=15));IS_SLA

lm.IS_SLA <- glmmTMB(SLA_M ~ IS + (1|SUBFAMILY), data=link)
summary(lm.IS_SLA)
Anova(lm.IS_SLA)
lm.IS_SLA_simres<-simulateResiduals(lm.IS_SLA);plot(lm.IS_SLA_simres)

IS_CN <- link %>% ggplot() +
  geom_smooth(data=link,aes(x=IS, y=CN_RATIO_PLANT), color='black', method='glm') +
  geom_point(data=link, aes(x=IS, y=CN_RATIO_PLANT, color=SUBFAMILY), size=5) +
  scale_color_viridis(discrete=T, option = 'plasma') + guides(color = guide_legend(override.aes=list(shape = 16))) +
  theme_bw(base_size = 25) +
  xlab("Log Grasshopper Incisor Strength (N)") +
  ylab("Weighted mean C:N Ratio in diet") +
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size=15));IS_CN

lm.IS_CN <- glmmTMB(CN_RATIO_PLANT ~ IS, data=link)
summary(lm.IS_CN)
Anova(lm.IS_CN)
lm.IS_CN_simres<-simulateResiduals(lm.IS_CN);plot(lm.IS_CN_simres)


#BODY VOLUME
BV_LDMC <- link %>% ggplot() +
  geom_smooth(data=link %>% filter(LDMC_M<420),aes(x=BV, y=LDMC_M), color='black', method='glm') +
  geom_point(data=link %>% filter(LDMC_M<420), aes(x=BV, y=LDMC_M, color=SUBFAMILY), size=5) +
  geom_point(data=link %>% filter(LDMC_M>420), aes(x=BV, y=LDMC_M, color=SUBFAMILY), size=5, shape=17) +
  scale_color_viridis(discrete=T, option = 'plasma') + guides(color = guide_legend(override.aes=list(shape = 16))) +
  theme_bw(base_size = 25) +
  xlab("Weighted Mean BV/Body Volume (mm/mm^3)") +
  ylab("Weighted mean LDMC in diet (mg/g)") +
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none", # remove legend
    legend.text = element_text(size=15));BV_LDMC

lm.BV_LDMC <- glmmTMB(LDMC_M ~ BV + IS  + (1|SUBFAMILY), data=link %>% filter(LDMC_M<420))
summary(lm.BV_LDMC)
Anova(lm.BV_LDMC)
r.squaredGLMM(lm.BV_LDMC)
check_collinearity(lm.BV_LDMC)
lm.BV_LDMC_simres<-simulateResiduals(lm.BV_LDMC);plot(lm.BV_LDMC_simres)

mydf <- ggpredict(lm.BV_LDMC)
plot(mydf)

ggsave("BV_LDMC.tiff", BV_LDMC, width=10, height=10, units="in", dpi=600, compression = "lzw")

BV_SLA <- link %>% ggplot() +
  geom_smooth(data=link,aes(x=BV, y=SLA_M), color='black', method='glm') +
  geom_point(data=link, aes(x=BV, y=SLA_M, color=SUBFAMILY), size=5) +
  scale_color_viridis(discrete=T, option = 'plasma') + guides(color = guide_legend(override.aes=list(shape = 16))) +
  theme_bw(base_size = 25) +
  xlab("Weighted Mean BV/Body Volume (mm/mm^3)") +
  ylab("Weighted mean SLA in diet (mm^2/mg)") +
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size=15));BV_SLA

lm.BV_SLA <- glmmTMB(SLA_M ~ BV + (1|SUBFAMILY), data=link)
summary(lm.BV_SLA)
Anova(lm.BV_SLA)
lm.BV_SLA_simres<-simulateResiduals(lm.BV_SLA);plot(lm.BV_SLA_simres)

BV_CN <- link %>% ggplot() +
  geom_smooth(data=link,aes(x=BV, y=CN_RATIO_PLANT), color='black', method='glm') +
  geom_point(data=link, aes(x=BV, y=CN_RATIO_PLANT, color=SUBFAMILY), size=5) +
  scale_color_viridis(discrete=T, option = 'plasma') + guides(color = guide_legend(override.aes=list(shape = 16))) +
  theme_bw(base_size = 25) +
  xlab("Log Grasshopper Body Volume (mm/mm^3)") +
  ylab("Weighted mean C:N Ratio in diet") +
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size=15));BV_CN

lm.BV_CN <- glmmTMB(CN_RATIO_PLANT ~ BV + (1|SUBFAMILY), data=link)
summary(lm.BV_CN)
Anova(lm.BV_CN)
lm.BV_CN_simres<-simulateResiduals(lm.BV_CN);plot(lm.BV_CN_simres)

#GRASSHOPPER C:N RATIO
lm.ghoppCN_LDMC <- glmmTMB(LDMC_M ~ CN_RATIO_GRASSHOPPER + (1|SUBFAMILY), data=link %>% filter(LDMC_M<420))
summary(lm.ghoppCN_LDMC)
Anova(lm.ghoppCN_LDMC)
r.squaredGLMM(lm.ghoppCN_LDMC)
lm.ghoppCN_LDMC_simres<-simulateResiduals(lm.ghoppCN_LDMC);plot(lm.ghoppCN_LDMC_simres)

hist(resid(lm.ghoppCN_LDMC))

CN_tab1 <- ggpredict(lm.ghoppCN_CN, terms = c("CN_RATIO_GRASSHOPPER [4:4.6, by=.01]"), 
                     ci.lvl = .95, type = "fe") %>% as.data.frame() 

ghoppCN_LDMC <- link %>% ggplot() +
  geom_smooth(data=link %>% filter(LDMC_M<420),aes(x=CN_RATIO_GRASSHOPPER, y=LDMC_M), color='black', method='glm') +
  geom_point(data=link %>% filter(LDMC_M<420), aes(x=CN_RATIO_GRASSHOPPER, y=LDMC_M, color=SUBFAMILY), size=5) +
  geom_point(data=link %>% filter(LDMC_M>420), aes(x=CN_RATIO_GRASSHOPPER, y=LDMC_M, color=SUBFAMILY), size=5, shape=17) +
  scale_color_viridis(discrete=T, option = 'plasma') + guides(color = guide_legend(override.aes=list(shape = 16))) +
  theme_bw(base_size = 25) +
  xlab("Grasshopper C:N Ratio") +
  ylab("Weighted mean LDMC in diet (mg/g)") +
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none", # remove legend
    legend.text = element_text(size=15));ghoppCN_LDMC



ggsave("ghoppCN_LDMC.tiff", ghoppCN_LDMC, width=10, height=10, units="in", dpi=600, compression = "lzw")

ghoppCN_SLA <- link %>% ggplot() +
  geom_point(data=link, aes(x=CN_RATIO_GRASSHOPPER, y=SLA_M, color=SUBFAMILY), size=8) +
  scale_color_viridis(discrete=T, option = 'plasma') + guides(color = guide_legend(override.aes=list(shape = 16))) +
  theme_bw(base_size = 25) +
  xlab("Grasshopper C:N Ratio") +
  ylab("Weighted mean SLA in diet (mm^2/mg)") +
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none", # remove legend
    legend.text = element_text(size=15));ghoppCN_SLA

ggsave("ghoppCN_SLA.tiff", ghoppCN_SLA, width=10, height=10, units="in", dpi=600, compression = "lzw")

lm.ghoppCN_SLA <- glmmTMB(SLA_M ~ CN_RATIO_GRASSHOPPER, data=link)
summary(lm.ghoppCN_SLA)
Anova(lm.ghoppCN_SLA)
lm.ghoppCN_SLA_simres<-simulateResiduals(lm.ghoppCN_SLA);plot(lm.ghoppCN_SLA_simres)

### CN ghop ~ CN plant ####
ghoppCN_CN <- link %>% ggplot() +
  geom_smooth(data=link ,aes(x=CN_RATIO_GRASSHOPPER, y=CN_RATIO_PLANT), color='black', method='glm') +
  geom_point(data=link, aes(x=CN_RATIO_GRASSHOPPER, y=CN_RATIO_PLANT, color=SUBFAMILY), size=8) +
  scale_color_viridis(discrete=T, option = 'plasma') + guides(color = guide_legend(override.aes=list(shape = 16))) +
  theme_bw(base_size = 25) +
  xlab("Grasshopper C:N Ratio") +
  ylab("Weighted mean Plant C:N Ratio in diet") +
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none", # remove legend
    legend.text = element_text(size=15));ghoppCN_CN

lm.ghoppCN_CN <- glmmTMB(CN_RATIO_PLANT ~ CN_RATIO_GRASSHOPPER + (1|SUBFAMILY), data=link)
summary(lm.ghoppCN_CN)
Anova(lm.ghoppCN_CN)
r.squaredGLMM(lm.ghoppCN_CN)
lm.ghoppCN_CN_simres<-simulateResiduals(lm.ghoppCN_CN);plot(lm.ghoppCN_CN_simres)



### IS and LDMC ####
ISLD_tab1 <- ggpredict(lm.IS_LDMC, terms = c("IS [2.5:4.6, by=.01]"), 
                     ci.lvl = .95, type = "fe") %>% as.data.frame() 

ghoppIS_LDMC <- link %>% ggplot() +
  geom_line(data=ISLD_tab1, aes(x=x,y=predicted), lwd=1.5, lty=1)+
  #geom_smooth(data=link, aes(x=CN_RATIO_GRASSHOPPER, y=CN_RATIO_PLANT, color=SUBFAMILY), se=F, method="lm")+
  geom_ribbon(data=ISLD_tab1, aes(x=x, ymin = conf.low, ymax = conf.high), alpha = .2) +
  geom_point(data=link %>% filter(LDMC_M<420), aes(x=IS, y=LDMC_M, color=SUBFAMILY), size=5) +
  geom_point(data=link %>% filter(LDMC_M>420), aes(x=IS, y=LDMC_M, color=SUBFAMILY), size=5, shape=17) +  
  scale_color_viridis(discrete=T, option = 'plasma') + guides(color = guide_legend(override.aes=list(shape = 16))) +
  theme_bw(base_size = 25) +
  xlab("Grasshopper Incisor strength") +
  ylab("Weighted mean Plant LDMC in diet") +
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none", # remove legend
    legend.text = element_text(size=15));ghoppIS_LDMC

### CN and LDMC ####
CN_tab1 <- ggpredict(lm.ghoppCN_CN, terms = c("CN_RATIO_GRASSHOPPER [4:4.6, by=.01]"), 
                     ci.lvl = .95, type = "fe") %>% as.data.frame() 

ghoppCN_CN <- link %>% ggplot() +
  geom_line(data=CN_tab1, aes(x=x,y=predicted), lwd=1.5, lty=1)+
  #geom_smooth(data=link, aes(x=CN_RATIO_GRASSHOPPER, y=CN_RATIO_PLANT, color=SUBFAMILY), se=F, method="lm")+
  geom_ribbon(data=CN_tab1, aes(x=x, ymin = conf.low, ymax = conf.high), alpha = .2) +
  geom_point(data=link, aes(x=CN_RATIO_GRASSHOPPER, y=CN_RATIO_PLANT, color=SUBFAMILY), size=8) +
  scale_color_viridis(discrete=T, option = 'plasma') + guides(color = guide_legend(override.aes=list(shape = 16))) +
  theme_bw(base_size = 25) +
  xlab("Grasshopper C:N Ratio") +
  ylab("Weighted mean Plant C:N Ratio in diet") +
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none", # remove legend
    legend.text = element_text(size=15));ghoppCN_CN


##COMMUNITY SES -- OLD APE VERSION TO CHECK ####
t8 <- t7[-6,]
ses_wing <- ses.mpd(comm3, dist(t8[,c(3)], method="euclidean"), abundance.weighted=T)
ses_wing
comm5$ses_wing <- ses_wing$mpd.obs.z

ses_IS <- ses.mpd(comm3, dist(t8[,c(1)], method="euclidean"),null.model="taxa.labels", abundance.weighted=T)
ses_IS
comm5$ses_IS <- ses_IS$mpd.obs.z
summary(ses_IS$mpd.obs.z)
summary(ses_IS$mpd.obs.p)
ses_IS %>% count(if_else(mpd.obs.p<.1, 1,0))

ses_feed <- ses.mpd(comm2, dist(t7[,5], method="euclidean"), abundance.weighted=T)
ses_feed
comm5$ses_feed <- ses_feed$mpd.obs.z

ses_cn <- ses.mpd(comm3, dist(t8[,c(2)], method="euclidean"), abundance.weighted=T)
ses_cn
comm5$ses_cn <- ses_cn$mpd.obs.z
summary(ses_cn$mpd.obs.z)
summary(ses_cn$mpd.obs.p)
ses_cn %>% count(if_else(mpd.obs.p>.1, 1,0))

#ENVIRONMENTAL VARIABLES EFFECT ON SES ####
sesIS <- ggplot(comm5, aes(x=TIME_OF_YEAR, y=ses_IS, color=TIME_OF_YEAR)) + 
  geom_boxplot(outlier.shape=NA, lwd=1)+geom_jitter(height=0, width=.2) + 
  geom_smooth(method = 'glm')+
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_color_viridis(option="G", discrete=T, direction=-1, end=.6)+
  scale_y_continuous(limits=c(-3,3))+
  theme_bw(base_size = 16)

sesCN <- ggplot(comm5, aes(x=TIME_OF_YEAR, y=ses_cn, color=TIME_OF_YEAR)) + 
  geom_boxplot(outlier.shape=NA, lwd=1)+geom_jitter(height=0, width=.2) + 
  geom_smooth(method = 'glm')+
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_color_viridis(option="G", discrete=T, direction=-1, end=.6)+
  scale_y_continuous(limits=c(-3,3))+
  theme_bw(base_size = 16)

sesWing <- ggplot(comm5, aes(x=TIME_OF_YEAR, y=ses_wing, color=TIME_OF_YEAR)) + 
  geom_boxplot(outlier.shape=NA, lwd=1)+geom_jitter(height=0, width=.2) + 
  geom_smooth(method = 'glm')+
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_color_viridis(option="G", discrete=T, direction=-1, end=.6)+
  scale_y_continuous(limits=c(-3,3))+
  theme_bw(base_size = 16)


sesplot <- sesIS+sesCN+sesWing+ plot_layout(guides='collect') & theme(legend.position='top')
ggsave("sesplot.tiff", sesplot, width=12, height=5, units="in", dpi=600, compression = "lzw")


ggplot(comm5, aes(x=TIME_OF_YEAR, y=ses_cn, color=TIME_OF_YEAR)) + 
  geom_boxplot()+geom_jitter(height=0, width=.2) + 
  geom_smooth(method = 'glm')+geom_hline(yintercept = 0, linetype="dashed")+
  theme_bw(base_size = 16)

ggplot(comm5, aes(x=TIME_SINCE_FIRE_DAYS, y=IS, color=TIME_OF_YEAR)) + 
  geom_point() + geom_smooth(method = 'glm') 

ggplot(comm5, aes(x=TIME_OF_YEAR, y=ses_wing, color=TIME_OF_YEAR)) + 
  geom_boxplot()+geom_jitter(height=0, width=.2) + 
  geom_smooth(method = 'glm')+geom_hline(yintercept = 0, linetype="dashed")+
  theme_bw(base_size = 16)

ggplot(comm5, aes(x=TIME_SINCE_FIRE_DAYS, y=ses_wing, color=TIME_OF_YEAR)) + 
  geom_point() + geom_smooth(method = 'glm', formula = y~x+I(x^2))

ggplot(comm5, aes(x=TIME_OF_YEAR, y=ses_feed, color=TIME_OF_YEAR)) + 
  geom_boxplot()+geom_jitter(height=0, width=.2) + 
  geom_smooth(method = 'glm')+geom_hline(yintercept = 0, linetype="dashed")+
  theme_bw(base_size = 16)

ggplot(comm5, aes(x=TIME_SINCE_FIRE_DAYS, y=ses_wing, color=TIME_OF_YEAR)) + 
  geom_point() + geom_smooth(method = 'glm', formula = y~x+I(x^2))

#FULL MODEL FOR SES
head(comm5)
comm_ses <- comm5 %>% select(TIME_OF_YEAR,SITE,ses_IS,ses_cn,ses_wing) %>% 
  group_by(SITE,TIME_OF_YEAR) %>% 
  pivot_longer(cols=3:5, names_to = "Trait", values_to = "SES")


lm.ses<-glmmTMB(SES~Trait*TIME_OF_YEAR+ (1|SITE/Trait), data=comm_ses)
Anova(lm.ses)
simulateResiduals(lm.ses, plot=T)
summary(lm.ses)
emmeans(lm.ses, pairwise~Trait, infer=T)

lm.ses_wing<-glmmTMB(ses_wing~TIME_OF_YEAR+(1|SITE), data=comm5)
Anova(lm.ses_wing)
check_collinearity(lm.ses_wing)
simulateResiduals(lm.ses_wing, plot=T)
summary(lm.ses_wing)
emmeans(lm.ses_wing, ~1, infer=T)


lm.ses_IS<-glmmTMB(ses_IS~TIME_OF_YEAR+(1|SITE), data=comm5)
Anova(lm.ses_IS)
check_collinearity(lm.ses_IS)
simulateResiduals(lm.ses_IS, plot=T)
summary(lm.ses_IS)
emmeans(lm.ses_IS, ~1, infer=T)

lm.ses_CN<-glmmTMB(ses_IS~TIME_OF_YEAR+(1|SITE), data=comm5)
Anova(lm.ses_IS)
check_collinearity(lm.ses_IS)
simulateResiduals(lm.ses_IS, plot=T)
summary(lm.ses_IS)
emmeans(lm.ses_IS, ~1, infer=T)

sesplot2 <- ggplot(comm_ses) + 
  geom_boxplot(aes(x=Trait, y=SES),
               outlier.shape=NA, lwd=1)+
  geom_jitter(aes(x=Trait, y=SES, color=TIME_OF_YEAR), height=0, width=.2) + 
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_color_viridis(option="G", discrete=T, direction=-1, end=.6)+
  scale_y_continuous(name="Standardize Effect Size")+
  theme_bw(base_size = 16)+
  theme(legend.position = "top")


ggsave("sesplot2.tiff", sesplot2, width=6, height=5, units="in", dpi=600, compression = "lzw")


#ENVIRONMENTAL VARIABLES EFFECT ON RAW RICHNESS
ggplot(comm5, aes(x=TIME_SINCE_FIRE_DAYS, y=richness)) + geom_point() + geom_smooth(method = 'glm') #Likely quadratic
ggplot(comm5, aes(x=FIRE_FREQUENCY, y=richness)) + geom_point() + geom_smooth() #Likely linear
#ggplot(comm5, aes(x=TEMP_M, y=richness, color=TIME_OF_YEAR)) + geom_point() + geom_smooth(method = 'glm') #Likely linear
ggplot(comm5, aes(x=TEMP_CV, y=richness)) + geom_point() + geom_smooth() #Likely linear
ggplot(comm5, aes(x=TIME_OF_YEAR, y=richness)) + geom_boxplot()

#Full model
lm.richness_all<-glmmTMB(richness~(I(TIME_SINCE_FIRE_DAYS^2)+TIME_SINCE_FIRE_DAYS+FIRE_FREQUENCY+TEMP_CV)*TIME_OF_YEAR+(1|SITE), data=comm5)
Anova(lm.richness_all)
check_collinearity(lm.richness_all)
lm.richness_all_simres<-simulateResiduals(lm.richness_all);plot(lm.richness_all_simres)
#Eliminate terms to pick the best model
lm.richness<-glmmTMB(richness~(I(TIME_SINCE_FIRE_DAYS^2)+TIME_SINCE_FIRE_DAYS+TEMP_CV)*TIME_OF_YEAR+FIRE_FREQUENCY+(1|SITE), data=comm5)
Anova(lm.richness)
lm.richness<-glmmTMB(richness~TEMP_CV*TIME_OF_YEAR+I(TIME_SINCE_FIRE_DAYS^2)+TIME_SINCE_FIRE_DAYS+FIRE_FREQUENCY+(1|SITE), data=comm5)
Anova(lm.richness)
r.squaredGLMM(lm.richness)
check_collinearity(lm.richness)
lm.richness_simres<-simulateResiduals(lm.richness);plot(lm.richness_simres)

#ENVIRONMENTAL VARIABLES EFFECT ON SHANNON RICHNESS
ggplot(comm5, aes(x=TIME_SINCE_FIRE_DAYS, y=shannon)) + geom_point() + geom_smooth(method = 'lm',formula = y~poly(x,2)) #Likely quadratic
ggplot(comm5, aes(x=FIRE_FREQUENCY, y=shannon)) + geom_point() + geom_smooth() #Likely linear
#ggplot(comm5, aes(x=TEMP_M, y=shannon)) + geom_point() + geom_smooth() #Quadratic or linear, try fitting specifics
#ggplot(comm5, aes(x=TEMP_M, y=shannon)) + geom_point() + geom_smooth(method = 'glm') #Seems likely
#ggplot(comm5, aes(x=TEMP_M, y=shannon)) + geom_point() + geom_smooth(method = 'lm',formula = y~poly(x,2)) #Maybe, fit linear to match raw richness
ggplot(comm5, aes(x=TEMP_CV, y=shannon)) + geom_point() + geom_smooth(method = 'glm') #Likely linear
ggplot(comm5, aes(x=TIME_OF_YEAR, y=shannon)) + geom_boxplot()

#Full model
lm.shannon_all<-glmmTMB(shannon~(I(TIME_SINCE_FIRE_DAYS^2)+TIME_SINCE_FIRE_DAYS+FIRE_FREQUENCY+TEMP_CV)*TIME_OF_YEAR+(1|SITE), data=comm5)
Anova(lm.shannon_all)
check_collinearity(lm.shannon_all)
lm.shannon_all_simres<-simulateResiduals(lm.shannon_all);plot(lm.shannon_all_simres)
#Eliminate terms to pick the best model
lm.shannon<-glmmTMB(shannon~(I(TIME_SINCE_FIRE_DAYS^2)+TIME_SINCE_FIRE_DAYS+TEMP_CV)*TIME_OF_YEAR+FIRE_FREQUENCY+(1|SITE), data=comm5)
Anova(lm.shannon)
lm.shannon<-glmmTMB(shannon~TEMP_CV*TIME_OF_YEAR+I(TIME_SINCE_FIRE_DAYS^2)+TIME_SINCE_FIRE_DAYS+FIRE_FREQUENCY+(1|SITE), data=comm5)
Anova(lm.shannon)
lm.shannon<-glmmTMB(shannon~TEMP_CV*TIME_OF_YEAR+I(TIME_SINCE_FIRE_DAYS^2)+TIME_SINCE_FIRE_DAYS+(1|SITE), data=comm5)
Anova(lm.shannon)
summary(lm.shannon)
r.squaredGLMM(lm.shannon)
check_collinearity(lm.shannon)
lm.shannon_simres<-simulateResiduals(lm.shannon);plot(lm.shannon_simres)

#Shannon Index Figures
shannon_tsfd<-comm5 %>% ggplot() + 
  geom_smooth(data=comm5,aes(x=exp(TIME_SINCE_FIRE_DAYS), y=shannon),
              method = 'glm',
              formula = y~poly(log(x),2),
              size=2,
              color='red') +
  geom_point(data=comm5,aes(x=exp(TIME_SINCE_FIRE_DAYS), y=shannon), size=4) +
  scale_x_continuous(trans = scales::log_trans(),
                     breaks = c(seq(50,3550,by=100)),
                     labels = c("50","150","","350","","550","","","","","","","","","","1550","","","","","","","","","","","","","","","","","","","","3550")) +
  geom_vline(aes(xintercept=350), color='black', linetype='dashed', size=2) +
  theme_bw(base_size = 25) +
  xlab("Time Since Fire (days)") +
  ylab("Shannon's Diversity Index") +
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank());shannon_tsfd

ggsave("shannon_tsfd.tiff", shannon_tsfd, width=10, height=10, units="in", dpi=600, compression = "lzw")

shannon_tcv<-comm5 %>% ggplot() + 
  geom_smooth(data=comm5,aes(x=TEMP_CV, y=shannon, color=TIME_OF_YEAR),
              method = 'glm',
              size=2) +
  geom_point(data=comm5,aes(x=TEMP_CV, y=shannon, color=TIME_OF_YEAR), size=4) +
  scale_color_manual(values=c('EARLY'='deepskyblue2','LATE'='darkblue'), labels=c('EARLY'='July','LATE'='September')) +
  theme_bw(base_size = 25) +
  xlab("Coefficient of Variance of Site Temp (°C)") +
  ylab("Shannon's Diversity Index") +
  labs(color='Time of Year') +
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank());shannon_tcv

ggsave("shannon_tcv.tiff", shannon_tcv, width=12.5, height=10, units="in", dpi=600, compression = "lzw")

shannon_figures<-shannon_tsfd+shannon_tcv +
  plot_annotation(tag_levels = 'A',tag_suffix = ')');shannon_figures
ggsave("shannon_figures.tiff", shannon_figures, width=22.5, height=10, units="in", dpi=600, compression = "lzw")

#ENVIRONMENTAL VARIABLES EFFECT ON INVERSE SIMPSONS RICHNESS
ggplot(comm5, aes(x=TIME_SINCE_FIRE_DAYS, y=simpsonsinv, color=TIME_OF_YEAR)) + geom_point() + geom_smooth(method = 'lm',formula = y~poly(x,2)) #Likely quadratic
ggplot(comm5, aes(x=FIRE_FREQUENCY, y=simpsonsinv)) + geom_point() + geom_smooth() #Likely linear
#ggplot(comm5, aes(x=TEMP_M, y=simpsonsinv)) + geom_point() + geom_smooth() #Same as Shannon, try linear and quadratic
#ggplot(comm5, aes(x=TEMP_M, y=simpsonsinv)) + geom_point() + geom_smooth(method = 'glm') #Seems likely
#ggplot(comm5, aes(x=TEMP_M, y=simpsonsinv)) + geom_point() + geom_smooth(method = 'lm',formula = y~poly(x,2)) #More likely quadratic, hard to say
ggplot(comm5, aes(x=TEMP_CV, y=simpsonsinv)) + geom_point() + geom_smooth() #Likely linear
ggplot(comm5, aes(x=TIME_OF_YEAR, y=simpsonsinv)) + geom_boxplot()

#Full model
lm.simpsonsinv_all<-glmmTMB(simpsonsinv~(I(TIME_SINCE_FIRE_DAYS^2)+TIME_SINCE_FIRE_DAYS+FIRE_FREQUENCY+TEMP_CV)*TIME_OF_YEAR+(1|SITE), data=comm5)
Anova(lm.simpsonsinv_all)
check_collinearity(lm.simpsonsinv_all)
lm.simpsonsinv_all_simres<-simulateResiduals(lm.simpsonsinv_all);plot(lm.simpsonsinv_all_simres)
#Eliminate terms to pick the best model
lm.simpsonsinv<-glmmTMB(simpsonsinv~(I(TIME_SINCE_FIRE_DAYS^2)+TIME_SINCE_FIRE_DAYS+FIRE_FREQUENCY)*TIME_OF_YEAR+TEMP_CV+(1|SITE), data=comm5)
Anova(lm.simpsonsinv)
lm.simpsonsinv<-glmmTMB(simpsonsinv~(I(TIME_SINCE_FIRE_DAYS^2)+TIME_SINCE_FIRE_DAYS)*TIME_OF_YEAR+FIRE_FREQUENCY+TEMP_CV+(1|SITE), data=comm5)
Anova(lm.simpsonsinv)
lm.simpsonsinv<-glmmTMB(simpsonsinv~(I(TIME_SINCE_FIRE_DAYS^2)+TIME_SINCE_FIRE_DAYS)*TIME_OF_YEAR+TEMP_CV+(1|SITE), data=comm5)
Anova(lm.simpsonsinv)
r.squaredGLMM(lm.simpsonsinv)
check_collinearity(lm.simpsonsinv)
lm.simpsonsinv_simres<-simulateResiduals(lm.simpsonsinv);plot(lm.simpsonsinv_simres)

#Inverse Simpson Figures
simpsonsinv_tsfd<-comm5 %>% ggplot() + 
  geom_smooth(data=comm5,aes(x=exp(TIME_SINCE_FIRE_DAYS), y=simpsonsinv, color=TIME_OF_YEAR),
              method = 'glm',
              formula = y~poly(log(x),2),
              size=2) +
  geom_point(data=comm5,aes(x=exp(TIME_SINCE_FIRE_DAYS), y=simpsonsinv, color=TIME_OF_YEAR), size=4) +
  scale_color_manual(values=c('EARLY'='#ffa600','LATE'='#00a1eb'), labels=c('EARLY'='July','LATE'='September')) +
  theme_bw(base_size = 25) +
  xlab("Time Since Fire (days)") +
  ylab("Inverse Simpson's Index") +
  labs(color='Time of Year') +
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank());simpsonsinv_tsfd

ggsave("simpsonsinv_tsfd.tiff", simpsonsinv_tsfd, width=12.5, height=10, units="in", dpi=600, compression = "lzw")
  
simpsonsinv_tempm<-comm5 %>% ggplot() + 
  geom_smooth(data=comm5,aes(x=TEMP_M, y=simpsonsinv),method = 'glm',formula = y~poly(x,2),size=2,color='forestgreen') +
  geom_point(data=comm5,aes(x=TEMP_M, y=simpsonsinv), size=4) +
  theme_bw(base_size = 25) +
  xlab("Mean Site Temperature (°C)") +
  ylab("Inverse Simpson's Index") +
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank());simpsonsinv_tempm

ggsave("simpsonsinv_tempm.tiff", simpsonsinv_tempm, width=10, height=10, units="in", dpi=600, compression = "lzw")

#CCA TO EVALUATE HOW ALL PREDICTORS EFFECT RESPONSES
full.dat<-comm5 %>% 
  select(SITE,IS:FEEDING_GUILD,TIME_SINCE_FIRE_DAYS,FIRE_FREQUENCY,TEMP_M,TEMP_CV,TIME_OF_YEAR) %>% 
  drop_na()%>% 
  as.data.frame();full.dat
traits<-full.dat %>% 
  select(IS:FEEDING_GUILD) %>% 
  as.data.frame();traits
env<-full.dat %>% 
  select(TIME_SINCE_FIRE_DAYS,FIRE_FREQUENCY,TEMP_M,TEMP_CV,TIME_OF_YEAR) %>% 
  as.data.frame();env
cca.res<-cca(traits,env)
summary(cca.res)

plot(cca.res)
cca.anova<-anova(cca.res);cca.anova
envfit_result <- envfit(cca.res, env);envfit_result

# Plot the CCA with significance labels
ordiplot(cca.res, display = "sites")
plot(envfit_result, p.max = 0.05, col = "red", add = TRUE, cex = 0.8)

# Extract significant axes
significant_axes <- envfit_result$vectors$signif
# Display the significant axes
significant_axes

vif.cca(cca.res)

#ENVIRONMENTAL VARIABLES EFFECT ON INCISOR STRENGTH
ggplot(comm5, aes(x=TIME_SINCE_FIRE_DAYS, y=IS)) + geom_point() + geom_smooth() #Likely linear
ggplot(comm5, aes(x=FIRE_FREQUENCY, y=IS)) + geom_point() + geom_smooth() #Likely linear
ggplot(comm5, aes(x=TEMP_M, y=IS)) + geom_point() + geom_smooth() #Likely linear
ggplot(comm5, aes(x=TEMP_CV, y=IS)) + geom_point() + geom_smooth() #Likely linear
ggplot(comm5, aes(x=TIME_OF_YEAR, y=IS)) + geom_boxplot()

#Full model
m.is_all<-glmmTMB(IS~(TIME_SINCE_FIRE_DAYS+FIRE_FREQUENCY+TEMP_CV)*TIME_OF_YEAR+(1|SITE), data=comm5)
Anova(m.is_all)
emmeans(m.is_all,pairwise~TIME_OF_YEAR,type='response',infer=T)
ggplot(comm5, aes(x=TEMP_CV, y=IS, color = TIME_OF_YEAR)) + geom_point() + geom_smooth(method = 'glm')
ggplot(comm5, aes(x=TIME_SINCE_FIRE_DAYS, y=IS, color = TIME_OF_YEAR)) + geom_point() + geom_smooth(method = 'glm')
check_collinearity(m.is_all)
m.is_all_simres<-simulateResiduals(m.is_all);plot(m.is_all_simres)
#Eliminate terms to pick the best model
m.is<-glmmTMB(IS~(TIME_SINCE_FIRE_DAYS+FIRE_FREQUENCY)*TIME_OF_YEAR+TEMP_CV+(1|SITE), data=comm5)
Anova(m.is)
m.is<-glmmTMB(IS~TIME_SINCE_FIRE_DAYS*TIME_OF_YEAR+FIRE_FREQUENCY+TEMP_CV+(1|SITE), data=comm5)
Anova(m.is)
m.is<-glmmTMB(IS~TIME_SINCE_FIRE_DAYS+FIRE_FREQUENCY+TEMP_CV+TIME_OF_YEAR+(1|SITE), data=comm5)
Anova(m.is)
m.is<-glmmTMB(IS~TIME_SINCE_FIRE_DAYS+TEMP_CV+TIME_OF_YEAR+(1|SITE), data=comm5)
Anova(m.is)
r.squaredGLMM(m.is)
emmeans(m.is,pairwise~TIME_OF_YEAR,type='response',infer=T)
check_collinearity(m.is)
m.is_simres<-simulateResiduals(m.is);plot(m.is_simres)

#ISvEnv Figures
IS_Env_tsfd<-comm5 %>% ggplot() + 
  geom_smooth(data=comm5,aes(x=exp(TIME_SINCE_FIRE_DAYS), y=IS),method='glm',size=2, formula = y~x, color='red') +
  geom_point(data=comm5,aes(x=exp(TIME_SINCE_FIRE_DAYS), y=IS), size=4) +
  scale_x_continuous(trans = scales::log_trans(),
                     breaks = c(seq(50,3550,by=100)),
                     labels = c("50","150","","","","550","","","","","","","","","","1550","","","","","","","","","","","","","","","","","","","","3550")) +
  geom_vline(aes(xintercept=350), color='black', linetype='dashed', size=2) +
  theme_bw(base_size = 25) +
  xlab("Time Since Fire (days)") +
  ylab("Weighted Mean Incisor Strength (N)") +
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank());IS_Env_tsfd

ggsave("IS_Env_tsfd.tiff", IS_Env_tsfd, width=10, height=10, units="in", dpi=600, compression = "lzw")

IS_Env_tempcv<-comm5 %>% ggplot() + 
  geom_smooth(data=comm5,aes(x=TEMP_CV, y=IS),method='glm',size=2,color='forestgreen') +
  geom_point(data=comm5,aes(x=TEMP_CV, y=IS), size=4) +
  theme_bw(base_size = 25) +
  xlab("Coefficient of Variance of Site Temperature (°C)") +
  ylab("Weighted Mean Incisor Strength (N)") +
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank());IS_Env_tempcv

ggsave("IS_Env_tempcv.tiff", IS_Env_tempcv, width=10, height=10, units="in", dpi=600, compression = "lzw")

IS_Env_ToY<-comm5 %>% ggplot() + 
  geom_boxplot(data=comm5,aes(x=TIME_OF_YEAR, y=IS, fill=TIME_OF_YEAR), width=0.5, outlier.color = NA) +
  geom_jitter(data=comm5,aes(x=TIME_OF_YEAR, y=IS), width=0.2, size=3) +
  scale_fill_manual(values=c('EARLY'='deepskyblue','LATE'='dodgerblue4')) +
  theme_bw(base_size = 25) +
  xlab("Time of Year") +
  ylab("Weighted Mean Incisor Strength (N)") +
  scale_x_discrete(labels = c('EARLY' = 'July', 'LATE' = 'September')) +
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.title.y = element_blank());IS_Env_ToY

ggsave("IS_Env_ToY.tiff", IS_Env_ToY, width=10, height=10, units="in", dpi=600, compression = "lzw")

IS_SES_ToY<-comm5 %>% ggplot() + 
  geom_hline(yintercept = 0, linetype="dashed", linewidth=1.25, color="grey")+
  geom_boxplot(data=comm5,aes(x=TIME_OF_YEAR, y=ses_IS, color=TIME_OF_YEAR), 
               width=0.5, outlier.color = NA, linewidth=1.25) +
  geom_jitter(data=comm5,aes(x=TIME_OF_YEAR, y=ses_IS, fill=TIME_OF_YEAR), 
              width=0.2, size=4, pch=21, stroke=2, color="black") +
  scale_color_manual(values=c('EARLY'='deepskyblue','LATE'='dodgerblue4')) +
  scale_fill_manual(values=c('EARLY'='deepskyblue','LATE'='dodgerblue4')) +
  theme_bw(base_size = 25) +
  xlab("Time of Year") +
  ylab("SES Incisor Strength") +
  scale_x_discrete(labels = c('EARLY' = 'July', 'LATE' = 'September')) +
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none");IS_SES_ToY

ggsave("IS_SES_ToY.tiff", IS_Env_ToY, width=10, height=10, units="in", dpi=600, compression = "lzw")

IS_Env_figures2<-IS_Env_tsfd+IS_SES_ToY+
  plot_annotation(tag_levels = 'A',tag_suffix = ')');IS_Env_figures2
ggsave("IS_Env_figures2.tiff", IS_Env_figures2, width=20, height=10, units="in", dpi=600, compression = "lzw")

#ENVIRONMENTAL VARIABLES EFFECT ON WINGSPAN
ggplot(comm5, aes(x=TIME_SINCE_FIRE_DAYS, y=WINGSPAN)) + geom_point() + geom_smooth() #Try linear and quadratic
ggplot(comm5, aes(x=TIME_SINCE_FIRE_DAYS, y=WINGSPAN)) + geom_point() + geom_smooth(method='glm')
#Could be true that immediately after a fire, ghopps with longer wings dispersed and only the surviving shortwinged ghopps remained. Then, as time goes on, the long winged individuals return.
ggplot(comm5, aes(x=TIME_SINCE_FIRE_DAYS, y=WINGSPAN)) + geom_point() + geom_smooth(method='lm',formula = y~poly(x,2))
#Could be true that the fire totally decimates the community, then in the following days only strong dispersing, long-winged grasshoppers are able to colonize so soon. Later, more short winged individuals return. Then as vegetation is not managed, only strong dispersers can navigate the habitat. 
#worth testing both theories.
ggplot(comm5, aes(x=FIRE_FREQUENCY, y=WINGSPAN)) + geom_point() + geom_smooth() #Likely linear
#ggplot(comm5, aes(x=TEMP_M, y=WINGSPAN)) + geom_point() + geom_smooth() #Try linear and quadratic
#ggplot(comm5, aes(x=TEMP_M, y=WINGSPAN)) + geom_point() + geom_smooth(method = 'glm')
#ggplot(comm5, aes(x=TEMP_M, y=WINGSPAN)) + geom_point() + geom_smooth(method = 'lm',formula = y~poly(x,2))
#Neither seem particularly significant. I'd fit a linear trend.
ggplot(comm5, aes(x=TEMP_CV, y=WINGSPAN)) + geom_point() + geom_smooth() #Likely linear
ggplot(comm5, aes(x=TIME_OF_YEAR, y=WINGSPAN)) + geom_boxplot()

#Full model
m.wingspan_all<-glmmTMB(WINGSPAN~(I(TIME_SINCE_FIRE_DAYS^2)+TIME_SINCE_FIRE_DAYS+FIRE_FREQUENCY+TEMP_CV)*TIME_OF_YEAR+(1|SITE), data=comm5)
Anova(m.wingspan_all)
emmeans(m.wingspan_all,pairwise~TIME_OF_YEAR,type='response',infer=T)
ggplot(comm5, aes(x=TIME_SINCE_FIRE_DAYS, y=WINGSPAN, color=TIME_OF_YEAR)) + geom_point() + geom_smooth(method='lm',formula = y~poly(x,2))
check_collinearity(m.wingspan_all)
m.wingspan_all_simres<-simulateResiduals(m.wingspan_all);plot(m.wingspan_all_simres)
#Eliminate terms to pick the best model
m.wingspan<-glmmTMB(WINGSPAN~(I(TIME_SINCE_FIRE_DAYS^2)+TIME_SINCE_FIRE_DAYS+TEMP_CV)*TIME_OF_YEAR+FIRE_FREQUENCY+(1|SITE), data=comm5)
Anova(m.wingspan)
m.wingspan<-glmmTMB(WINGSPAN~(I(TIME_SINCE_FIRE_DAYS^2)+TIME_SINCE_FIRE_DAYS)*TIME_OF_YEAR+FIRE_FREQUENCY+TEMP_CV+(1|SITE), data=comm5)
Anova(m.wingspan)
m.wingspan<-glmmTMB(WINGSPAN~(I(TIME_SINCE_FIRE_DAYS^2)+TIME_SINCE_FIRE_DAYS)*TIME_OF_YEAR+TEMP_CV+(1|SITE), data=comm5)
Anova(m.wingspan)
m.wingspan<-glmmTMB(WINGSPAN~(I(TIME_SINCE_FIRE_DAYS^2)+TIME_SINCE_FIRE_DAYS)*TIME_OF_YEAR+(1|SITE), data=comm5)
Anova(m.wingspan)
r.squaredGLMM(m.wingspan)
check_collinearity(m.wingspan)
m.wingspan_simres<-simulateResiduals(m.wingspan);plot(m.wingspan_simres)

#WingspanvsEnv Figure
Wingspan_Env_tsfd<-comm5 %>% ggplot() + 
  geom_smooth(data=comm5,aes(x=exp(TIME_SINCE_FIRE_DAYS), y=WINGSPAN, color=TIME_OF_YEAR),
              method='glm',
              formula = y~poly(log(x),2),
              size=2) +
  geom_point(data=comm5,aes(x=exp(TIME_SINCE_FIRE_DAYS), y=WINGSPAN, color=TIME_OF_YEAR), size=4) +
  scale_x_continuous(trans = scales::log_trans(),
                     breaks = c(seq(50,3550,by=100)),
                     labels = c("50","150","","","","550","","","","","","","","","","1550","","","","","","","","","","","","","","","","","","","","3550")) +
  scale_color_manual(values=c('EARLY'='deepskyblue2','LATE'='darkblue'), labels=c('EARLY'='July','LATE'='September')) +
  geom_vline(aes(xintercept=350), color='black', linetype='dashed', size=2) +
  theme_bw(base_size = 25) +
  xlab("Time Since Fire (days)") +
  ylab("Weighted Mean Wingspan/Body Volume (mm/mm^3)") +
  labs(color='Time of Year') +
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank());Wingspan_Env_tsfd

ggsave("Wingspan_Env_tsfd.tiff", Wingspan_Env_tsfd, width=12.5, height=10, units="in", dpi=600, compression = "lzw")

Wingspan_Env_tsfd_half<-comm5 %>% filter(TIME_OF_YEAR!='EARLY') %>% ggplot() + 
  geom_smooth(data=comm5,aes(x=exp(TIME_SINCE_FIRE_DAYS), y=WINGSPAN),
              color='darkblue',
              method='glm',
              formula = y~poly(log(x),2),
              size=2) +
  geom_point(data=comm5,aes(x=exp(TIME_SINCE_FIRE_DAYS), y=WINGSPAN), size=4,color='darkblue') +
  scale_x_continuous(trans = scales::log_trans(),
                     breaks = c(seq(50,3550,by=100)),
                     labels = c("50","150","","","","550","","","","","","","","","","1550","","","","","","","","","","","","","","","","","","","","3550")) +
  #scale_color_manual(values=c('EARLY'='deepskyblue2','LATE'='darkblue'), labels=c('EARLY'='July','LATE'='September')) +
  geom_vline(aes(xintercept=350), color='black', linetype='dashed', size=2) +
  theme_bw(base_size = 25) +
  xlab("Time Since Fire (days)") +
  ylab("Weighted Mean Wingspan/Body Volume (mm/mm^3)") +
  labs(color='Time of Year') +
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank());Wingspan_Env_tsfd_half

ggsave("Wingspan_Env_tsfd_half.tiff", Wingspan_Env_tsfd_half, width=10, height=10, units="in", dpi=600, compression = "lzw")

Wingspan_SES_tsfd<-comm5 %>% ggplot() +
  geom_hline(yintercept = 0, linetype="dashed", linewidth=1.5, color="grey")+
  geom_smooth(data=comm5,aes(x=exp(TIME_SINCE_FIRE_DAYS), y=ses_wing, color=TIME_OF_YEAR),
              method='glm',
              formula = y~poly(log(x),2),
              size=2) +
  geom_point(data=comm5,aes(x=exp(TIME_SINCE_FIRE_DAYS), y=ses_wing, color=TIME_OF_YEAR), size=4) +
  scale_x_continuous(trans = scales::log_trans(),
                     breaks = c(seq(50,3550,by=100)),
                     labels = c("50","150","","","","550","","","","","","","","","","1550","","","","","","","","","","","","","","","","","","","","3550")) +
  scale_color_manual(values=c('EARLY'='deepskyblue2','LATE'='darkblue'), labels=c('EARLY'='July','LATE'='September')) +
  geom_vline(aes(xintercept=350), color='black', linetype='dashed', size=2) +
  theme_bw(base_size = 25) +
  xlab("Time Since Fire (days)") +
  ylab("SES Wingspan/Body Volume") +
  labs(color='Time of Year') +
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank());Wingspan_SES_tsfd

Comm_figures2<-IS_Env_tsfd+IS_SES_ToY+Wingspan_Env_tsfd+Wingspan_SES_tsfd+
  plot_annotation(tag_levels = 'A',tag_suffix = ')')+plot_layout(guides = "collect");Comm_figures2

ggsave("Comm_figures2.tiff", Comm_figures2, width=25, height=20, units="in", dpi=600, compression = "lzw")

#ENVIRONMENTAL VARIABLES EFFECT ON HERBIVORE C:N RATIO
ggplot(comm5, aes(x=TIME_SINCE_FIRE_DAYS, y=CN_RATIO)) + geom_point() + geom_smooth() #Likely linear
ggplot(comm5, aes(x=FIRE_FREQUENCY, y=CN_RATIO)) + geom_point() + geom_smooth() #Likely linear
#ggplot(comm5, aes(x=TEMP_M, y=CN_RATIO)) + geom_point() + geom_smooth() #Likely quadratic
ggplot(comm5, aes(x=TEMP_CV, y=CN_RATIO)) + geom_point() + geom_smooth() #Likely quadratic
ggplot(comm5, aes(x=TIME_OF_YEAR, y=CN_RATIO)) + geom_boxplot()

#Full model
m.CN_all<-glmmTMB(CN_RATIO~(TIME_SINCE_FIRE_DAYS+FIRE_FREQUENCY+I(TEMP_CV^2)+TEMP_CV)*TIME_OF_YEAR+(1|SITE), data=comm5)
Anova(m.CN_all)
emmeans(m.CN_all,pairwise~TIME_OF_YEAR,type='response',infer=T)
check_collinearity(m.CN_all)
m.CN_all_simres<-simulateResiduals(m.CN_all);plot(m.CN_all_simres)
#Eliminate terms to pick the best model
m.CN<-glmmTMB(CN_RATIO~(FIRE_FREQUENCY+I(TEMP_CV^2)+TEMP_CV)*TIME_OF_YEAR+TIME_SINCE_FIRE_DAYS+(1|SITE), data=comm5)
Anova(m.CN)
m.CN<-glmmTMB(CN_RATIO~FIRE_FREQUENCY*TIME_OF_YEAR+TIME_SINCE_FIRE_DAYS+I(TEMP_CV^2)+TEMP_CV+(1|SITE), data=comm5)
Anova(m.CN)
m.CN<-glmmTMB(CN_RATIO~TIME_SINCE_FIRE_DAYS+FIRE_FREQUENCY+I(TEMP_CV^2)+TEMP_CV+TIME_OF_YEAR+(1|SITE), data=comm5)
Anova(m.CN)
m.CN<-glmmTMB(CN_RATIO~TIME_SINCE_FIRE_DAYS+FIRE_FREQUENCY+TEMP_CV+TIME_OF_YEAR+(1|SITE), data=comm5)
Anova(m.CN)
m.CN<-glmmTMB(CN_RATIO~TIME_SINCE_FIRE_DAYS+TEMP_CV+TIME_OF_YEAR+(1|SITE), data=comm5)
Anova(m.CN)
m.CN<-glmmTMB(CN_RATIO~TEMP_CV+TIME_OF_YEAR+(1|SITE), data=comm5)
Anova(m.CN)
m.CN<-glmmTMB(CN_RATIO~TIME_OF_YEAR+(1|SITE), data=comm5)
Anova(m.CN)
r.squaredGLMM(m.CN)
emmeans(m.CN,pairwise~TIME_OF_YEAR,type='response',infer=T)
check_collinearity(m.CN)
m.CN_simres<-simulateResiduals(m.CN);plot(m.CN_simres)

#Grasshopper C:N Ratio Figure
CN_Env_ToY<-comm5 %>% ggplot() + 
  geom_boxplot(data=comm5,aes(x=TIME_OF_YEAR, y=CN_RATIO, fill=TIME_OF_YEAR), width=0.5, outlier.color = NA) +
  geom_jitter(data=comm5,aes(x=TIME_OF_YEAR, y=CN_RATIO), width=0.2) +
  scale_fill_manual(values = c('EARLY'='deepskyblue2','LATE'='darkblue')) +
  theme_bw(base_size = 25) +
  xlab("Time of Year") +
  ylab("Weighted Mean Grasshopper C:N Ratio") +
  scale_x_discrete(labels = c('EARLY' = 'July', 'LATE' = 'September')) +
  theme(
    axis.text = element_text(color = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none");CN_Env_ToY

ggsave("CN_Env_ToY.tiff", CN_Env_ToY, width=10, height=10, units="in", dpi=600, compression = "lzw")

#ENVIRONMENTAL VARIABLES EFFECT ON ANTENNAL LENGTH
ggplot(comm5, aes(x=TIME_SINCE_FIRE_DAYS, y=ANTENNAL_LENGTH)) + geom_point() + geom_smooth() #Likely linear
ggplot(comm5, aes(x=FIRE_FREQUENCY, y=ANTENNAL_LENGTH)) + geom_point() + geom_smooth() #Likely linear
ggplot(comm5, aes(x=TEMP_M, y=ANTENNAL_LENGTH)) + geom_point() + geom_smooth() #Likely linear
ggplot(comm5, aes(x=TEMP_CV, y=ANTENNAL_LENGTH)) + geom_point() + geom_smooth() #Likely linear
ggplot(comm5, aes(x=TIME_OF_YEAR, y=ANTENNAL_LENGTH)) + geom_boxplot()

#Full model
m.AL_all<-glmmTMB(ANTENNAL_LENGTH~(TIME_SINCE_FIRE_DAYS+FIRE_FREQUENCY+TEMP_M+TEMP_CV)*TIME_OF_YEAR, data=comm5)
Anova(m.AL_all)
ggplot(comm5, aes(x=TIME_SINCE_FIRE_DAYS, y=ANTENNAL_LENGTH)) + geom_point() + geom_smooth(method = 'glm')
ggplot(comm5, aes(x=FIRE_FREQUENCY, y=ANTENNAL_LENGTH)) + geom_point() + geom_smooth(method = 'glm')
check_collinearity(m.AL_all)
m.AL_all_simres<-simulateResiduals(m.AL_all);plot(m.AL_all_simres)
#Eliminate terms to pick the best model
m.AL<-glmmTMB(ANTENNAL_LENGTH~TIME_SINCE_FIRE_DAYS+FIRE_FREQUENCY, data=comm5)
Anova(m.AL)
check_collinearity(m.AL)
m.AL_simres<-simulateResiduals(m.AL);plot(m.AL_simres)

#ENVIRONMENTAL VARIABLES EFFECT ON HERBIVORE FEEDING GUILD
ggplot(comm5, aes(x=TIME_SINCE_FIRE_DAYS, y=FEEDING_GUILD)) + geom_point() + geom_smooth() #Likely linear
ggplot(comm5, aes(x=FIRE_FREQUENCY, y=FEEDING_GUILD)) + geom_point() + geom_smooth() #Likely linear
#ggplot(comm5, aes(x=TEMP_M, y=FEEDING_GUILD)) + geom_point() + geom_smooth() #Likely linear
ggplot(comm5, aes(x=TEMP_CV, y=FEEDING_GUILD)) + geom_point() + geom_smooth() #Likely linear
ggplot(comm5, aes(x=TIME_OF_YEAR, y=FEEDING_GUILD)) + geom_boxplot()

#Full model
m.guild_all<-glmmTMB(FEEDING_GUILD~(TIME_SINCE_FIRE_DAYS+FIRE_FREQUENCY+TEMP_CV)*TIME_OF_YEAR+(1|SITE), data=comm5)
Anova(m.guild_all)
check_collinearity(m.guild_all)
m.guild_all_simres<-simulateResiduals(m.guild_all);plot(m.guild_all_simres)
#Eliminate terms to pick the best model
m.guild<-glmmTMB(FEEDING_GUILD~(TIME_SINCE_FIRE_DAYS+TEMP_CV)*TIME_OF_YEAR+FIRE_FREQUENCY+(1|SITE), data=comm5)
Anova(m.guild)
m.guild<-glmmTMB(FEEDING_GUILD~TEMP_CV*TIME_OF_YEAR+TIME_SINCE_FIRE_DAYS+FIRE_FREQUENCY+(1|SITE), data=comm5)
Anova(m.guild)
m.guild<-glmmTMB(FEEDING_GUILD~TIME_SINCE_FIRE_DAYS+FIRE_FREQUENCY+TEMP_CV+TIME_OF_YEAR+(1|SITE), data=comm5)
Anova(m.guild)
m.guild<-glmmTMB(FEEDING_GUILD~TIME_SINCE_FIRE_DAYS+FIRE_FREQUENCY+TIME_OF_YEAR+(1|SITE), data=comm5)
Anova(m.guild)
m.guild<-glmmTMB(FEEDING_GUILD~TIME_SINCE_FIRE_DAYS+FIRE_FREQUENCY+(1|SITE), data=comm5)
Anova(m.guild)
r.squaredGLMM(m.guild)
check_collinearity(m.guild)
m.guild_simres<-simulateResiduals(m.guild);plot(m.guild_simres)

guild_Env_tsfd<-comm5 %>% ggplot() + 
  geom_smooth(data=comm5,aes(x=exp(TIME_SINCE_FIRE_DAYS), y=FEEDING_GUILD),method='glm',size=2, color='red') +
  geom_point(data=comm5,aes(x=exp(TIME_SINCE_FIRE_DAYS), y=FEEDING_GUILD), size=4) +
  scale_x_continuous(trans = scales::log_trans(),
                     breaks = c(seq(50,3550,by=100)),
                     labels = c("50","150","","","","550","","","","","","","","","","1550","","","","","","","","","","","","","","","","","","","","3550")) +
  geom_vline(aes(xintercept=350), color='black', linetype='dashed', size=2) +
  theme_bw(base_size = 25) +
  xlab("Time Since Fire (days)") +
  ylab("Herbivore Feeding Guild") +
  theme(
    axis.text = element_text(color = 'black'),
    axis.title.y = element_text(margin = margin(r = 75)),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank());guild_Env_tsfd

ggsave("guild_Env_tsfd.tiff", guild_Env_tsfd, width=11, height=10, units="in", dpi=600, compression = "lzw")

guild_Env_ff<-comm5 %>% ggplot() + 
  geom_smooth(data=comm5,aes(x=FIRE_FREQUENCY, y=FEEDING_GUILD),method='glm',size=2,color='darkorange') +
  geom_point(data=comm5,aes(x=FIRE_FREQUENCY, y=FEEDING_GUILD), size=4) +
  theme_bw(base_size = 25) +
  xlab("Fire Frequency") +
  ylab("Herbivore Feeding Guild") +
  theme(
    axis.text = element_text(color = 'black'),
    #axis.title.y = element_text(margin = margin(r = 75)),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank());guild_Env_ff

ggsave("guild_Env_ff.tiff", guild_Env_ff, width=11, height=10, units="in", dpi=600, compression = "lzw")

guild_Env_figures<-guild_Env_tsfd+guild_Env_ff+
  plot_annotation(tag_levels = 'A',tag_suffix = ')');guild_Env_figures
ggsave("guild_Env_figures.tiff", guild_Env_figures, width=22, height=10, units="in", dpi=600, compression = "lzw")

## Community SES Analysis ####
ses_bs <- ses.mpd(mitespp, dist(size, method="euclidean"), abundance.weighted=T)

ses_bs # prints off the Starndardize effect size results (observed, randomized, z-value, etc) for each site

mite1$ses.bodysize <- ses_bs$mpd.obs.z  ## extract ses values for each site.

# SPP TRAIT ANALYSES ####
head(comm5)
trait1 <- t4 %>% select(-CN_RATIO)
head(trait1)

dat_long <- comm5 %>% rename(IS_CWM=IS, WINGSPAN_CWM=WINGSPAN, ANTENNAL_CWM=ANTENNAL_LENGTH) %>% pivot_longer(cols=3:25, names_to = "SPECIES", values_to = "Abund")

dat_long1 <- full_join(dat_long,trait1) %>% drop_na(BV)

ggplot(dat_long1 %>% filter(SPECIES!='ODAP'), aes(x=TIME_SINCE_FIRE_DAYS, y=Abund)) + 
  geom_point()+
  geom_smooth(method = "glm", se = F, method.args = list(family = "poisson"))+
  facet_wrap(~SPECIES) + theme_bw()

## SPP 4TH CORNER ANALYSIS ####
commn <- dat_long1 %>% filter(SPECIES!='ODAP') %>% select(SITE,TIME_OF_YEAR,SPECIES,Abund) %>%   
  pivot_wider(names_from = SPECIES, values_from = Abund)
commn1 <- commn[3:16] %>% select(sort(names(.)))
trait2 <- trait1 %>% filter(SPECIES!='ODAP') %>% remove_rownames %>% column_to_rownames(var="SPECIES")
env <- comm5 %>% select(SITE, TIME_OF_YEAR, TIME_SINCE_FIRE_DAYS, TEMP_M) %>% as.data.frame()

fit <- traitglm(commn1,env[2:4],trait2[1:4], family="negative.binomial", method="manyglm")

fit$fourth
plot(fit)

anova(fit, nBoot = 100)
summary(fit, nBoot=100)

a        = max( abs(fit$fourth.corner) )
colort   = colorRampPalette(c("#FDE725FF","white","#440154FF")) 
plot.4th = levelplot(t(as.matrix(fit$fourth.corner)), xlab="Environmental Variables",
                     ylab="Species traits", col.regions=colort(100), at=seq(-a, a, length=100),
                     scales = list( x= list(rot = 45)))
print(plot.4th)

## SPP SUMMARY FIGURES ####
dat_summ <- dat_long1 %>% filter(SPECIES!='ODAP') %>% group_by(SPECIES) %>% 
  summarise(TIME_SINCE_FIRE_wm=weighted.mean(TIME_SINCE_FIRE_DAYS, Abund),
            TEMP_M_wm = weighted.mean(TEMP_M, Abund),
            BV=mean(BV), IS=mean(IS),
            WINGSPAN=mean(WINGSPAN), 
            FEMUR_AREA=mean(FEMUR_AREA))

ggplot(dat_summ, aes(x=TEMP_M_wm, y=BV)) + 
  geom_point() + geom_smooth(method="lm")
