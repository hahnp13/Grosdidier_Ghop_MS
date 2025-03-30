library(tidyverse)
library(vegan) #Community Ecology
library(car) #Anova
library(DHARMa) #Model Testing: QQPlots, Predicted vs. Residuals, etc.
library(emmeans) #Least-Squared Means
library(glmmTMB) #Generalized linear mixed models
library(viridis) #Viridis Color Package
library(GGally) #GGpairs
library(performance) #Check for collinearity using vif for glmmTMB

library(nlstools)
library(bbmle)
library(nlme)
library(FSA)
library(lmtest)
library(lme4)
library(lmerTest)
library(emmeans)
library(tidyverse)
library(minpack.lm)

#COMMUNITY SAMPLING ANALYSIS####
comm<-read.csv('Data/Grasshopper ID_nymphs ided_05202023.csv') %>% filter(HABITAT!='F')
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
  select(-c(GENUS,HABITAT,SEX,AGE,IDed,MORPH.CODE,SPP.NOTES)) %>% 
  group_by(YEAR,TIME.OF.YEAR,SITE,SPECIES) %>% 
  summarise(X..INDV=sum(X..INDV)) %>% 
  pivot_wider(names_from = c(YEAR,TIME.OF.YEAR,SPECIES), values_from = X..INDV) %>% 
  replace(is.na(.),0)

# write.csv(comm0check, "Outputs/Ordway_AcacarAptsph.csv")
head(comm0check)


## SET VAGUE STARTING PARAMETERS AND ESTIMATE PARAMETERS
## USE 'nlsLM' from package:minpack.lm TO ESTIMATE STARTING PARAMETERS
a_start=.01
b_start=1
ls1.skA <- nlsLM(Y2022_LATE_APSP ~ (a * Y2022_EARLY_APSP) / (1 + b * Y2022_EARLY_APSP),
                data=comm0check, start=list(b=b_start,a=a_start))
summary(ls1.skA)

ls1.skA <- nls(Y2022_LATE_APSP ~ (a * Y2022_EARLY_APSP) / (1 + b * Y2022_EARLY_APSP),
                 data=comm0check, start=list(b=b_start,a=a_start))
summary(ls1.skA)


ls2.skA <- nls(Y2022_LATE_ACCA ~ (a * Y2022_EARLY_ACCA) / (1 + b * Y2022_EARLY_ACCA),
                 data=comm0check, start=list(b=b_start,a=a_start))
summary(ls2.skA)



plot(nlsResiduals(ls1.skA)) ## examine residuals. Most look okay. Variance increases at higher fitted values,
# but residuals typically centered around zero so models seem to estimate reasonable parameters



## 2022 Early vs Late season ####
ggplot(comm0check, aes(x=Y2022_EARLY_APSP, y=Y2022_LATE_APSP)) +
  geom_point()+
  geom_smooth(method="nls", formula=y~(a*x)/(1+b*x), 
              method.args=list(start=c(a=.1,b=5)), se=F)+
  theme_bw(base_size = 20)

ggplot(comm0check, aes(x=Y2022_EARLY_ACCA, y=Y2022_LATE_ACCA)) +
  geom_point()+
  geom_smooth(method="nls", formula=y~(a*x)/(1+b*x), 
              method.args=list(start=c(a=.1,b=1)), se=F)+
  theme_bw(base_size = 20)

ggplot(comm0check, aes(x=Y2021_EARLY_AMMY, y=Y2021_LATE_AMMY)) +
  geom_point()+
  geom_smooth(method="nls", formula= y ~ SSlogis(x, Asym, xmid, scal), se=F)+
  theme_bw(base_size = 20)

ggplot(comm0check, aes(x=Y2021_EARLY_MERO, y=Y2021_LATE_MERO)) +
  geom_point()+
  geom_smooth(method="nls", formula= y ~ SSlogis(x, Asym, xmid, scal), se=F)+
  theme_bw(base_size = 20)

ggplot(comm0check, aes(x=Y2021_EARLY_ODAP, y=Y2021_LATE_ODAP)) +
  geom_point()+
  geom_smooth(method="nls", formula= y ~ SSlogis(x, Asym, xmid, scal), 
              method.args = list(start = c(Asym = 1, xmid = .5, scal = 1)), 
              se=F)+
  theme_bw(base_size = 20)

## between years -- (POSSIBLY NOT USEFUL)

ggplot(comm0check, aes(x=Y2021_LATE_APSP, y=Y2022_EARLY_APSP)) +
  geom_point()+
  geom_smooth(method="glm", formula=y~x, method.args=list(family=poisson))+
  theme_bw(base_size = 20)


ggplot(comm0check, aes(x=Y2021_LATE_ACCA, y=Y2022_LATE_ACCA)) +
  geom_point()+
  geom_smooth(method="glm", formula=y~x, method.args=list(family=poisson))+
  theme_bw(base_size = 20)

ggplot(comm0check, aes(x=Y2021_EARLY_AMMY, y=Y2022_EARLY_AMMY)) +
  geom_point()+
  geom_smooth(method="glm", formula=y~x, method.args=list(family=poisson))+
  theme_bw(base_size = 20)

ggplot(comm0check, aes(x=Y2021_LATE_MERO, y=Y2022_EARLY_MERO)) +
  geom_point()+
  geom_smooth(method="glm", formula=y~x, method.args=list(family=poisson))+
  theme_bw(base_size = 20)

ggplot(comm0check, aes(x=Y2021_LATE_ODAP, y=Y2022_EARLY_ODAP)) +
  geom_point()+
  geom_smooth(method="glm", formula=y~x, method.args=list(family=poisson))+
  theme_bw(base_size = 20)


### chat gpt example
# Sample data (replace with actual data)
earlyabund <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
lateabund <- c(8, 15, 22, 28, 32, 35, 37, 38, 39, 40)

# Beverton-Holt model function
beverton_holt <- function(x, a, b) {
  (a * x) / (1 + b * x)
}

# Fit the model using non-linear least squares (nls)
fit <- nls(lateabund ~ beverton_holt(earlyabund, a, b),
           start = list(a = max(lateabund), b = 0.01),
           data = data.frame(earlyabund, lateabund))

# Print summary of the model
summary(fit)

# Plot the data and fitted curve
plot(earlyabund, lateabund, pch = 16, col = "blue", xlab = "Early Abundance", ylab = "Late Abundance", main = "Beverton-Holt Model Fit")
curve((coef(fit)["a"] * x) / (1 + coef(fit)["b"] * x), add = TRUE, col = "red", lwd = 2)
legend("bottomright", legend = c("Observed Data", "Fitted Model"), col = c("blue", "red"), pch = c(16, NA), lty = c(NA, 1))

