#Descriptive statistics and unlivable models for cases and controls only

library(naniar)
library(viridis)
library(visdat)
library(tidyverse)
library(Hmisc)
library(rms)
library(pROC)
library(Rcpp)
library(grid)
library(gridExtra)

rm(list=ls())
load("N:\\DataAnalysis\\CrossSectional\\CasesCX.RData")
load("N:\\DataAnalysis\\CrossSectional\\ControlsCX.RData")

d<-rbind(case_cx, con_cx)
d$case = droplevels(d$case)
d %>% group_by(case) %>% summarise(n_distinct(id))
table(d$case)
#recode outcomes 
d<- d %>%  mutate(MM = ifelse(case=="Myeloma", 1, 0)) %>% 
  relocate(MM, .after=case)
table(d$MM)
colnames(d)

#summarise missing data
vis_dat(d)
gg_miss_var(d) + theme_bw() 

#number with missing values for each variable
d_miss <- d %>% group_by(case) %>% miss_var_summary() %>% filter(!(variable %in% c("id", "age", "sex", "MM")))

#heat map for missing data 
ggplot(d_miss, aes(x=case, y=variable, fill=pct_miss)) +
  geom_tile() +
  scale_fill_viridis(discrete=FALSE, limits=c(0, 50), breaks=seq(0, 50, by=10)) +
  guides(fill=guide_colourbar(title= "% missing")) +
  ylab("Blood Test")

rm(bpd)

#boxplots all blood tests cases and controls
ggplot (data=d, aes(y=Albumin, x=case)) + geom_boxplot() + labs(x="")
ggplot (data=d, aes(y=ALP, x=case)) + geom_boxplot() + labs(x="")
ggplot (data=d, aes(y=ALT, x=case)) + geom_boxplot() + labs(x="")
ggplot (data=d, aes(y=Basophil, x=case)) + geom_boxplot() + labs(x="")
ggplot (data=d, aes(y=Calcium, x=case)) + geom_boxplot() + labs(x="")
ggplot (data=d, aes(y=CRP, x=case)) + geom_boxplot() + labs(x="")
ggplot (data=d, aes(y=Creatinine, x=case)) + geom_boxplot() + labs(x="")
ggplot (data=d, aes(y=Eosinophil, x=case)) + geom_boxplot() + labs(x="")
ggplot (data=d, aes(y=Haemoglobin, x=case)) + geom_boxplot() + labs(x="")
ggplot (data=d, aes(y=Lymphocyte, x=case)) + geom_boxplot() + labs(x="")
ggplot (data=d, aes(y=MCV, x=case)) + geom_boxplot() + labs(x="")
ggplot (data=d, aes(y=Monocyte, x=case)) + geom_boxplot() + labs(x="")
ggplot (data=d, aes(y=Neutrophil, x=case)) + geom_boxplot() + labs(x="")
ggplot (data=d, aes(y=Platelets, x=case)) + geom_boxplot() + labs(x="")
ggplot (data=d, aes(y=WCC, x=case)) + geom_boxplot() + labs(x="")

#plots of continuous variables 
ggplot (data=d, aes(age)) + geom_histogram(aes(y=..density..)) + facet_wrap(~case)
ggplot (data=d, aes(Albumin)) + geom_histogram(aes(y=..density..)) + facet_wrap(~case)
ggplot (data=d, aes(ALP)) + geom_histogram(aes(y=..density..)) + facet_wrap(~case)
ggplot (data=d, aes(ALT)) + geom_histogram(aes(y=..density..)) + facet_wrap(~case)
ggplot (data=d, aes(Basophil)) + geom_histogram(aes(y=..density..)) + facet_wrap(~case)
ggplot (data=d, aes(Calcium)) + geom_histogram(aes(y=..density..)) + facet_wrap(case)
ggplot (data=d, aes(Creatinine)) + geom_histogram(aes(y=..density..)) + facet_wrap(~case)
ggplot (data=d, aes(CRP)) + geom_histogram(aes(y=..density..)) + facet_wrap(~case)
ggplot (data=d, aes(Eosinophil)) + geom_histogram(aes(y=..density..)) + facet_wrap(~case)
ggplot (data=d, aes(Haemoglobin)) + geom_histogram(aes(y=..density..)) + facet_wrap(~case)
ggplot (data=d, aes(Lymphocyte)) + geom_histogram(aes(y=..density..)) + facet_wrap(~case)
ggplot (data=d, aes(MCV)) + geom_histogram(aes(y=..density..)) + facet_wrap(~case)
ggplot (data=d, aes(Monocyte)) + geom_histogram(aes(y=..density..)) + facet_wrap(~case)
ggplot (data=d, aes(Neutrophil)) + geom_histogram(aes(y=..density..)) + facet_wrap(~case)
ggplot (data=d, aes(Platelets)) + geom_histogram(aes(y=..density..)) + facet_wrap(~case)
ggplot (data=d, aes(WCC)) + geom_histogram(aes(y=..density..)) + facet_wrap(~case)

#descriptive stats 
#table of summary stats for cases and controls
Hmisc::summaryM(age + sex + Albumin + ALP + ALT  + Basophil + Calcium +
  Creatinine + CRP + Eosinophil + Haemoglobin + Lymphocyte + MCV +
  Monocyte + Neutrophil + Platelets + WCC ~ MM, 
  data=d, overall=TRUE, continuous=5)

#summary stats for each test by outcome
tapply(d$Albumin, d$MM, summary)
tapply(d$ALP, d$MM, summary)
tapply(d$ALT, d$MM, summary)
tapply(d$Basophil, d$MM, summary)
tapply(d$Creatinine, d$MM, summary)
tapply(d$CRP, d$MM, summary)
tapply(d$Eosinophil, d$MM, summary)
tapply(d$Haemoglobin, d$MM, summary)
tapply(d$Lymphocyte, d$MM, summary)
tapply(d$MCV, d$MM, summary)
tapply(d$Monocyte, d$MM, summary)
tapply(d$Neutrophil, d$MM, summary)
tapply(d$Platelets, d$MM, summary)
tapply(d$WCC, d$MM, summary)

#age and sex
tapply(d$age, d$MM, summary)
table(d$sex, d$MM, useNA="ifany")

#use scale function to centre all continuous variables round mean
d$age_c<-scale(d$age, center=T, scale=F)
d$Albumin_c<-scale(d$Albumin, center=T, scale=F)
d$ALP_c<-scale(d$ALP, center=T, scale=F)
d$ALT_c<-scale(d$ALT, center=T, scale=F)
d$Basophil_c<-scale(d$Basophil, center=T, scale=F)
d$Calcium_c<-scale(d$Calcium, center=T, scale=F)
d$Creatinine_c<-scale(d$Creatinine, center=T, scale=F)
d$CRP_c<-scale(d$CRP, center=T, scale=F)
d$Eosinophil_c<-scale(d$Eosinophil, center=T, scale=F)
d$Haemoglobin_c<-scale(d$Haemoglobin, center=T, scale=F)
d$Lymphocyte_c<-scale(d$Lymphocyte, center=T, scale=F)
d$MCV_c<-scale(d$MCV, center=T, scale=F)
d$Monocyte_c<-scale(d$Monocyte, center=T, scale=F)
d$Neutrophil_c<-scale(d$Neutrophil, center=T, scale=F)
d$Platelets_c<-scale(d$Platelets, center=T, scale=F)
d$WCC_c<-scale(d$WCC, center=T, scale=F)

#univaraible inspection of relationship between predictors and outcomes
ggplot(d, aes(x=age_c, y=MM)) +
  Hmisc::histSpikeg(MM ~ age_c, lowess = TRUE, data=d) +
  labs (x="Age (centred)", y="Prob ")

ggplot(d, aes(x=Albumin_c, y=MM)) +
  Hmisc::histSpikeg(MM ~ Albumin_c, lowess = TRUE, data=d) +
  labs (x="Albumin (centred)", y="Prob ")

ggplot(d, aes(x=ALP_c, y=MM)) +
  Hmisc::histSpikeg(MM ~ ALP_c, lowess = TRUE, data=d) +
  labs (x="ALP (centred)", y="Prob ")

ggplot(d, aes(x=ALT_c, y=MM)) +
  Hmisc::histSpikeg(MM ~ ALT_c, lowess = TRUE, data=d) +
  labs (x="ALT (centred)", y="Prob ")

ggplot(d, aes(x=Basophil_c, y=MM)) +
  Hmisc::histSpikeg(MM ~ Basophil_c, lowess = TRUE, data=d) +
  labs (x="Basophil (centred)", y="Prob ")

ggplot(d, aes(x=Calcium_c, y=MM)) +
  Hmisc::histSpikeg(MM ~ Calcium_c, lowess = TRUE, data=d) +
  labs (x="Calcium (centred)", y="Prob ")

ggplot(d, aes(x=Creatinine_c, y=MM)) +
  Hmisc::histSpikeg(MM ~ Creatinine_c, lowess = TRUE, data=d) +
  labs (x="Creatinine (centred)", y="Prob ")

ggplot(d, aes(x=CRP_c, y=MM)) +
  Hmisc::histSpikeg(MM ~ CRP_c, lowess = TRUE, data=d) +
  labs (x="CRP (centred)", y="Prob ")

ggplot(d, aes(x=Eosinophil_c, y=MM)) +
  Hmisc::histSpikeg(MM ~ Eosinophil_c, lowess = TRUE, data=d) +
  labs (x="Eosinophil (centred)", y="Prob ")

ggplot(d, aes(x=Haemoglobin_c, y=MM)) +
  Hmisc::histSpikeg(MM ~ Haemoglobin_c, lowess = TRUE, data=d) +
  labs (x="Haemoglobin (centred)", y="Prob ")

ggplot(d, aes(x=Lymphocyte_c, y=MM)) +
  Hmisc::histSpikeg(MM ~ Lymphocyte_c, lowess = TRUE, data=d) +
  labs (x="Lymphocyte (centred)", y="Prob ")

ggplot(d, aes(x=MCV_c, y=MM)) +
  Hmisc::histSpikeg(MM ~ MCV_c, lowess = TRUE, data=d) +
  labs (x="MCV (centred)", y="Prob ")

ggplot(d, aes(x=Monocyte_c, y=MM)) +
  Hmisc::histSpikeg(MM ~ Monocyte_c, lowess = TRUE, data=d) +
  labs (x="Monocyte (centred)", y="Prob ")

ggplot(d, aes(x=Neutrophil_c, y=MM)) +
  Hmisc::histSpikeg(MM ~ Neutrophil_c, lowess = TRUE, data=d) +
  labs (x="Neutrophil (centred)", y="Prob ")

ggplot(d, aes(x=Platelets_c, y=MM)) +
  Hmisc::histSpikeg(MM ~ Platelets_c, lowess = TRUE, data=d) +
  labs (x="Platelets (centred)", y="Prob ")

ggplot(d, aes(x=WCC_c, y=MM)) +
  Hmisc::histSpikeg(MM ~ WCC_c, lowess = TRUE, data=d) +
  labs (x="WCC (centred)", y="Prob ")


###########################################
#Models for continuous variables - all using 3 spline knot points

#Albumin, ALP, ALT, Basophil, Calcium, Creatinine, CRP, Eosinophil, 
#Haemoglobin, Lymphocyes, MCV, Monocytes, Neutrophil, Platelets, WCC

#define datadist
dd<-datadist(d)
options(datadist="dd")

#age
agemod <- lrm(MM ~ rcs(age, 3), data=d, x=TRUE, y=TRUE)
print(agemod)
ggplot(Predict(agemod, age, fun=plogis))

#sex
sexmod <- lrm(MM ~ sex, data=d, x=TRUE, y=TRUE)
print(sexmod)
ggplot(Predict(sexmod, sex, fun=plogis))

#albumin
albmod <- lrm(MM ~ rcs(Albumin, 3), data=d, x=TRUE, y=TRUE)
print(alb)
ggplot(Predict(alb, Albumin, fun=plogis))

#ALP
alpmod <- lrm(MM ~ rcs(ALP, 3), data=d, x=TRUE, y=TRUE)
print(alpmod)
ggplot(Predict(alpmod, ALP, fun=plogis))

#ALT
altmod <- lrm(MM ~ rcs(ALT, 3), data=d, x=TRUE, y=TRUE)
print(altmod)
ggplot(Predict(altmod, ALT, fun=plogis))

#Basophil
basmod <- lrm(MM ~ rcs(Basophil, 3), data=d, x=TRUE, y=TRUE)
print(basmod)
ggplot(Predict(basmod, Basophil, fun=plogis))

#calcium
calmod <- lrm(MM ~ rcs(Calcium, 3), data=d, x=TRUE, y=TRUE)
print(calmod)
ggplot(Predict(calmod, Calcium, fun=plogis))

#creatinine
crmod <- lrm(MM ~ rcs(Creatinine, 3), data=d, x=TRUE, y=TRUE)
print(crmod)
ggplot(Predict(crmod, Creatinine, fun=plogis))

#CRP
crpmod <- lrm(MM ~ rcs(CRP, 3), data=d, x=TRUE, y=TRUE)
print(crpmod)
ggplot(Predict(crpmod, CRP, fun=plogis))

#Eosinophil
eosmod <- lrm(MM ~ rcs(Eosinophil, 3), data=d, x=TRUE, y=TRUE)
print(eosmod)
ggplot(Predict(eosmod, Eosinophil, fun=plogis))

#Haemoglobin
hbmod <- lrm(MM ~ rcs(Haemoglobin, 3), data=d, x=TRUE, y=TRUE)
print(hbmod)
ggplot(Predict(hbmod, Haemoglobin, fun=plogis))

#Lymphocyte
lymod <- lrm(MM ~ rcs(Lymphocyte, 3), data=d, x=TRUE, y=TRUE)
print(lymod)
ggplot(Predict(lymod, Lymphocyte, fun=plogis))

#Monocyte
mcvmod <- lrm(MM ~ rcs(MCV, 3), data=d, x=TRUE, y=TRUE)
print(mcvmod)
ggplot(Predict(mcvmod, MCV, fun=plogis))

#Monocyte
monmod <- lrm(MM ~ rcs(Monocyte, 3), data=d, x=TRUE, y=TRUE)
print(monmod)
ggplot(Predict(monmod, Monocyte, fun=plogis))

#Neutrophil
ntmod <- lrm(MM ~ rcs(Neutrophil, 3), data=d, x=TRUE, y=TRUE)
print(ntmod)
ggplot(Predict(ntmod, Neutrophil, fun=plogis))

#Platelets
ptmod <- lrm(MM ~ rcs(Platelets, 3), data=d, x=TRUE, y=TRUE)
print(ptmod)
ggplot(Predict(ptmod, Platelets, fun=plogis))

#WCC
wccmod <- lrm(MM ~ rcs(WCC, 3), data=d, x=TRUE, y=TRUE)
print(wccmod)
ggplot(Predict(wccmod, WCC, fun=plogis))
