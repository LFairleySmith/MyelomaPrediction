#models including MGUS as controls group

library(tidyverse)
library(mice)
library(psfmi)
library(reshape2)
library(miceadds)
library(naniar)
library(viridis)
library(visdat)
library(rms)
library(grid)
library(gridExtra)
library(pROC)
library(dcurves)
library(rmda)
library(ggpubr)

rm(list=ls())
load("N:\\SharedFolder\\DataAnalysis\\CrossSectional\\CasesCX.RData")
load("N:\\SharedFolder\\DataAnalysis\\CrossSectional\\ControlsCX.RData")
load("N:\\SharedFolder\\DataAnalysis\\CrossSectional\\MGUSCX.RData")

d<-rbind(case_cx, con_cx, mgus_cx)
d %>% group_by(case) %>% summarise(n_distinct(id))
table(d$case)
d<- d %>%  mutate(MM = ifelse(case=="Myeloma", 1, 0)) %>% 
  relocate(MM, .after=case)
table(d$MM, d$case)

#box plots over all blood test results 
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

####impute data  
#only keep relevant variables for this analysis - drop id, case 
d<- d %>% select(-id) %>%  rename(Sex=sex, Age=age)
d<- d %>%  mutate(MM = ifelse(case=="Myeloma", 1, 0)) %>% 
  relocate(MM, .after=case)
table(d$MM, d$case)

#Check missing data  
gg_miss_var(d) + theme_bw()
#and impute missing values, 20 imputations
#need to amend predictor matrix to exclude MM and use case in imputation model
d_imp<-mice(d, m=5, maxit=0)
pred<-d_imp$predictorMatrix
pred[, "MM"] <-0
pred 

d_imp<-mice(d, m=20, maxit=1000, predictorMatrix = pred)
d_imp

plot(d_imp)

#get dataframe with imputed values
impD_wMGUS<-complete(d_imp, action="long", include=FALSE)

#save these to use later 
saveRDS(d_imp, file="ImpMIDS_wMGUS")
save(impD_wMGUS, file="Imp_long_wMGUS.RData")

#load for analysis
load("Imp_long_wMGUS.RData")
d_imp<-readRDS(file="ImpMIDS_wMGUS")

#denisty plots for selected variables (or see code below using ggplot)
densityplot(d_imp, ~Albumin + ALP + ALT)
densityplot(d_imp, ~Basophil + Calcium + Creatinine)
densityplot(d_imp, ~CRP + Eosinophil + Haemoglobin)
densityplot(d_imp, ~Lymphocyte + MCV + Monocyte)
densityplot(d_imp, ~Neutrophil + Platelets + WCC)

#Define datadist
dd<-datadist(d)
options(datadist="dd")

#imputation all predictors 
mod1<- fit.mult.impute(MM ~ rcs(Age, 3) + Sex + rcs(Albumin, 3) + rcs(ALP, 3) + rcs(ALT, 3) + rcs(Basophil, 3) + 
                       rcs(Calcium, 3) + rcs(Creatinine, 3) + rcs(CRP, 3) + rcs(Eosinophil, 3) + rcs(Haemoglobin, 3) + 
                       rcs(Lymphocyte, 3) + rcs(MCV, 3) + rcs(Monocyte, 3) + rcs(Neutrophil, 3) + 
                       rcs(Platelets, 3) + rcs(WCC, 3),
                     fitter=lrm, xtrans=d_imp, x=TRUE, y=TRUE)
mod1
specs(mod1)
summary(mod1)

#extract coefficients from model 
mgus_coef<-coef(mod1)
mgus_se<-sqrt(diag(vcov(mod1)))
MGUS<- data.frame(mgus_coef, mgus_se)
write.csv(MGUS, "N:\\SharedFolder\\DataAnalysis\\CrossSectional\\ModelCoeffs\\MOD_MGUS.csv", row.names = TRUE)

#CHeck linearity of variables
plot(anova(mod1))
anova(mod1)

#plots all variables in model 
ggplot(Predict(mod1, fun=plogis), adj.subtitle=FALSE,flipxdiscrete =FALSE )

#using rms validate function
val <- validate(mod1, B=1000)
val
#c index optimum corrected 
0.5*(val[1, 5]+1)

#predicted probabilites - to use in diagnostic stats and calibration plot 
pred_prob <-predict(mod1, type="fitted")
pred_LP <-predict(mod1, type="lp")
mm<- d$MM
summary(pred_prob)
summary(pred_LP)
case<-d$case

c<-data.frame(mm, case, pred_prob, pred_LP) 

#plot of precited probs for each of three groups
ggplot(c, aes(x=pred_prob)) +
  geom_histogram(aes(y=..density..), bins=100) +
  facet_wrap(~case) +
  labs(x = "Predcited probability", 
       y="Denisty %",
       title="Predicted probabity by case status") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limits=c(0, 10), breaks=seq(0, 10, 2.5)) 
  

#calibration plot 
#create 10 risk groups and get observed and expected values
g1<-c %>%  
  mutate(bin=ntile(pred_prob, 10)) %>% 
  group_by(bin) %>% 
  mutate(n = n(), 
         exp = mean(pred_prob), 
         obs = mean(as.numeric(mm)),
         se_obs = sqrt((obs* (1-obs))/n),
         ul_obs = obs + 1.96 * se_obs,
         ll_obs = obs - 1.96 * se_obs) %>% 
  ungroup()

#calibration plot
cp <- ggplot(data=g1, aes(x=exp, y=obs, ymin=ll_obs, ymax=ul_obs)) +
  geom_pointrange(size=0.5, color="black") +
  scale_y_continuous(limits=c(0, 1), breaks = seq(0, 1, by=0.1)) +
  scale_x_continuous(limits=c(0, 1), breaks = seq(0, 1, by=0.1)) +
  geom_abline() +
  geom_smooth(aes(x=pred_prob, y=as.numeric(mm)),
              color="red", se=FALSE, method="loess") +
  xlab("") +
  ylab("Observed Probability") +
  theme_minimal() +
  theme(aspect.ratio = 1) 
cp

#distribution plot
hp<-ggplot(g1, aes(x=pred_prob)) +
  geom_histogram(fill= "black", bins=200) +
  scale_x_continuous(limits=c(0, 1), breaks = seq(0, 1, by=0.1)) +
  xlab("Predcited probability") +
  ylab("") +
  theme_minimal()  +
  theme(panel.grid.minor = element_blank())

hp

g<- arrangeGrob(cp, hp, respect= TRUE, heights = c(1, 0.25), ncol=1)
grid.newpage()
grid.draw(g)

#model diagnostics at varying probability thresholds
roc_mgus<-roc(c$mm, c$pred_prob)
#check thresholds to include 
ci.coords(roc_mgus, x=c(0.1, 0.2, 0.4), input="threshold", ret=c("sensitivity", "specificity"))

#save as dataframe to export results
rocthres_mgus<-tibble(threshold=seq(0.1, 0.5, by=0.05)) %>% 
  mutate(
    data=threshold %>% map(~ {
      res<- roc_mgus %>% ci.coords(x=.x, ret=c("sensitivity", "specificity"))
      list(
        sens_lci = res$sensitivity[[1]],
        sens = res$sensitivity[[2]],
        sens_uci = res$sensitivity[[3]],
        spec_lci = res$specificity[[1]],
        spec = res$specificity[[2]],
        spec_uci = res$specificity[[3]]
      )
    })
  ) %>% 
  unnest_wider(data)
write.csv(rocthres_mgus, "N:\\SharedFolder\\DataAnalysis\\CrossSectional\\SensSpec\\MGUSMOD.csv", row.names = TRUE)

##########################
#to get psfmi to run remove case variable and recode sex as a a binary 
i_wMGUS<- impD_wMGUS %>% select(-case)  %>% 
  mutate(sex = ifelse(Sex=="Female", 1, 0)) %>% 
  select(-Sex)

#model validation using psfmi 
pool_m1<- psfmi_lr(data=i_wMGUS, nimp=20, impvar=".imp", 
                   formula = MM ~ rcs(Age, 3) + sex + rcs(Albumin, 3) + rcs(ALP, 3) + rcs(ALT, 3) + rcs(Basophil, 3) + 
                     rcs(Calcium, 3) + rcs(Creatinine, 3) + rcs(CRP, 3) + rcs(Eosinophil, 3) + rcs(Haemoglobin, 3) + 
                     rcs(Lymphocyte, 3) + rcs(MCV, 3) + rcs(Monocyte, 3) + rcs(Neutrophil, 3) + 
                     rcs(Platelets, 3) + rcs(WCC, 3), method="D1")
m1_perf<- pool_performance(data=i_wMGUS, nimp=20, impvar=".imp", 
                           formula=pool_m1$formula_final,
                           cal.plot=TRUE, plot.method = "mean", 
                           groups_cal = 10)
m1_perf

#internal validation
set.seed(233)
MGUS_val_boot <- psfmi_validate(pool_m1, val_method="MI_boot",
                              data_orig=d, nboot=200,
                              nimp_mice=20, miceImp=miceImp)
MGUS_val_boot$stats_val
MGUS_val_boot$intercept_test
#save output from bootstrap samples to get means and CIs
MGUSboot<-as.data.frame(MGUS_val_boot$res_boot)
save(MGUSboot, file="MGUSboot.RData")
save(MGUS_val_boot, file="MGUSvalbootMI.RData")
load(file="MGUSboot.RData")
load(file="MGUSvalbootMI.RData")

mean(MGUSboot$Slope)
quantile(MGUSboot$Slope, probs=c(0.025, 0.975))

MGUSboot$auc<-0.8517 #apparent AUC
MGUSboot$opt<-MGUSboot$ROC_app-MGUSboot$ROC_test
MGUSboot$auc_cor<-MGUSboot$auc-MGUSboot$opt
mean(MGUSboot$auc_cor)
quantile(MGUSboot$auc_cor, probs=c(0.025, 0.975))


#Decision curve analysis 
dca(mm~pred_prob, data=c, 
          as_probability="pred_prob", 
          label=list(pred_prob= "MGUS model")) %>% 
  plot(smooth=TRUE) +
  ggplot2::labs(x="Threshold probability")


#number of MGUS at different thresholds
d_mgus<- c %>%  filter(case=="MGUS") %>% select(pred_prob)
myseq<-seq(0.1, 0.5, by=0.05)
mgus_pos<-data.frame(threshold = myseq, 
                do.call(rbind, lapply(myseq, function(x) colSums(d_mgus>=x))))

#################################
#Model all ages, MGUS as controls excluding calcium as predictor
mnoc<- fit.mult.impute(MM ~ rcs(Age, 3) + Sex + rcs(Albumin, 3) + rcs(ALP, 3) + rcs(ALT, 3) + rcs(Basophil, 3) + 
                         rcs(Creatinine, 3) + rcs(CRP, 3) + rcs(Eosinophil, 3) + rcs(Haemoglobin, 3) + 
                         rcs(Lymphocyte, 3) + rcs(MCV, 3) + rcs(Monocyte, 3) + rcs(Neutrophil, 3) + 
                         rcs(Platelets, 3) + rcs(WCC, 3),
                       fitter=lrm, xtrans=d_imp, x=TRUE, y=TRUE)
mnoc
specs(mnoc)
summary(mnoc)

#extract coefficients from model 
mgusNoCa_coef<-coef(mnoc)
mgusNoCa_se<-sqrt(diag(vcov(mnoc)))
MGUSNoCa<- data.frame(mgusNoCa_coef, mgusNoCa_se)
write.csv(MGUSNoCa, "N:\\SharedFolder\\DataAnalysis\\CrossSectional\\ModelCoeffs\\MOD_MGUSNoCa.csv", row.names = TRUE)


#CHeck linearity of variables
plot(anova(mnoc))
anova(mnoc)

#plots all variables in model 
ggplot(Predict(mnoc, fun=plogis), adj.subtitle=FALSE, flipxdiscrete =FALSE )

#using rms validate function
val <- validate(mnoc, B=1000)
val
#c index optimum corrected 
0.5*(val[1, 5]+1)

#predicted probabilites - to use in diagnostic stats and calibration plot 
pred_probmnoc <-predict(mnoc, type="fitted")
pred_LPmnoc <-predict(mnoc, type="lp")
mm<- d$MM
summary(pred_probmnoc)
summary(pred_LPmnoc)
case<-d$case

cmnoc<-data.frame(mm, case, pred_probmnoc, pred_LPmnoc) 

#plot of precited probs for each of three groups
ggplot(cmnoc, aes(x=pred_probmnoc)) +
  geom_histogram(aes(y=..density..), bins=100) +
  facet_wrap(~case) +
  labs(x = "Predcited probability", 
       y="Denisty %",
       title="Predicted probabity by case status") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limits=c(0, 10), breaks=seq(0, 10, 2.5)) 


#calibration plot 
#create 10 risk groups and get observed and expected values
g1<-cmnoc %>%  
  mutate(bin=ntile(pred_probmnoc, 10)) %>% 
  group_by(bin) %>% 
  mutate(n = n(), 
         exp = mean(pred_probmnoc), 
         obs = mean(as.numeric(mm)),
         se_obs = sqrt((obs* (1-obs))/n),
         ul_obs = obs + 1.96 * se_obs,
         ll_obs = obs - 1.96 * se_obs) %>% 
  ungroup()

#calibration plot
cp <- ggplot(data=g1, aes(x=exp, y=obs, ymin=ll_obs, ymax=ul_obs)) +
  geom_pointrange(size=0.5, color="black") +
  scale_y_continuous(limits=c(0, 1), breaks = seq(0, 1, by=0.1)) +
  scale_x_continuous(limits=c(0, 1), breaks = seq(0, 1, by=0.1)) +
  geom_abline() +
  geom_smooth(aes(x=pred_probmnoc, y=as.numeric(mm)),
              color="red", se=FALSE, method="loess") +
  xlab("") +
  ylab("Observed Probability") +
  theme_minimal() +
  theme(aspect.ratio = 1) 
cp

#distribution plot
hp<-ggplot(g1, aes(x=pred_probmnoc)) +
  geom_histogram(fill= "black", bins=200) +
  scale_x_continuous(limits=c(0, 1), breaks = seq(0, 1, by=0.1)) +
  xlab("Predcited probability") +
  ylab("") +
  theme_minimal()  +
  theme(panel.grid.minor = element_blank())

hp

g<- arrangeGrob(cp, hp, respect= TRUE, heights = c(1, 0.25), ncol=1)
grid.newpage()
grid.draw(g)

#model diagnostics at varying probability thresholds
roc_mgusnoca<-roc(cmnoc$mm, cmnoc$pred_probmnoc)
#save as dataframe to export results
rocthres_mnoc<-tibble(threshold=seq(0.1, 0.5, by=0.05)) %>% 
  mutate(
    data=threshold %>% map(~ {
      res<- roc_mgusnoca %>% ci.coords(x=.x, ret=c("sensitivity", "specificity"))
      list(
        sens_lci = res$sensitivity[[1]],
        sens = res$sensitivity[[2]],
        sens_uci = res$sensitivity[[3]],
        spec_lci = res$specificity[[1]],
        spec = res$specificity[[2]],
        spec_uci = res$specificity[[3]]
      )
    })
  ) %>% 
  unnest_wider(data)
write.csv(rocthres_mnoc, "N:\\SharedFolder\\DataAnalysis\\CrossSectional\\SensSpec\\MGUSNoCA.csv", row.names = TRUE)


#model validation using psfmi 
pool_mnoc<- psfmi_lr(data=i_wMGUS, nimp=20, impvar=".imp", 
                   formula = MM ~ rcs(Age, 3) + sex + rcs(Albumin, 3) + rcs(ALP, 3) + rcs(ALT, 3) + rcs(Basophil, 3) + 
                     rcs(Creatinine, 3) + rcs(CRP, 3) + rcs(Eosinophil, 3) + rcs(Haemoglobin, 3) + 
                     rcs(Lymphocyte, 3) + rcs(MCV, 3) + rcs(Monocyte, 3) + rcs(Neutrophil, 3) + 
                     rcs(Platelets, 3) + rcs(WCC, 3), method="D1")
m1_perf<- pool_performance(data=i_wMGUS, nimp=20, impvar=".imp", 
                           formula=pool_mnoc$formula_final,
                           cal.plot=TRUE, plot.method = "mean", 
                           groups_cal = 10)
m1_perf

#internal validation
set.seed(233)
MGUS_val_boot <- psfmi_validate(pool_mnoc, val_method="MI_boot",
                                data_orig=d, nboot=200,
                                nimp_mice=20, miceImp=miceImp)
MGUS_val_boot$stats_val
MGUS_val_boot$intercept_test
#save output from bootstrap samples to get means and CIs
MGUSboot<-as.data.frame(MGUS_val_boot$res_boot)
mean(MGUSboot$Slope)
quantile(MGUSboot$Slope, probs=c(0.025, 0.975))

MGUSboot$auc<-0.832 #apparent AUC
MGUSboot$opt<-MGUSboot$ROC_app-MGUSboot$ROC_test
MGUSboot$auc_cor<-MGUSboot$auc-MGUSboot$opt
mean(MGUSboot$auc_cor)
quantile(MGUSboot$auc_cor, probs=c(0.025, 0.975))


####################################
#FURTHER MODELS FOr 60+ only (MGUS as controls)
#model 1 includes all 15 predcitors
#model 2 excludes calcium as a predcitors
####################

#Remove <60s
d <- d %>%
  mutate(age_cat=cut(Age, breaks=c(30, 40, 50, 60, 70, 80, 90, 100), 
                     right=FALSE, 
                     label=c("30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90-99")))
table(d$age_cat, d$case)
table(d$age_cat, d$MM)

#removes under 60s
d60<-d %>% filter(Age>=60)
table(d60$MM)
table(d60$case)
#recode sex as binary for models
d60<- d60 %>% mutate(sex = ifelse(Sex=="Female", 1, 0)) 

#need to define datadist
dd<-datadist(d60)
options(datadist="dd")

#also remove from d_imp
d_imp60<- d_imp %>% filter(Age>=60)
impD60<-complete(d_imp60, action="long", include=FALSE)
summary(impD60$Age)

####MODEL WITH ALL PREDICTORS
m_60<- fit.mult.impute(MM ~ rcs(Age, 3) + Sex + rcs(Albumin, 3) + rcs(ALP, 3) + rcs(ALT, 3) + rcs(Basophil, 3) + 
                         rcs(Calcium, 3) + rcs(Creatinine, 3) + rcs(CRP, 3) + rcs(Eosinophil, 3) + rcs(Haemoglobin, 3) + 
                         rcs(Lymphocyte, 3) + rcs(MCV, 3) + rcs(Monocyte, 3) + rcs(Neutrophil, 3) + 
                         rcs(Platelets, 3) + rcs(WCC, 3),
                       fitter=lrm, xtrans=d_imp60, x=TRUE, y=TRUE)
m_60
specs(m_60)
plot(anova(m_60))

#extract coefficients from model 
mgusO60_coef<-coef(m_60)
mgusO60_se<-sqrt(diag(vcov(m_60)))
MGUSO60<- data.frame(mgusO60_coef, mgusO60_se)
write.csv(MGUSO60, "N:\\SharedFolder\\DataAnalysis\\CrossSectional\\ModelCoeffs\\MOD_MGUSOver60.csv", row.names = TRUE)


#plots all variables in model 
ggplot(Predict(m_60, fun=plogis), adj.subtitle=FALSE, flipxdiscrete =FALSE )

#using rms validate function
val <- validate(m_60, B=1000)
val
#c index optimum corrected 
0.5*(val[1, 5]+1)

#model validation using psfmi 
#to get psfmi to run remove case variable and recode sex as a a binary 
i_60<- impD60 %>% select(-case)  %>% 
  mutate(sex = ifelse(Sex=="Female", 1, 0)) %>% 
  select(-Sex)

#use psfmi to get model estimates 
pool_m60<- psfmi_lr(data=i_60, nimp=20, impvar=".imp", 
                    formula = MM ~ rcs(Age, 3) + sex + rcs(Albumin, 3) + rcs(ALP, 3) + rcs(ALT, 3) + rcs(Basophil, 3) + 
                      rcs(Calcium, 3) + rcs(Creatinine, 3) + rcs(CRP, 3) + rcs(Eosinophil, 3) + rcs(Haemoglobin, 3) + 
                      rcs(Lymphocyte, 3) + rcs(MCV, 3) + rcs(Monocyte, 3) + rcs(Neutrophil, 3) + 
                      rcs(Platelets, 3) + rcs(WCC, 3), method="D1")
pool_m60$RR_model
#pooled performance stats
m60_perf<- pool_performance(data=i_60, nimp=20, impvar=".imp", 
                            formula=pool_m60$formula_final,
                            cal.plot=TRUE, plot.method = "mean", 
                            groups_cal = 10)
m60_perf

#internal validation
set.seed(233)
rm60_MI_boot <- psfmi_validate(pool_m60, val_method="MI_boot",
                               data_orig=d60, nboot=200,
                               nimp_mice=20, miceImp=miceImp)
rm60_MI_boot$stats_val
rm60_MI_boot$intercept_test
Mboot<-as.data.frame(rm60_MI_boot$res_boot)
mean(Mboot$Slope)
quantile(Mboot$Slope, probs=c(0.025, 0.975))

Mboot$auc<-0.864
Mboot$opt<-Mboot$ROC_app-Mboot$ROC_test
Mboot$auc_cor<-Mboot$auc-Mboot$opt
mean(Mboot$auc_cor)
quantile(Mboot$auc_cor, probs=c(0.025, 0.975))

#predicted probabilites - to use in diagnostic stats and calibration plot (add to dataset C with LP and PP from full model)
pred_prob60 <-predict(m_60, type="fitted")
pred_LP60 <-predict(m_60, type="lp")
summary(pred_prob60)
summary(pred_LP60)
mm<- d60$MM

c60<-data.frame(mm, pred_prob60, pred_LP60) 

#calibration plot 
#create 10 risk groups and get observed and expected values
g1<-c60 %>%  
  mutate(bin=ntile(pred_prob60, 10)) %>% 
  group_by(bin) %>% 
  mutate(n = n(), 
         exp = mean(pred_prob60), 
         obs = mean(as.numeric(mm)),
         se_obs = sqrt((obs* (1-obs))/n),
         ul_obs = obs + 1.96 * se_obs,
         ll_obs = obs - 1.96 * se_obs) %>% 
  ungroup()

#calibration plot
cp <- ggplot(data=g1, aes(x=exp, y=obs, ymin=ll_obs, ymax=ul_obs)) +
  geom_pointrange(size=0.5, color="black") +
  scale_y_continuous(limits=c(0, 1), breaks = seq(0, 1, by=0.1)) +
  scale_x_continuous(limits=c(0, 1), breaks = seq(0, 1, by=0.1)) +
  geom_abline() +
  geom_smooth(aes(x=pred_prob60, y=as.numeric(mm)),
              color="red", se=FALSE, method="loess") +
  xlab("") +
  ylab("Observed Probability") +
  theme_minimal() +
  theme(aspect.ratio = 1) 
cp

#distribution plot
hp<-ggplot(g1, aes(x=pred_prob60)) +
  geom_histogram(fill= "black", bins=200) +
  scale_x_continuous(limits=c(0, 1), breaks = seq(0, 1, by=0.1)) +
  xlab("Predcited probability") +
  ylab("") +
  theme_minimal()  +
  theme(panel.grid.minor = element_blank())

hp

g<- arrangeGrob(cp, hp, respect= TRUE, heights = c(1, 0.25), ncol=1)
grid.newpage()
grid.draw(g)

#back to back histogram - counts
hp2<-ggplot(g1, aes(x=pred_prob60)) +
  geom_histogram(data=subset(g1, mm==0), 
                 aes(y=..count.., fill="Control"), bins=100) +
  geom_histogram(data=subset(g1, mm==1), 
                 aes(y=- ..count.., fill="MM"), bins=100) +
  labs(x = "Predcited probability", 
       y="", 
       fill="") +
  theme_minimal() +
  scale_y_continuous(limits=c(-25, 150), breaks= c(-25, 0, 50, 100, 150)) +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  scale_fill_manual(values=c("grey", "black")) +
  theme(legend.position = c(0.85, 0.8), 
        legend.background=element_rect(fill="white", color="black"))  
hp2    


#model diagnostics at varying probability thresholds (0.1, 0.2. 0.4)
roc_mgus60<-roc(c60$mm, c60$pred_prob60)
roc_mgus60
#save as dataframe to export results
rocthres_mgus60<-tibble(threshold=seq(0.1, 0.5, by=0.05)) %>% 
  mutate(
    data=threshold %>% map(~ {
      res<- roc_mgus60 %>% ci.coords(x=.x, ret=c("sensitivity", "specificity"))
      list(
        sens_lci = res$sensitivity[[1]],
        sens = res$sensitivity[[2]],
        sens_uci = res$sensitivity[[3]],
        spec_lci = res$specificity[[1]],
        spec = res$specificity[[2]],
        spec_uci = res$specificity[[3]]
      )
    })
  ) %>% 
  unnest_wider(data)
write.csv(rocthres_mgus60, "N:\\SharedFolder\\DataAnalysis\\CrossSectional\\SensSpec\\MGUS60.csv", row.names = TRUE)

####MODEL EXCLUDING CALCIUM
mnoc_60<- fit.mult.impute(MM ~ rcs(Age, 3) + Sex + rcs(Albumin, 3) + rcs(ALP, 3) + rcs(ALT, 3) + rcs(Basophil, 3) + 
                         rcs(Creatinine, 3) + rcs(CRP, 3) + rcs(Eosinophil, 3) + rcs(Haemoglobin, 3) + 
                         rcs(Lymphocyte, 3) + rcs(MCV, 3) + rcs(Monocyte, 3) + rcs(Neutrophil, 3) + 
                         rcs(Platelets, 3) + rcs(WCC, 3),
                       fitter=lrm, xtrans=d_imp60, x=TRUE, y=TRUE)
mnoc_60
specs(mnoc_60)
plot(anova(mnoc_60))

#extract coefficients from model 
mgusO60NoCa_coef<-coef(mnoc_60)
mgusO60NoCa_se<-sqrt(diag(vcov(mnoc_60)))
MGUSO60NoCa<- data.frame(mgusO60NoCa_coef, mgusO60NoCa_se)
write.csv(MGUSO60NoCa, "N:\\SharedFolder\\DataAnalysis\\CrossSectional\\ModelCoeffs\\MOD_MGUSOver60NoCa.csv", row.names = TRUE)

#plots all variables in model 
ggplot(Predict(mnoc_60, fun=plogis), adj.subtitle=FALSE, flipxdiscrete =FALSE )

#using rms validate function
val <- validate(mnoc_60, B=1000)
val
#c index optimum corrected 
0.5*(val[1, 5]+1)

#use psfmi to get model estimates 
pool_m60<- psfmi_lr(data=i_60, nimp=20, impvar=".imp", 
                    formula = MM ~ rcs(Age, 3) + sex + rcs(Albumin, 3) + rcs(ALP, 3) + rcs(ALT, 3) + rcs(Basophil, 3) + 
                      rcs(Creatinine, 3) + rcs(CRP, 3) + rcs(Eosinophil, 3) + rcs(Haemoglobin, 3) + 
                      rcs(Lymphocyte, 3) + rcs(MCV, 3) + rcs(Monocyte, 3) + rcs(Neutrophil, 3) + 
                      rcs(Platelets, 3) + rcs(WCC, 3), method="D1")
pool_m60$RR_model
#pooled performance stats
m60_perf<- pool_performance(data=i_60, nimp=20, impvar=".imp", 
                            formula=pool_m60$formula_final,
                            cal.plot=TRUE, plot.method = "mean", 
                            groups_cal = 10)
m60_perf

#internal validation
set.seed(233)
rmnoc60_MI_boot <- psfmi_validate(pool_m60, val_method="MI_boot",
                               data_orig=d60, nboot=200,
                               nimp_mice=20, miceImp=miceImp)
rmnoc60_MI_boot$stats_val
rmnoc60_MI_boot$intercept_test
Mboot<-as.data.frame(rmnoc60_MI_boot$res_boot)
mean(Mboot$Slope)
quantile(Mboot$Slope, probs=c(0.025, 0.975))

Mboot$auc<-0.846
Mboot$opt<-Mboot$ROC_app-Mboot$ROC_test
Mboot$auc_cor<-Mboot$auc-Mboot$opt
mean(Mboot$auc_cor)
quantile(Mboot$auc_cor, probs=c(0.025, 0.975))

#predicted probabilites - to use in diagnostic stats and calibration plot (add to dataset C with LP and PP from full model)
pred_prob60noc <-predict(mnoc_60, type="fitted")
pred_LP60noc <-predict(mnoc_60, type="lp")
summary(pred_prob60noc)
summary(pred_LP60noc)
mm<- d60$MM

c60noc<-data.frame(mm, pred_prob60noc, pred_LP60noc) 

#calibration plot 
#create 10 risk groups and get observed and expected values
g1<-c60noc %>%  
  mutate(bin=ntile(pred_prob60noc, 10)) %>% 
  group_by(bin) %>% 
  mutate(n = n(), 
         exp = mean(pred_prob60noc), 
         obs = mean(as.numeric(mm)),
         se_obs = sqrt((obs* (1-obs))/n),
         ul_obs = obs + 1.96 * se_obs,
         ll_obs = obs - 1.96 * se_obs) %>% 
  ungroup()

#calibration plot
cp <- ggplot(data=g1, aes(x=exp, y=obs, ymin=ll_obs, ymax=ul_obs)) +
  geom_pointrange(size=0.5, color="black") +
  scale_y_continuous(limits=c(0, 1), breaks = seq(0, 1, by=0.1)) +
  scale_x_continuous(limits=c(0, 1), breaks = seq(0, 1, by=0.1)) +
  geom_abline() +
  geom_smooth(aes(x=pred_prob60noc, y=as.numeric(mm)),
              color="red", se=FALSE, method="loess") +
  xlab("") +
  ylab("Observed Probability") +
  theme_minimal() +
  theme(aspect.ratio = 1) 
cp

#distribution plot
hp<-ggplot(g1, aes(x=pred_prob60noc)) +
  geom_histogram(fill= "black", bins=200) +
  scale_x_continuous(limits=c(0, 1), breaks = seq(0, 1, by=0.1)) +
  xlab("Predcited probability") +
  ylab("") +
  theme_minimal()  +
  theme(panel.grid.minor = element_blank())
hp

g<- arrangeGrob(cp, hp, respect= TRUE, heights = c(1, 0.25), ncol=1)
grid.newpage()
grid.draw(g)

#back to back histogram - counts
hp2<-ggplot(g1, aes(x=pred_prob60noc)) +
  geom_histogram(data=subset(g1, mm==0), 
                 aes(y=..count.., fill="Control"), bins=100) +
  geom_histogram(data=subset(g1, mm==1), 
                 aes(y=- ..count.., fill="MM"), bins=100) +
  labs(x = "Predcited probability", 
       y="", 
       fill="") +
  theme_minimal() +
  scale_y_continuous(limits=c(-25, 150), breaks= c(-25, 0, 50, 100, 150)) +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  scale_fill_manual(values=c("grey", "black")) +
  theme(legend.position = c(0.85, 0.8), 
        legend.background=element_rect(fill="white", color="black"))  
hp2    


#model diagnostics at varying probability thresholds (0.1, 0.2. 0.4)
roc_mgus60noCa<-roc(c60noc$mm, c60noc$pred_prob60noc)
roc_mgus60noCa
#save as dataframe to export results
rocthres_m60noc<-tibble(threshold=seq(0.1, 0.5, by=0.05)) %>% 
  mutate(
    data=threshold %>% map(~ {
      res<- roc_mgus60noCa %>% ci.coords(x=.x, ret=c("sensitivity", "specificity"))
      list(
        sens_lci = res$sensitivity[[1]],
        sens = res$sensitivity[[2]],
        sens_uci = res$sensitivity[[3]],
        spec_lci = res$specificity[[1]],
        spec = res$specificity[[2]],
        spec_uci = res$specificity[[3]]
      )
    })
  ) %>% 
  unnest_wider(data)
write.csv(rocthres_m60noc, "N:\\SharedFolder\\DataAnalysis\\CrossSectional\\SensSpec\\MGUS60NoCA.csv", row.names = TRUE)

