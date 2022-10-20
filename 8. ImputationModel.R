#Mutiple imputation using mice
#run model using psfmi
#both full model and using backwards selection

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
setwd("N:\\SharedFolder\\DataAnalysis\\CrossSectional\\Imputation")
load("N:\\SharedFolder\\DataAnalysis\\CrossSectional\\CasesCX.RData")
load("N:\\SharedFolder\\DataAnalysis\\CrossSectional\\ControlsCX.RData")

d<-rbind(case_cx, con_cx)
d$case = droplevels(d$case)
d %>% group_by(case) %>% summarise(n_distinct(id))
table(d$case)
#recode outcomes 
d<- d %>%  mutate(MM = ifelse(case=="Myeloma", 1, 0)) %>% 
  relocate(MM, .after=case)
table(d$MM)
colnames(d)
#only keep relevant variables for this analysis - drop id, case 
d<- d %>% select(-c(id, case)) %>%  rename(Sex=sex, Age=age) %>% select(MM, Age, Sex, everything())
#define datadist
dd<-datadist(d)
options(datadist="dd")


#Check missing data  
vis_dat(d)
gg_miss_var(d) + theme_bw()

#and impute missing values, 20 imputations
d_imp<-mice(d, m=20, maxit=1000)
plot(d_imp)
d_imp
d_imp$method

#get dataframe with imputed values
impD<-complete(d_imp, action="long", include=FALSE)

#save these to use later 
saveRDS(d_imp, file="ImpMIDS")
save(impD, file="Imp_long.RData")

#load for analysis
load("Imp_long.RData")
d_imp<-readRDS(file="ImpMIDS")

#denisty plots for blood test variables 
densityplot(d_imp, ~Albumin + ALP + ALT)
densityplot(d_imp, ~Basophil + Calcium + Creatinine)
densityplot(d_imp, ~CRP + Eosinophil + Haemoglobin)
densityplot(d_imp, ~Lymphocyte + MCV + Monocyte)
densityplot(d_imp, ~Neutrophil + Platelets + WCC)
 
#compare distributions of observed and imputed data across variables using ggplot
all_D <- complete(d_imp, action= "long", include=TRUE)
imp<- melt(all_D, c(".imp", ".id"))
imp$Imputed <- ifelse(imp$".imp"==0, "Observed", "Imputed")
imp<- imp %>% filter(! (variable %in% c("MM", "sex", "age")))
ggplot(data=imp, aes(x=value, group=.imp, color=Imputed)) +
  stat_density(geom="path", position="identity") +
  facet_wrap(~variable, ncol=3, scales="free")

###### Modeling with imputed data
######
#imputation using fit.mult.impute from rms package 
m2<- fit.mult.impute(MM ~ rcs(Age, 3) + Sex + rcs(Albumin, 3) + rcs(ALP, 3) + rcs(ALT, 3) + rcs(Basophil, 3) + 
                       rcs(Calcium, 3) + rcs(Creatinine, 3) + rcs(CRP, 3) + rcs(Eosinophil, 3) + rcs(Haemoglobin, 3) + 
                       rcs(Lymphocyte, 3) + rcs(MCV, 3) + rcs(Monocyte, 3) + rcs(Neutrophil, 3) + 
                       rcs(Platelets, 3) + rcs(WCC, 3),
                     fitter=lrm, xtrans=d_imp, x=TRUE, y=TRUE)
m2
specs(m2)
round(exp(m2$coef), 3)
round(exp(confint.default(m2)), 3)

#extract coefficients from model 
FM_coef<-coef(m2)
FM_se<-sqrt(diag(vcov(m2)))
FM<- data.frame(FM_coef, FM_se)
write.csv(FM, "N:\\SharedFolder\\DataAnalysis\\CrossSectional\\ModelCoeffs\\FullMOD.csv", row.names = TRUE)

# save model 
saveRDS(m2, file="FullModel")

#Check linearity of variables
plot(anova(m2))
anova(m2)

#plots all variables in model 
ggplot(Predict(m2, fun=plogis), adj.subtitle=FALSE,flipxdiscrete =FALSE )
#separate plots and then combine - to use in paper
p1 <- ggplot(Predict(m2, Albumin, fun=plogis), adj.subtitle=FALSE, ylab="Predicted probability", ylim=c(0, 1)) +
  theme_bw(base_size=7 )
p2 <- ggplot(Predict(m2, ALP, fun=plogis), adj.subtitle=FALSE, ylab="Predicted probability", ylim=c(0, 1)) +
  theme_bw(base_size=7 )
p3 <- ggplot(Predict(m2, ALT, fun=plogis), adj.subtitle=FALSE, ylab="Predicted probability", ylim=c(0, 1)) +
  theme_bw(base_size=7 )
p4 <- ggplot(Predict(m2, Basophil, fun=plogis), adj.subtitle=FALSE, ylab="Predicted probability", ylim=c(0, 1)) +
  theme_bw(base_size=7 )
p5 <- ggplot(Predict(m2, Calcium, fun=plogis), adj.subtitle=FALSE, ylab="Predicted probability", ylim=c(0, 1)) +
  theme_bw(base_size=7 )
p6 <- ggplot(Predict(m2, Creatinine, fun=plogis), adj.subtitle=FALSE, ylab="Predicted probability", ylim=c(0, 1)) +
  theme_bw(base_size=7 )
p7 <- ggplot(Predict(m2, CRP, fun=plogis), adj.subtitle=FALSE, ylab="Predicted probability", ylim=c(0, 1)) +
  theme_bw(base_size=7 )
p8 <- ggplot(Predict(m2, Eosinophil, fun=plogis), adj.subtitle=FALSE, ylab="Predicted probability", ylim=c(0, 1)) +
  theme_bw(base_size=7 )
p9 <- ggplot(Predict(m2, Haemoglobin, fun=plogis), adj.subtitle=FALSE, ylab="Predicted probability", ylim=c(0, 1)) +
  theme_bw(base_size=7 )
p10 <- ggplot(Predict(m2, Lymphocyte, fun=plogis), adj.subtitle=FALSE, ylab="Predicted probability", ylim=c(0, 1)) +
  theme_bw(base_size=7 )
p11 <- ggplot(Predict(m2, MCV, fun=plogis), adj.subtitle=FALSE, ylab="Predicted probability", ylim=c(0, 1)) +
  theme_bw(base_size=7 )
p12 <- ggplot(Predict(m2, Monocyte, fun=plogis), adj.subtitle=FALSE, ylab="Predicted probability", ylim=c(0, 1)) +
  theme_bw(base_size=7 )
p13 <- ggplot(Predict(m2, Neutrophil, fun=plogis), adj.subtitle=FALSE, ylab="Predicted probability", ylim=c(0, 1), xlim=c(2, 15)) +
  theme_bw(base_size=7 )
p14 <- ggplot(Predict(m2, Platelets, fun=plogis), adj.subtitle=FALSE, ylab="Predicted probability", ylim=c(0, 1)) +
  theme_bw(base_size=7 )
p15 <- ggplot(Predict(m2, WCC, fun=plogis), adj.subtitle=FALSE, ylab="Predicted probability", ylim=c(0, 1), xlim=c(2, 15)) +
  theme_bw(base_size=7 )
cowplot::plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, 
                   nrow=5, ncol=3)
pp_fm<- cowplot::plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, 
                      nrow=5, ncol=3)
pp_fm
#save plot specifying height and width dimensions
ggsave("N:\\SharedFolder\\DataAnalysis\\CrossSectional\\PredProbBT.tiff", plot=pp_fm, dpi=300, width=8.2, height=6.9, unit="in")

#predicted probabilites - to use in diagnostic stats and calibration plot 
pred_prob <-predict(m2, type="fitted")
pred_LP <-predict(m2, type="lp")
mm<- d$MM
summary(pred_prob)
summary(pred_LP)

c<-data.frame(mm, pred_prob, pred_LP) 

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

#calibration plot - black dashed line
cp <- ggplot(data=g1, aes(x=exp, y=obs, ymin=ll_obs, ymax=ul_obs)) +
  geom_pointrange(size=0.5, color="black") +
  scale_y_continuous(limits=c(0, 1), breaks = seq(0, 1, by=0.1)) +
  scale_x_continuous(limits=c(0, 1), breaks = seq(0, 1, by=0.1)) +
  geom_abline() +
  geom_smooth(aes(x=pred_prob, y=as.numeric(mm)),
              color="black", se=FALSE, method="loess", linetype="dashed") +
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
ggsave("N:\\SharedFolder\\DataAnalysis\\CrossSectional\\CalPlot_black.tiff", plot=g, dpi=300)

#back to back histogram - counts
hp2<-ggplot(g1, aes(x=pred_prob)) +
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
ggsave("N:\\SharedFolder\\DataAnalysis\\CrossSectional\\PredProbcaseControl.tiff", plot=hp2, dpi=300)


#model diagnostics at varying probability thresholds
roc_imp<-roc(c$mm, c$pred_prob)
#save as dataframe to export results
rocthres_imp<-tibble(threshold=seq(0.1, 0.5, by=0.05)) %>% 
  mutate(
    data=threshold %>% map(~ {
      res<- roc_imp %>% ci.coords(x=.x, ret=c("sensitivity", "specificity"))
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
write.csv(rocthres_imp, "N:\\SharedFolder\\DataAnalysis\\CrossSectional\\SensSpec\\FullMOD.csv", row.names = TRUE)


#model validation using psfmi 
#use psfmi to get model estimates 
pool_m1<- psfmi_lr(data=impD, nimp=20, impvar=".imp", 
                   formula = MM ~ rcs(Age, 3) + Sex + rcs(Albumin, 3) + rcs(ALP, 3) + rcs(ALT, 3) + rcs(Basophil, 3) + 
                     rcs(Calcium, 3) + rcs(Creatinine, 3) + rcs(CRP, 3) + rcs(Eosinophil, 3) + rcs(Haemoglobin, 3) + 
                     rcs(Lymphocyte, 3) + rcs(MCV, 3) + rcs(Monocyte, 3) + rcs(Neutrophil, 3) + 
                     rcs(Platelets, 3) + rcs(WCC, 3), method="D1")
pool_m1$RR_model
#save model coefficients as dataframe (to export for results tables)
m1_coeff<-pool_m1$RR_model
m1_coeff<- m1_coeff %>% map_df(as_tibble)
write.csv(m1_coeff, "N:\\SharedFolder\\DataAnalysis\\CrossSectional\\FullIMPmod_coeff.csv", row.names = TRUE)
pool_m1$formula_final

m1_perf<- pool_performance(data=impD, nimp=20, impvar=".imp", 
                           formula=pool_m1$formula_final,
                           cal.plot=TRUE, plot.method = "mean", 
                           groups_cal = 10)
m1_perf

#internal validation
set.seed(233)
res_MI_boot <- psfmi_validate(pool_m1, val_method="MI_boot",
                              data_orig=d, nboot=1000,
                         nimp_mice=20, miceImp=miceImp)
res_MI_boot$stats_val
res_MI_boot$intercept_test
res_MI_boot
#save output from bootstrap samples to get means and CIs
M1boot<-as.data.frame(res_MI_boot$res_boot)
save(M1boot, file="M1boot.RData")
save(res_MI_boot, file="Mi_bootM1.RData")
load(file="Mi_bootM1.RData")
load(file="M1boot.RData")

mean(M1boot$Slope)
quantile(M1boot$Slope, probs=c(0.025, 0.975))

M1boot$auc<-0.8699
M1boot$opt<-M1boot$ROC_app-M1boot$ROC_test
M1boot$auc_cor<-M1boot$auc-M1boot$opt
mean(M1boot$auc_cor)
quantile(M1boot$auc_cor, probs=c(0.025, 0.975))


#using rms validate function
val <- validate(m2, B=1000)
val
#c index optimum corrected 
0.5*(val[1, 5]+1)

#Decision curve analysis 
dcares<-decision_curve(mm~pred_prob, family=binomial(link="logit"),
                       policy=c("opt-in", "opt-out"), fitted.risk=TRUE,
                       thresholds=seq(0, 1, by=0.01), confidence.intervals = 0.95,
                       bootstraps = 1000, data=c)
plot_decision_curve(dcares, confidence.intervals = TRUE,
                    legend.position = c("topright"), cost.benefit.axis = FALSE,
                    col=c("blue"))
summary(dcares, measure=c("NB"))

dcap<-dca(mm~pred_prob, data=c, 
    as_probability="pred_prob", 
    label=list(pred_prob= "Full model")) %>% 
  plot(smooth=TRUE) +
  ggplot2::labs(x="Threshold probability")
dcap
ggsave("N:\\SharedFolder\\DataAnalysis\\CrossSectional\\DCAplot.tiff", plot=dcap, dpi=300)


##########################################################
#MODEL excluding calcium as predictor
m_noc<- fit.mult.impute(MM ~ rcs(Age, 3) + Sex + rcs(Albumin, 3) + rcs(ALP, 3) + rcs(ALT, 3) + rcs(Basophil, 3) + 
                       rcs(Creatinine, 3) + rcs(CRP, 3) + rcs(Eosinophil, 3) + rcs(Haemoglobin, 3) + 
                       rcs(Lymphocyte, 3) + rcs(MCV, 3) + rcs(Monocyte, 3) + rcs(Neutrophil, 3) + 
                       rcs(Platelets, 3) + rcs(WCC, 3),
                     fitter=lrm, xtrans=d_imp, x=TRUE, y=TRUE)
m_noc
specs(m_noc)
plot(anova(m_noc))

#extract coefficients from model 
NoCa_coef<-coef(m_noc)
NoCa_se<-sqrt(diag(vcov(m_noc)))
NoCa<- data.frame(NoCa_coef, NoCa_se)
write.csv(NoCa, "N:\\SharedFolder\\DataAnalysis\\CrossSectional\\ModelCoeffs\\MODNoCa.csv", row.names = TRUE)

#plots all variables in model 
ggplot(Predict(m_noc, fun=plogis), adj.subtitle=FALSE, flipxdiscrete =FALSE )

#using rms validate function
val <- validate(m_noc, B=1000)
val
#c index optimum corrected 
0.5*(val[1, 5]+1)

#model validation using psfmi 
#use psfmi to get model estimates 
pool_m_noc<- psfmi_lr(data=impD, nimp=20, impvar=".imp", 
                   formula = MM ~ rcs(Age, 3) + Sex + rcs(Albumin, 3) + rcs(ALP, 3) + rcs(ALT, 3) + rcs(Basophil, 3) + 
                     rcs(Creatinine, 3) + rcs(CRP, 3) + rcs(Eosinophil, 3) + rcs(Haemoglobin, 3) + 
                     rcs(Lymphocyte, 3) + rcs(MCV, 3) + rcs(Monocyte, 3) + rcs(Neutrophil, 3) + 
                     rcs(Platelets, 3) + rcs(WCC, 3), method="D1")
pool_m_noc$RR_model
#pooled performance stats
m1_perf<- pool_performance(data=impD, nimp=20, impvar=".imp", 
                           formula=pool_m_noc$formula_final,
                           cal.plot=TRUE, plot.method = "mean", 
                           groups_cal = 10)
m1_perf

#internal validation
set.seed(233)
res_MI_boot <- psfmi_validate(pool_m_noc, val_method="MI_boot",
                              data_orig=d, nboot=200,
                              nimp_mice=20, miceImp=miceImp)
res_MI_boot$stats_val
res_MI_boot$intercept_test
Mboot<-as.data.frame(res_MI_boot$res_boot)
mean(Mboot$Slope)
quantile(Mboot$Slope, probs=c(0.025, 0.975))

Mboot$auc<-0.8467
Mboot$opt<-Mboot$ROC_app-Mboot$ROC_test
Mboot$auc_cor<-Mboot$auc-Mboot$opt
mean(Mboot$auc_cor)
quantile(Mboot$auc_cor, probs=c(0.025, 0.975))

#predicted probabilites - to use in diagnostic stats and calibration plot (add to dataset C with LP and PP from full model)
c$pred_prob2 <-predict(m_noc, type="fitted")
c$pred_LP2 <-predict(m_noc, type="lp")
summary(c$pred_prob2)
summary(c$pred_LP2)

#calibration plot 
#create 10 risk groups and get observed and expected values
g1<-c %>%  
  mutate(bin=ntile(pred_prob2, 10)) %>% 
  group_by(bin) %>% 
  mutate(n = n(), 
         exp = mean(pred_prob2), 
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

#back to back histogram - counts
hp2<-ggplot(g1, aes(x=pred_prob)) +
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
roc_noC<-roc(c$mm, c$pred_prob2)
roc_noC
#save as dataframe to export results
rocthres_noC<-tibble(threshold=seq(0.1, 0.5, by=0.05)) %>% 
  mutate(
    data=threshold %>% map(~ {
      res<- roc_noC %>% ci.coords(x=.x, ret=c("sensitivity", "specificity"))
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
write.csv(rocthres_noC, "N:\\SharedFolder\\DataAnalysis\\CrossSectional\\SensSpec\\NoCaMOD.csv", row.names = TRUE)


#Decision curve analysis 
dcares<-decision_curve(mm~pred_prob, family=binomial(link="logit"),
                       policy=c("opt-in", "opt-out"), fitted.risk=TRUE,
                       thresholds=seq(0, 1, by=0.01), confidence.intervals = 0.95,
                       bootstraps = 1000, data=c)
plot_decision_curve(dcares, confidence.intervals = TRUE,
                    legend.position = c("topright"), cost.benefit.axis = FALSE,
                    col=c("blue"))
summary(dcares, measure=c("NB"))
#Decision curve 
dca(mm~ pred_prob2, data=c2, 
    label=list(pred_prob_bw = "No Ca model")) %>% 
  plot(smooth=TRUE) +
  ggplot2::labs(x="Threshold probability")

#DCA plot with both models - full model and no calcium model 
dca(mm~pred_prob + pred_prob2, data=c, 
    label=list(pred_prob= "Full model", pred_prob2 = "No Ca model")) %>% 
  plot(smooth=TRUE) +
  ggplot2::labs(x="Threshold probability")

####################
#MODEL excluding <60s
d <- d %>%
  mutate(age_cat=cut(Age, breaks=c(30, 40, 50, 60, 70, 80, 90, 100), 
                     right=FALSE, 
                     label=c("30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90-99")))
table(d$age_cat, d$MM)
#removes under 60s
d60<-d %>% filter(Age>=60)
table(d60$MM)

#need to define datadist
dd<-datadist(d60)
options(datadist="dd")

#also remove from d_imp
d_imp60<- d_imp %>% filter(Age>=60)
impD60<-complete(d_imp60, action="long", include=FALSE)
summary(impD60$Age)

m_60<- fit.mult.impute(MM ~ rcs(Age, 3) + Sex + rcs(Albumin, 3) + rcs(ALP, 3) + rcs(ALT, 3) + rcs(Basophil, 3) + 
                          rcs(Calcium, 3) + rcs(Creatinine, 3) + rcs(CRP, 3) + rcs(Eosinophil, 3) + rcs(Haemoglobin, 3) + 
                          rcs(Lymphocyte, 3) + rcs(MCV, 3) + rcs(Monocyte, 3) + rcs(Neutrophil, 3) + 
                          rcs(Platelets, 3) + rcs(WCC, 3),
                        fitter=lrm, xtrans=d_imp60, x=TRUE, y=TRUE)
m_60
specs(m_60)
plot(anova(m_60))

#extract coefficients from model 
O60_coef<-coef(m_60)
O60_se<-sqrt(diag(vcov(m_60)))
O60<- data.frame(O60_coef, O60_se)
write.csv(O60, "N:\\SharedFolder\\DataAnalysis\\CrossSectional\\ModelCoeffs\\MODOver60.csv", row.names = TRUE)


#plots all variables in model 
ggplot(Predict(m_60, fun=plogis), adj.subtitle=FALSE, flipxdiscrete =FALSE )

#using rms validate function
val <- validate(m_60, B=1000)
val
#c index optimum corrected 
0.5*(val[1, 5]+1)

#model validation using psfmi 
#use psfmi to get model estimates 
pool_m60<- psfmi_lr(data=impD60, nimp=20, impvar=".imp", 
                      formula = MM ~ rcs(Age, 3) + Sex + rcs(Albumin, 3) + rcs(ALP, 3) + rcs(ALT, 3) + rcs(Basophil, 3) + 
                        rcs(Calcium, 3) + rcs(Creatinine, 3) + rcs(CRP, 3) + rcs(Eosinophil, 3) + rcs(Haemoglobin, 3) + 
                        rcs(Lymphocyte, 3) + rcs(MCV, 3) + rcs(Monocyte, 3) + rcs(Neutrophil, 3) + 
                        rcs(Platelets, 3) + rcs(WCC, 3), method="D1")
pool_m60$RR_model
#pooled performance stats
m60_perf<- pool_performance(data=impD60, nimp=20, impvar=".imp", 
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

Mboot$auc<-0.884
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
roc_60<-roc(c60$mm, c60$pred_prob60)
roc_60
#save as dataframe to export results
rocthres_60<-tibble(threshold=seq(0.1, 0.5, by=0.05)) %>% 
  mutate(
    data=threshold %>% map(~ {
      res<- roc_60 %>% ci.coords(x=.x, ret=c("sensitivity", "specificity"))
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
write.csv(rocthres_60, "N:\\SharedFolder\\DataAnalysis\\CrossSectional\\SensSpec\\OVER60MOD.csv", row.names = TRUE)


#Decision curve analysis 
dca(mm~ pred_prob60, data=c60, 
    label=list(pred_prob60 = "Model 60+ only")) %>% 
  plot(smooth=TRUE) +
  ggplot2::labs(x="Threshold probability")

