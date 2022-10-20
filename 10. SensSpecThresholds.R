#code to run on model sens and spec to get PPV, NPV and number needed to test to detect 1 case for different model thresholds
#across a range of thresholds 0.1 to 0.5 and different prevalence estimates

library(tidyverse)
library(openxlsx)

rm(list=ls())

#data sens and spec for all models
s<-read.xlsx("ModelComparison_allthresholds.xlsx")

#use 3 different prevalence estimates
#60 per 100,000
pop<-100000
cases<-60
prev<-cases/pop

s60<- s %>% mutate(
    PPV=(sens*prev)/(sens*prev +(1-spec)*(1-prev)),
    NPV=(spec*(1-prev))/(((1-sens)*prev) + (spec*(1-prev))),
    TP=cases*sens,
    FN=cases-TP,
    TN=(pop-cases)*spec,
    FP=pop-cases-TN,
    reflex=TP+FP,
    NND=reflex/TP)
write.xlsx(s60, file="DiagStatsPrev60.xlsx", colNames = TRUE, rowNames = TRUE)

#15 per 100,000
pop<-100000
cases<-15
prev<-cases/pop

s15<- s %>% mutate(
  PPV=(sens*prev)/(sens*prev +(1-spec)*(1-prev)),
  NPV=(spec*(1-prev))/(((1-sens)*prev) + (spec*(1-prev))),
  TP=cases*sens,
  FN=cases-TP,
  TN=(pop-cases)*spec,
  FP=pop-cases-TN,
  reflex=TP+FP,
  NND=reflex/TP)
write.xlsx(s15, file="DiagStatsPrev15.xlsx", colNames = TRUE, rowNames = TRUE)

#10 per 100,000
pop<-100000
cases<-10
prev<-cases/pop

s10<- s %>% mutate(
  PPV=(sens*prev)/(sens*prev +(1-spec)*(1-prev)),
  NPV=(spec*(1-prev))/(((1-sens)*prev) + (spec*(1-prev))),
  TP=cases*sens,
  FN=cases-TP,
  TN=(pop-cases)*spec,
  FP=pop-cases-TN,
  reflex=TP+FP,
  NND=reflex/TP)
write.xlsx(s10, file="DiagStatsPrev10.xlsx", colNames = TRUE, rowNames = TRUE)


#combine all and plot
s60$prev<-60
s15$prev<-15
s10$prev<-10

allP<-rbind(s60, s15, s10)
table(allP$Model)

#order models 
allP$Model<-as.factor(allP$Model)
table(allP$Model)
allP$Model<-factor(allP$Model, levels=c("Full ", "noCa", "over60", "MGUS", "MGUSnoCa", "MGUSover60", "MGUSover60noCa"))
table(allP$Model)

#plots across thresholds comparing different models and metrics

#including all prevalence
#number needed to diagnose
#add labels to prev
prev.labs<-c ("Prevalence = 10", "Prevalence = 15", "Prevalence = 60")
names(prev.labs) <- c("10", "15", "60")

# Plots to include in paper 
#number needed to diagnose 1 case
ggplot(d=allP, aes(x=Threshold, y=NND, color=Model)) +
  geom_point() +
  facet_grid(.~prev, labeller=as_labeller(prev.labs)) +
  ggtitle("Number needed to test to diagnose 1 case")
ggsave(dpi=300, file="NNDpaperplot.tiff")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#these plots have not been included 
#number to reflex test 
ggplot(d=allP, aes(x=Threshold, y=reflex, color=Model)) +
  geom_point() +
  facet_grid(.~prev, labeller=as_labeller(prev.labs)) +
  ggtitle("Total to reflex test")

#number of cancers detected 
ggplot(d=allP, aes(x=Threshold, y=TP, color=Model)) +
  geom_point() +
  facet_grid(.~prev, labeller=as_labeller(prev.labs)) +
  ggtitle("Number of cancers detected")

#number of cancers missed
ggplot(d=allP, aes(x=Threshold, y=FN, color=Model)) +
  geom_point() +
  facet_grid(.~prev, labeller=as_labeller(prev.labs)) +
  ggtitle("Number of cancers missed")

