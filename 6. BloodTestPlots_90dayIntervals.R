#Trends plots over time means and 95C%Is in 90 day intervals 

library(tidyverse)
library(ggpubr)


#clear R environment
rm(list=ls())

#load full data 
load("N:\\SharedFolder\\DataAnalysis\\FormattedData\\AllTestsCombined.RData")
#create labels and groups for 90day intervals 
molab<-c("-1800", "-1710", "-1620", "-1530", "-1440", 
          "-1350", "-1260", "-1170", "-1080", "-990", 
          "-900", "-810", "-720", "-630", "-540", "-450",
          "-360", "-270", "-180", "-90")

tc<-tc %>%
  mutate(days90=cut(days_prior, breaks=seq(-1800, 0, by=90),
                         right=TRUE, label=molab))
table(tc$days90, tc$case)

#only keep for 15 blood tests used in main modeling 
table(tc$test)
tc<- tc %>% mutate(test= recode(test, "ALk-Phos"="ALP", "Mean Cell Volume"="MCV"))

btnames<-c("Albumin", "ALP", "ALT", "Basophil", "Calcium", "Creatinine", "CRP", "Eosinophil", 
           "Haemoglobin", "Lymphocyte", "MCV", "Monocyte", "Neutrophil", 
           "Platelets", "WCC") 


#Generate summary data to extract including Myeloma, MGUS and Controls
# means (SE) per 90 day period for each test 
m90<- tc %>% 
  filter(test %in% btnames) %>% 
  mutate(res=as.numeric(res))  %>% 
  group_by(test, case, days90) %>% 
  summarise (mean90=mean(res, na.rm=T), 
             sd=sd(res, na.rm=T),
             count=n())  %>%
  mutate(se=sd/sqrt(count),
         lci=mean90-1.96*se, 
         uci=mean90+1.96*se) %>% 
  as.data.frame()
m90$case<-factor(m90$case)

#exclude MGUS cohort to plot
m90<- m90 %>% filter(!case=="MGUS") 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#plot with error bars over time for each blood test 
#Myeloma and controls only
#Need to set scale for each plot

#Albumin
albm<- ggplot(data=subset(m90, test=="Albumin"), aes(x=days90, y=mean90, color=case)) +
  geom_point() +
  geom_line(aes(group=case)) +
  geom_errorbar(aes(ymin=lci, ymax=uci), width=0.2) +
  labs(title="Albumin", x="Time prior to diagnosis") +
  scale_y_continuous(name = "Albumin", limits=c(30, 45), breaks=seq(30, 45, 5)) +
  theme(axis.text.x = element_text(angle=90)) +
  theme(legend.title = element_blank(),
        legend.position = c(0.85, 0.8), 
      legend.background=element_rect(fill="white", color="black"),
        aspect.ratio =0.5) +
  scale_color_brewer(palette="Set1")
albm
ggsave("N:\\SharedFolder\\DataAnalysis\\Trajectories\\Mean90dayPlots\\Formatted\\Albumin.tiff", plot=albm, dpi=300)

#ALP
alpm<- ggplot(data=subset(m90, test=="ALP"), aes(x=days90, y=mean90, color=case)) +
  geom_point() +
  geom_line(aes(group=case)) +
  geom_errorbar(aes(ymin=lci, ymax=uci), width=0.2) +
  labs(title="ALP", x="Time prior to diagnosis") +
  scale_y_continuous(name = "ALP", limits=c(150, 350)) +
  theme(axis.text.x = element_text(angle=90)) +
  theme(legend.title = element_blank(),
        legend.position = c(0.85, 0.8), 
        legend.background=element_rect(fill="white", color="black"),
                aspect.ratio =0.5) +
  scale_color_brewer(palette="Set1")
alpm
ggsave("N:\\SharedFolder\\DataAnalysis\\Trajectories\\Mean90dayPlots\\Formatted\\ALP.tiff", plot=alpm, dpi=300)

#ALT
altm<- ggplot(data=subset(m90, test=="ALT"), aes(x=days90, y=mean90, color=case)) +
  geom_point() +
  geom_line(aes(group=case)) +
  geom_errorbar(aes(ymin=lci, ymax=uci), width=0.2) +
  labs(title="ALT", x="Time prior to diagnosis") +
  scale_y_continuous(name = "ALT", limits=c(0, 100)) +
  theme(axis.text.x = element_text(angle=90)) +
  theme(legend.title = element_blank(),
        legend.position = c(0.85, 0.8), 
        legend.background=element_rect(fill="white", color="black"),
        aspect.ratio =0.5) +
  scale_color_brewer(palette="Set1")
altm
ggsave("N:\\SharedFolder\\DataAnalysis\\Trajectories\\Mean90dayPlots\\Formatted\\ALT.tiff", plot=altm, dpi=300)

#Basophils
basm<- ggplot(data=subset(m90, test=="Basophil"), aes(x=days90, y=mean90, color=case)) +
  geom_point() +
  geom_line(aes(group=case)) +
  geom_errorbar(aes(ymin=lci, ymax=uci), width=0.2) +
  labs(title="Basophil", x="Time prior to diagnosis") +
  scale_y_continuous(name = "Basophil", limits=c(0, 0.1)) +
  theme(axis.text.x = element_text(angle=90)) +
  theme(legend.title = element_blank(),
        legend.position = c(0.85, 0.8), 
        legend.background=element_rect(fill="white", color="black"),
        aspect.ratio =0.5) +
  scale_color_brewer(palette="Set1")
basm
ggsave("N:\\SharedFolder\\DataAnalysis\\Trajectories\\Mean90dayPlots\\Formatted\\Basophil.tiff", plot=basm, dpi=300)

#Calcium
calm<- ggplot(data=subset(m90, test=="Calcium"), aes(x=days90, y=mean90, color=case)) +
  geom_point() +
  geom_line(aes(group=case)) +
  geom_errorbar(aes(ymin=lci, ymax=uci), width=0.2) +
  labs(title="Calcium", x="Time prior to diagnosis") +
  scale_y_continuous(name = "Calcium", limits=c(2, 3), breaks=seq(2, 3, 0.5)) +
  theme(axis.text.x = element_text(angle=90)) +
  theme(legend.title = element_blank(),
        legend.position = c(0.85, 0.8), 
        legend.background=element_rect(fill="white", color="black"),
        aspect.ratio =0.5) +
  scale_color_brewer(palette="Set1")
calm
ggsave("N:\\SharedFolder\\DataAnalysis\\Trajectories\\Mean90dayPlots\\Formatted\\Calcium.tiff", plot=calm, dpi=300)

#Creatinine
creatm<- ggplot(data=subset(m90, test=="Creatinine"), aes(x=days90, y=mean90, color=case)) +
  geom_point() +
  geom_line(aes(group=case)) +
  geom_errorbar(aes(ymin=lci, ymax=uci), width=0.2) +
  labs(title="Creatinine", x="Time prior to diagnosis") +
  scale_y_continuous(name = "Creatinine", limits=c(50, 300), breaks=seq(50, 300, 50)) +
  theme(axis.text.x = element_text(angle=90)) +
  theme(legend.title = element_blank(),
        legend.position = c(0.85, 0.8), 
        legend.background=element_rect(fill="white", color="black"),
        aspect.ratio =0.5) +
  scale_color_brewer(palette="Set1")
creatm
ggsave("N:\\SharedFolder\\DataAnalysis\\Trajectories\\Mean90dayPlots\\Formatted\\Creatinine.tiff", plot=creatm, dpi=300)

#CRP
crpm<- ggplot(data=subset(m90, test=="CRP"), aes(x=days90, y=mean90, color=case)) +
  geom_point() +
  geom_line(aes(group=case)) +
  geom_errorbar(aes(ymin=lci, ymax=uci), width=0.2) +
  labs(title="CRP", x="Time prior to diagnosis") +
  scale_y_continuous(name = "CRP", limits=c(0, 120)) +
  theme(axis.text.x = element_text(angle=90)) +
  theme(legend.title = element_blank(),
        legend.position = c(0.85, 0.8), 
        legend.background=element_rect(fill="white", color="black"),
              aspect.ratio =0.5) +
  scale_color_brewer(palette="Set1")
crpm
ggsave("N:\\SharedFolder\\DataAnalysis\\Trajectories\\Mean90dayPlots\\Formatted\\CRP.tiff", plot=crpm, dpi=300)

#Eosinophiol
eosm<- ggplot(data=subset(m90, test=="Eosinophil"), aes(x=days90, y=mean90, color=case)) +
  geom_point() +
  geom_line(aes(group=case)) +
  geom_errorbar(aes(ymin=lci, ymax=uci), width=0.2) +
  labs(title="Eosinophil", x="Time prior to diagnosis") +
  scale_y_continuous(name = "Eosinophil", limits=c(0, 0.3)) +
  theme(axis.text.x = element_text(angle=90)) +
  theme(legend.title = element_blank(),
        legend.position = c(0.85, 0.8), 
        legend.background=element_rect(fill="white", color="black"),
        aspect.ratio =0.5) +
  scale_color_brewer(palette="Set1")
eosm
ggsave("N:\\SharedFolder\\DataAnalysis\\Trajectories\\Mean90dayPlots\\Formatted\\Eosinophil.tiff", plot=eosm, dpi=300)

#Haemoglobin
hbm<- ggplot(data=subset(m90, test=="Haemoglobin"), aes(x=days90, y=mean90, color=case)) +
  geom_point() +
  geom_line(aes(group=case)) +
  geom_errorbar(aes(ymin=lci, ymax=uci), width=0.2) +
  labs(title="Haemoglobin", x="Time prior to diagnosis") +
  scale_y_continuous(name = "Haemoglobin", limits=c(100, 150), breaks=seq(100, 150, 10)) +
  theme(axis.text.x = element_text(angle=90)) +
  theme(legend.title = element_blank(),
        legend.position = c(0.85, 0.8), 
        legend.background=element_rect(fill="white", color="black"),
        aspect.ratio =0.5) +
  scale_color_brewer(palette="Set1")
hbm
ggsave("N:\\SharedFolder\\DataAnalysis\\Trajectories\\Mean90dayPlots\\Formatted\\Haemoglobin.tiff", plot=hbm, dpi=300)

#Lymphocyte
lymphm<- ggplot(data=subset(m90, test=="Lymphocyte"), aes(x=days90, y=mean90, color=case)) +
  geom_point() +
  geom_line(aes(group=case)) +
  geom_errorbar(aes(ymin=lci, ymax=uci), width=0.2) +
  labs(title="Lymphocytes", x="Time prior to diagnosis") +
  scale_y_continuous(name = "Lymphocyte", limits=c(1, 2.5), breaks=seq(1, 2.5, 0.5)) +
  theme(axis.text.x = element_text(angle=90)) +
  theme(legend.title = element_blank(),
        legend.position = c(0.85, 0.8), 
        legend.background=element_rect(fill="white", color="black"),
        aspect.ratio =0.5) +
  scale_color_brewer(palette="Set1")
lymphm
ggsave("N:\\SharedFolder\\DataAnalysis\\Trajectories\\Mean90dayPlots\\Formatted\\Lymphocytes.tiff", plot=lymphm, dpi=300)

#MCV
mcvm<- ggplot(data=subset(m90, test=="MCV"), aes(x=days90, y=mean90, color=case)) +
  geom_point() +
  geom_line(aes(group=case)) +
  geom_errorbar(aes(ymin=lci, ymax=uci), width=0.2) +
  labs(title="MCV", x="Time prior to diagnosis") +
  scale_y_continuous(name = "MCV", limits=c(80, 100), breaks=seq(80, 100, 10)) +
  theme(axis.text.x = element_text(angle=90)) +
  theme(legend.title = element_blank(),
        legend.position = c(0.85, 0.4), 
        legend.background=element_rect(fill="white", color="black"),
        aspect.ratio =0.5) +
  scale_color_brewer(palette="Set1")
mcvm
ggsave("N:\\SharedFolder\\DataAnalysis\\Trajectories\\Mean90dayPlots\\Formatted\\MCV.tiff", plot=mcvm, dpi=300)

#Monocyte
monom<- ggplot(data=subset(m90, test=="Monocyte"), aes(x=days90, y=mean90, color=case)) +
  geom_point() +
  geom_line(aes(group=case)) +
  geom_errorbar(aes(ymin=lci, ymax=uci), width=0.2) +
  labs(title="Monocytes", x="Time prior to diagnosis") +
  scale_y_continuous(name = "Monocyte", limits=c(0, 1), breaks=seq(0, 1, 0.25)) +
  theme(axis.text.x = element_text(angle=90)) +
  theme(legend.title = element_blank(),
        legend.position = c(0.85, 0.8), 
        legend.background=element_rect(fill="white", color="black"),
        aspect.ratio =0.5) +
  scale_color_brewer(palette="Set1")
monom
ggsave("N:\\SharedFolder\\DataAnalysis\\Trajectories\\Mean90dayPlots\\Formatted\\Monocytes.tiff", plot=monom, dpi=300)

#Neutrophil
neutm<- ggplot(data=subset(m90, test=="Neutrophil"), aes(x=days90, y=mean90, color=case)) +
  geom_point() +
  geom_line(aes(group=case)) +
  geom_errorbar(aes(ymin=lci, ymax=uci), width=0.2) +
  labs(title="Neutrophil", x="Time prior to diagnosis") +
  scale_y_continuous(name = "Neutrophil", limits=c(0, 10)) +
  theme(axis.text.x = element_text(angle=90)) +
  theme(legend.title = element_blank(),
        legend.position = c(0.85, 0.8), 
        legend.background=element_rect(fill="white", color="black"),
        aspect.ratio =0.5) +
  scale_color_brewer(palette="Set1")
neutm
ggsave("N:\\SharedFolder\\DataAnalysis\\Trajectories\\Mean90dayPlots\\Formatted\\Neutrophil.tiff", plot=neutm, dpi=300)

#Platelets
pltm<- ggplot(data=subset(m90, test=="Platelets"), aes(x=days90, y=mean90, color=case)) +
  geom_point() +
  geom_line(aes(group=case)) +
  geom_errorbar(aes(ymin=lci, ymax=uci), width=0.2) +
  labs(title="Platelets", x="Time prior to diagnosis") +
  scale_y_continuous(name = "Platelets", limits=c(200, 450), breaks=seq(200, 450, 50)) +
  theme(axis.text.x = element_text(angle=90)) +
  theme(legend.title = element_blank(),
        legend.position = c(0.85, 0.8), 
        legend.background=element_rect(fill="white", color="black"),
        aspect.ratio =0.5) +
  scale_color_brewer(palette="Set1")
pltm
ggsave("N:\\SharedFolder\\DataAnalysis\\Trajectories\\Mean90dayPlots\\Formatted\\Platelets.tiff", plot=pltm, dpi=300)

#WCC
wccm<- ggplot(data=subset(m90, test=="WCC"), aes(x=days90, y=mean90, color=case)) +
  geom_point() +
  geom_line(aes(group=case)) +
  geom_errorbar(aes(ymin=lci, ymax=uci), width=0.2) +
  labs(title="WCC", x="Time prior to diagnosis") +
  scale_y_continuous(name = "WCC", limits=c(5, 12), breaks=seq(5, 12, 1)) +
  theme(axis.text.x = element_text(angle=90)) +
  theme(legend.title = element_blank(),
        legend.position = c(0.85, 0.8), 
        legend.background=element_rect(fill="white", color="black"),
        aspect.ratio =0.5) +
  scale_color_brewer(palette="Set1")
wccm
ggsave("N:\\SharedFolder\\DataAnalysis\\Trajectories\\Mean90dayPlots\\Formatted\\WCC.tiff", plot=wccm, dpi=300)
