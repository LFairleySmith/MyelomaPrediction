#Total number of tests and median per patient
#Two time periods included: 1. Up to 5 years prior to index date, 
#                           2. Up to 1-year prior to index date 
#Number of tests in 90 day intervals to create heatmaps

library(tidyverse)

#clear R environment
rm(list=ls())

#load full data 
load("N:\\DataAnalysis\\FormattedData\\AllTestsCombined.RData")
#number with at least one blood test
tc %>% group_by(case) %>% summarise(n_distinct(id))

#create table with total number of tests and patients per test in past 5-years
bt_5yr <- tc %>%
  group_by(id, test) %>% 
  mutate(count=row_number(), totbt=n()) %>%
  filter(count==1) %>%   
  group_by(test, case) %>% 
  summarise(n_test=sum(totbt), n_pat=n_distinct(id), medbt=median(totbt),
            p25bt=quantile(totbt, probs=c(0.25)), p75bt=quantile(totbt, probs=c(0.75))) %>%
  arrange(test, case, by.group=TRUE)
write.csv(bt_5yr, "N:\\DataAnalysis\\FormattedData\\test_5yrs.csv", row.names = TRUE)

#create table with total number of tests and patients per test in past year
bt_1yr <- tc %>%
  filter(days_prior>=-365) %>% 
  group_by(id, test) %>% 
  mutate(count=row_number(), totbt=n()) %>%
  filter(count==1) %>%   
  group_by(test, case) %>% 
  summarise(n_test=sum(totbt), n_pat=n_distinct(id), medbt=median(totbt),
            p25bt=quantile(totbt, probs=c(0.25)), p75bt=quantile(totbt, probs=c(0.75))) %>%
  arrange(test, case, by.group=TRUE)
write.csv(bt_1yr, "N:\\DataAnalysis\\FormattedData\\test_1yr.csv", row.names = TRUE)

######Number of tests in 90 day intervals
#exclude myeloma specific tests (SFLC, B2)
table(tc$test)
tc<- tc %>% 
  filter(!test %in% c("Beta2microglobulin", "Free Kappa LC", "Free Lambda LC", "Serum Free LC"))  
table(tc$test)

summary(tc$days_prior)
molab<-c("-1800", "-1710", "-1620", "-1530", "-1440", 
          "-1350", "-1260", "-1170", "-1080", "-990", 
          "-900", "-810", "-720", "-630", "-540", "-450",
          "-360", "-270", "-180", "-90")

tc<-tc %>%
  mutate(days90=cut(days_prior, breaks=seq(-1800, 0, by=90),
                         right=TRUE, label=molab))
table(tc$days90, tc$case)

#generate 90days totals for each patient 
b_tot<- tc %>%
  group_by(id, test, days90) %>% 
  mutate(count=row_number(), totbt=n()) %>%
  filter(count==1) %>% 
  group_by(test, days90, case) %>% 
  summarise(n_test=sum(totbt), n_pat=n_distinct(id)) %>% 
  arrange(test, case, days90, by.group=TRUE)

summary(subset(b_tot, case=="Myeloma")$n_test)
summary(subset(b_tot, case=="MGUS")$n_test)
summary(subset(b_tot, case=="Control")$n_test)

summary(subset(b_tot, case=="Myeloma")$n_pat)
summary(subset(b_tot, case=="MGUS")$n_pat)
summary(subset(b_tot, case=="Control")$n_pat)

#heatmap plot - needs 4 one for cases and controls, no. tests and no pateints
library(viridis)

#myeloma number of tests
p1<- ggplot(subset(b_tot, case=="Myeloma"), aes(x=days90, ,y=test, fill=n_test)) +
  geom_tile() +
  scale_fill_viridis(discrete=FALSE, limits=c(0, 2500), 
                     breaks=seq(0, 2500, by=1000))+
  coord_fixed()+
  guides(fill=guide_colorbar(title="Total number of tests")) +
  labs(title="Myeloma - Total number of tests", x="Days prior to diagnosis", 
       y="Blood Test") +
  theme(axis.text.x = element_text(angle=45, size=7.5))
p1

#Controls -  number of tests
p2<-ggplot(subset(b_tot, case=="Control"), aes(x=days90, ,y=test, fill=n_test)) +
  geom_tile() +
  scale_fill_viridis(discrete=FALSE, limits=c(0, 12500), 
                     breaks=seq(0, 12500, by=2500))+
  coord_fixed()+
  guides(fill=guide_colorbar(title="Total number of tests")) +
  labs(title="Controls - Total number of tests", x="Days prior to diagnosis", 
       y="Blood Test") +
  theme(axis.text.x = element_text(angle=45, size=7.5))
p2

#myeloma number of patients
p3<- ggplot(subset(b_tot, case=="Myeloma"), aes(x=days90, ,y=test, fill=n_pat)) +
  geom_tile() +
  scale_fill_viridis(discrete=FALSE, limits=c(0, 400), 
                     breaks=seq(0, 400, by=50))+
  coord_fixed()+
  guides(fill=guide_colorbar(title="Total number of tests")) +
  labs(title="Myeloma - number of patients with at least 1 test", x="Days prior to diagnosis", 
       y="Blood Test") +
  theme(axis.text.x = element_text(angle=45, size=7.5))
p3

#Controls  number of patients
p4<- ggplot(subset(b_tot, case=="Control"), aes(x=days90, ,y=test, fill=n_pat)) +
  geom_tile() +
  scale_fill_viridis(discrete=FALSE, limits=c(0, 2500), 
                     breaks=seq(0, 2000, by=500))+
  coord_fixed()+
  guides(fill=guide_colorbar(title="Total number of tests")) +
  labs(title="Controls - number of patients with at least 1 test", x="Days prior to diagnosis", 
       y="Blood Test") +
  theme(axis.text.x = element_text(angle=45, size=7.5))
p4

pdf("N:\\DataAnalysis\\FormattedData\\NoTestsHeatMaps.pdf")
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()

#90 day interval plots

#exclude categorical variables and calculate means per 90 day period
m90<- tc %>% 
  filter(!test %in% c("Derived Fibrinogen", "eGFR")) %>% 
  filter(!case=="MGUS") %>% 
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
table(m90$test)
table(m90$case)

plot_list=list()
for (x in unique(m90$test)) {
  p= ggplot(m90[m90$test==x,],
               aes(x=days90, y=mean90, color=case)) +
          geom_point() +
          geom_line(aes(group=case)) +
          theme(axis.text.x = element_text(angle=45, size=7.5)) +
          labs(title=x, 
               x="Days prior to diagnosis", 
               y=x)
  plot_list[[x]]=p
}

#save plots as tiff files 
for (x in unique(m90$test)) {
  file_name=paste("plot_", x, ".tiff", sep="")
  tiff(file_name)
  print(plot_list[[x]])
  dev.off()
}  
  
