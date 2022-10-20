#descriptive stats for study population - only those with blood tests 

#Load libraries
library(tidyverse)
library(Hmisc)

#clear R environment
rm(list=ls())

#load data 
load("N:\\DataAnalysis\\FormattedData\\Myeloma.RData")
load("N:\\DataAnalysis\\FormattedData\\MGUS.RData")
load("N:\\DataAnalysis\\FormattedData\\Controls.RData")


#load cleaned blood test results and only keep those with at least one results 
#load full data 
load("N:\\DataAnalysis\\FormattedData\\AllTestsCombined.RData")
#get ids separately for myelomas, MGUS and Controls 
myids<- tc %>% filter(case=="Myeloma") %>% select(id) %>% distinct()
mgusids<- tc %>% filter(case=="MGUS") %>% select(id) %>% distinct()
conids<- tc %>% filter(case=="Control") %>% select(id) %>% distinct()
table(tc$case)
tc %>% group_by(case) %>% summarise(n_distinct(id)) 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Myeloma cohort
my_samp <- myids %>% left_join(pat)

#recode ethnicity to have same groups for cases and controls
my_samp$ethnicity<-as.factor(my_samp$ethnicity)
my_samp <- my_samp %>% mutate(ethnicity = na_if(ethnicity, "Not Stated"))
my_samp$ethnicity<-as.factor(my_samp$ethnicity)
my_samp$ethnicity = droplevels(my_samp$ethnicity)
my_samp <- my_samp %>% 
  mutate(ethnicity = recode(ethnicity, 
                            "White - British or Irish" = "WB", 
                            "White - Any other White background" = "Other", 
                            "Mixed"= "Other", 
                            "Other Ethnic Groups" = "Other", 
                            "Black or Black British - African,Caribbean or Other" = "Black",
                            "Asian or Asian British - Indian or Pakistani" = "SA") ) 
table(my_samp$ethnicity, exclude=NULL)
etab2<-table(my_samp$ethnicity, exclude=NULL)
etp2<-my_samp$ethnicity %>% table(exclude=NULL) %>% prop.table() %>% `*`(100) %>% round(2)
cbind(etab2, etp2)
#tab sex
stab<-table(my_samp$sex, exclude=NULL)
stp<-my_samp$sex %>% table() %>% prop.table() %>% `*`(100) %>% round(2)
cbind(stab, stp)
#summarise age
summary(my_samp$age)
mean(my_samp$age)
sd(my_samp$age)
ggplot(data=my_samp, aes(age)) +
  geom_histogram(breaks=seq(35, 100, by=5))+
  labs(title="Age distribution", x="Age at diagnosis", y="Count") +
  ylim(c(0, 100))


#create age groups in 10 year age bands
summary(my_samp$age)
my_samp <- my_samp %>%
  mutate(age_cat=cut(age, breaks=c(30, 40, 50, 60, 70, 80, 90, 100), 
                     right=FALSE, 
                     label=c("30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90-99")))
table(my_samp$age_cat)

#tabulate performance status NAs
ptab<-table(my_samp$performance_lab, exclude=NULL)
ptp<-my_samp$performance_lab %>% table() %>% prop.table() %>% `*`(100) %>% round(2)
cbind(ptab, ptp)
#plot imd deciles
ggplot(data=my_samp, aes(imd2015)) +
  geom_bar() +
  labs(title="IMD 2015", y="Count")
ggplot(data=my_samp, aes(imd2019)) +
  geom_bar() +
  labs(title="IMD 2019", y="Count")
#tabulte IMD
table(my_samp$imd2015, exclude=NULL)
table(my_samp$imd2019, exclude=NULL)

#numbers by sex and age group 
my_samp %>% 
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                   right=FALSE, 
                   label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  select(sex, age_cat) %>% table()
my_samp %>% select(sex) %>% table()

#check units and reference ranges within each test
refrange <- pat_wbloods %>%
  group_by(test_name, units, ref_range) %>%
  summarise(n_test=n(), n_pat=n_distinct(id)) %>%
  arrange(desc(test_name), by.group=TRUE)

#summarise department
dept <- pat_wbloods %>%
  group_by(dept) %>%
  summarise(n_dept=n()) %>%
  arrange(desc(n_dept), by.group=TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#MGUS cohort
mgus_samp <- mgusids %>% left_join(mgus_pat)

dim(mgus_samp)
describe(mgus_samp)
#recode ethnic groups to be the same for cases and controls
mgus_samp$ethnicity<-as.factor(mgus_samp$ethnicity)
mgus_samp <- mgus_samp %>% 
  mutate(ethnicity = na_if(ethnicity, "Not Stated")) %>% 
  mutate(ethnicity = na_if(ethnicity, "Not Known")) 
mgus_samp$ethnicity<-as.factor(mgus_samp$ethnicity)
mgus_samp$ethnicity = droplevels(mgus_samp$ethnicity)
mgus_samp <- mgus_samp %>% 
  mutate(ethnicity = recode(ethnicity, 
                            "White - British, Irish or Other" = "WB", 
                            "Other or Mixed"= "Other", 
                            "Black or Black British - African,Caribbean or Other" = "Black",
                            "Asian or Asian British - Indian, Pakistani or Other" = "SA") ) 
table(mgus_samp$ethnicity, exclude=NULL)
etab2<-table(mgus_samp$ethnicity, exclude=NULL)
etp2<-mgus_samp$ethnicity %>% table(exclude=NULL) %>% prop.table() %>% `*`(100) %>% round(2)
cbind(etab2, etp2)
#tab sex
stab<-table(mgus_samp$sex, exclude=NULL)
stp<-mgus_samp$sex %>% table() %>% prop.table() %>% `*`(100) %>% round(2)
cbind(stab, stp)
#summarise age
summary(mgus_samp$age)
mean(mgus_samp$age)
sd(mgus_samp$age)
ggplot(data=mgus_samp, aes(age)) +
  geom_histogram(breaks=seq(35, 100, by=5))+
  labs(title="Age distribution", x="Age at diagnosis", y="Count") +
  ylim(c(0, 100))

#create age groups in 10 year age bands
summary(mgus_samp$age)
mgus_samp <- mgus_samp %>%
  mutate(age_cat=cut(age, breaks=c(30, 40, 50, 60, 70, 80, 90, 100), 
                     right=FALSE, 
                     label=c("30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90-99")))
table(mgus_samp$age_cat)

#plot imd deciles
ggplot(data=mgus_samp, aes(imd2015)) +
  geom_bar() +
  labs(title="IMD 2015", y="Count")
ggplot(data=mgus_samp, aes(imd2019)) +
  geom_bar() +
  labs(title="IMD 2019", y="Count")
#tabulte IMD
table(mgus_samp$imd2015, exclude=NULL)
table(mgus_samp$imd2019, exclude=NULL)

#numbers by sex and age group 
mgus_samp %>% 
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  select(sex, age_cat) %>% table()
mgus_samp %>% select(sex) %>% table()

#check units and reference ranges within each test
refrange_mgus <- mgus_pat_wbloods %>%
  group_by(test_name, units, ref_range) %>%
  summarise(n_test=n(), n_pat=n_distinct(id)) %>%
  arrange(desc(test_name), by.group=TRUE)

#summarise department
dept_mgus <- mgus_pat_wbloods %>%
  group_by(dept) %>%
  summarise(n_dept=n()) %>%
  arrange(desc(n_dept), by.group=TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Control cohort
con_samp <- conids %>% left_join(con)

#recode ethnicity  
con_samp$ethnicity<-as.factor(con_samp$ethnicity)
con_samp <- con_samp %>% 
  mutate(ethnicity = na_if(ethnicity, "Not Stated")) %>% 
  mutate(ethnicity = na_if(ethnicity, "Not Known")) 
con_samp$ethnicity<-as.factor(con_samp$ethnicity)
con_samp$ethnicity = droplevels(con_samp$ethnicity)
con_samp <- con_samp %>% 
  mutate(ethnicity = recode(ethnicity, 
                            "White - British or Irish" = "WB", 
                            "White - Any other White background" = "Other", 
                            "Mixed or Other Ethnic Groups" = "Other", 
                            "Black or Black British" = "Black",
                            "Asian or Asian British" = "SA") ) 
table(con_samp$ethnicity, exclude=NULL)
etab2<-table(con_samp$ethnicity, exclude=NULL)
etp2<-con_samp$ethnicity %>% table(exclude = NULL) %>% prop.table() %>% `*`(100) %>% round(2)
cbind(etab2, etp2)

#tab sex
stab<-table(con_samp$sex, exclude=NULL)
stp<-con_samp$sex %>% table() %>% prop.table() %>% `*`(100) %>% round(2)
cbind(stab, stp)
#summarise age
summary(con_samp$age)
mean(con_samp$age)
sd(con_samp$age)
ggplot(data=con_samp, aes(age)) +
  geom_histogram(breaks=seq(35, 100, by=5)) +
  labs(title="Age distribution", x="Age at diagnosis", y="Count") 
#create age groups in 10 year age bands
summary(con_samp$age)
con_samp <- con_samp %>%
  mutate(age_cat=cut(age, breaks=c(30, 40, 50, 60, 70, 80, 90, 100), 
                     right=FALSE, 
                     label=c("30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90-99")))
table(con_samp$age_cat)
#plot imd deciles
ggplot(data=con_samp, aes(imd2015)) +
  geom_bar() +
  labs(title="IMD 2015", y="Count")
ggplot(data=con_samp, aes(imd2019)) +
  geom_bar() +
  labs(title="IMD 2019", y="Count")
#tabulte IMD
table(con_samp$imd2015, exclude=NULL)
table(con_samp$imd2019, exclude=NULL)

#numbers by sex and age group 
con_samp %>% 
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  select(sex, age_cat) %>% table()
con_samp %>% select(sex) %>% table()

#check units and reference ranges within each test
refrange_con <- con_wbloods %>%
  group_by(test_name, units, ref_range) %>%
  summarise(n_test=n(), n_pat=n_distinct(id)) %>%
  arrange(desc(test_name), by.group=TRUE)

#summarise department
dept_con <- con_wbloods %>%
  group_by(dept) %>%
  summarise(n_dept=n()) %>%
  arrange(desc(n_dept), by.group=TRUE)

########################Co-morbidities

#myeloma
load("N:\\DataAnalysis\\FormattedData\\ComorbidMyeloma.Rdata")

#merge with patient data
my_comb<-left_join(my_samp, cm_my)
summarise(my_comb, n_distinct(id)) 
#replace missing with 0 as no recorded comorbidities
my_comb <- my_comb %>%
  mutate_at(vars(diab, liver, can, metscan, obese, dem, para, 
                 stroke, hyp, renal, mi, copd, chf, pvd, rh, totcm), 
            ~ replace_na(., 0))
#table of total no. comorbidities
cmtab<-table(my_comb$totcm)
t2<-prop.table(cmtab)
cbind(cmtab, t2)

#select comorbidity data only to create sums and %s
c1<- my_comb %>% 
  select("diab", "liver", "can", "metscan", "obese", "dem", "para", 
         "stroke", "hyp", "renal", "mi", "copd", "chf", "pvd", "rh")
#total with each condition
apply(c1, 2, sum)
#% with each condition
apply(c1, 2, mean)

#MGUS
load("N:\\DataAnalysis\\FormattedData\\ComorbidMGUS.Rdata")

#merge with patient data
mgus_comb<-left_join(mgus_samp, cm_mgus)
summarise(mgus_comb, n_distinct(id)) 
#replace missing with 0 as no recorded comorbidities
mgus_comb <- mgus_comb %>%
  mutate_at(vars(diab, liver, can, metscan, obese, dem, para, 
                 stroke, hyp, renal, mi, copd, chf, pvd, rh, totcm), 
            ~ replace_na(., 0))
#table of total no. comorbidities
cmtab<-table(mgus_comb$totcm)
t2<-prop.table(cmtab)
cbind(cmtab, t2)

#select comorbidity data only to create sums and %s
c1<- mgus_comb %>% 
  select("diab", "liver", "can", "metscan", "obese", "dem", "para", 
         "stroke", "hyp", "renal", "mi", "copd", "chf", "pvd", "rh")
#total with each condition
apply(c1, 2, sum)
#% with each condition
apply(c1, 2, mean)

#Controls
load("N:\\DataAnalysis\\FormattedData\\ComorbidControl.Rdata")

#merge with patient data
con_comb<-left_join(con_samp, cm_con)
summarise(con_comb, n_distinct(id)) 
#replace missing with 0 as no recorded comorbidities
con_comb <- con_comb %>%
  mutate_at(vars(diab, liver, obese, dem, para, 
                 stroke, hyp, renal, mi, copd, chf, pvd, rh, totcm), 
            ~ replace_na(., 0))
#table of total no. comorbidities
cmtab<-table(con_comb$totcm)
t2<-prop.table(cmtab)
cbind(cmtab, t2)
#select comorbidity data only to create sums and %s
c1<- con_comb %>% 
  select("diab", "liver", "obese", "dem", "para", 
         "stroke", "hyp", "renal", "mi", "copd", "chf", "pvd", "rh")
#total with each condition
apply(c1, 2, sum)
#% with each condition
apply(c1, 2, mean)

