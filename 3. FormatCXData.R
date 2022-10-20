###Extract cross sectional data for cases, controls and MGUS
#define index date as date of SEP/Paraprotein/Immunofixation
#use test closest to index data
#one test for each ID
#only keep 15 test to include in the model 
#reshape wide 

library(tidyverse)
library(Hmisc)

#clear R environment
rm(list=ls())

#save 
load("N:\\SharedFolder\\DataAnalysis\\FormattedData\\AllTestsCombined.RData")

tc %>% group_by(case) %>% summarise(n_distinct(id))
tc %>% group_by(case) %>% summarise(n())

#only keep 15 tests to include in models 
testkeep<- c("Albumin", "ALk-Phos", "ALT", "Basophil", "Calcium", "Creatinine", "CRP", "Eosinophil", 
             "Haemoglobin", "Lymphocyte", "Mean Cell Volume", "Monocyte", "Neutrophil", 
             "Platelets", "WCC")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Myeloma cases
cases<- tc %>% filter (case=="Myeloma")
table(cases$test)
table(cases$dept)

cases <-cases %>%  filter((test %in% testkeep))
table(cases$test)
table(cases$dept)

#filter tests based on index date based on SEP/Paraprotein/Immunofixation
load("N:\\SharedFolder\\DataAnalysis\\CrossSectional\\sepday.Rdata")

#Only keep those with valid paraprotein 
cases<- cases %>% inner_join(sepday)
n_distinct(cases$id)
#only keep tests on or before day of paraprotein test 
cases$days<- cases$days_prior*-1
cases<-cases %>% filter(days>sepday)
n_distinct(cases$id)

#counts of IDs with at least 1 test 
counts<- cases %>% group_by(test) %>% summarise(n_distinct(id))

#keep last measurement for each test (some have more than 1 test result on the same day)
keep <- cases %>% 
  group_by(test, id) %>% 
  filter(days_prior==max(days_prior)) %>% 
  mutate(count=row_number(), tot=n())
table(keep$count)
table(keep$tot)
#check those with more than 1 test on the same day
dup<-keep %>% filter(tot>=2)
table(dup$dept)

#only keep first test for duplicates 
keep <- keep %>% filter(count==1)

#need to check days prior remove those >2 years  from diagnosis 
hist(keep$days_prior)
summary(keep$days_prior)
#check those exlcuded as tests >2 years 
ty<-keep %>% filter(days_prior<=-730)
keep <- keep %>% filter(days_prior>-730)
ggplot(data=keep, aes(x=days_prior)) +
  geom_histogram() +
  facet_wrap(~test) 
n_distinct(keep$id)


#all data checks complete 
rm("counts", "dup")

#summarise dept for cases in cross sectional analysis
table(keep$dept)
deptMMCX <- keep %>%
  group_by(dept) %>%
  summarise(n_dept=n()) %>%
  arrange(desc(n_dept), by.group=TRUE)
write.csv(deptMMCX, "N:\\SharedFolder\\DataAnalysis\\CrossSectional\\MMCX_Dept.csv", row.names = TRUE)

#have 1 row per id so use wide format using pivot wider 
case_cx <- keep %>% select(id, case, age, sex, test, res, days_prior) %>% 
  pivot_wider(id_cols=c(id, sex, age, case), names_from=test, values_from=c(res, days_prior)) %>% ungroup() 
#rename columns - remove "res_" prefix 
colnames(case_cx) <-gsub("res_", "", colnames(case_cx))
#convert character vars to numeric
numnames<- c("Albumin", "ALk-Phos", "ALT", "Basophil", "Calcium", "Creatinine", "CRP", 
            "Eosinophil", "Haemoglobin", "Lymphocyte", "Mean Cell Volume", "Monocyte", "Neutrophil", 
            "Platelets", "WCC")
case_cx <-case_cx %>% mutate_at(vars(numnames), as.numeric) 

#rename variabls with spaces in 
case_cx <-case_cx %>%  rename(ALP="ALk-Phos", MCV="Mean Cell Volume")

#remove if any ids have all missing values 
btnames<-c("Albumin", "ALP", "ALT", "Basophil", "Calcium", "Creatinine", "CRP", "Eosinophil", 
           "Haemoglobin", "Lymphocyte", "MCV", "Monocyte", "Neutrophil", 
           "Platelets", "WCC") 
case_cx<- case_cx %>% select(id, age, sex, case, all_of(btnames))

#save as dataframes to reload later 
save(case_cx, file="N:\\SharedFolder\\DataAnalysis\\CrossSectional\\CasesCX.RData")

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##CONTROLS

con<- tc %>% filter (case=="Control")
table(con$test)
table(con$dept)

con <-con %>%  filter(test %in% testkeep)
table(con$test)
table(con$dept)

#counts of IDs with at least 1 test 
counts<- con %>% group_by(test) %>% summarise(n_distinct(id))

#keep last measurement for each test (some have more than 1 test result on the same day)
keep <- con %>% 
  group_by(test, id) %>% 
  filter(days_prior==max(days_prior)) %>% 
  mutate(count=row_number(), tot=n())
table(keep$count)
table(keep$tot)
#check those with more than 1 test on the same day
dup<-keep %>% filter(tot>=2)
table(dup$dept)

#only keep first test for duplicates 
keep <- keep %>% filter(count==1)

#need to check days prior remove those >2 year index date
hist(keep$days_prior)
k1<- keep %>% filter(days_prior<=-720)
table(k1$test)
ggplot(data=k1, aes(x=days_prior)) +
  geom_histogram() +
  facet_wrap(~test) 

#summarise dept for controls in cross sectional analysis
table(keep$dept)
deptconCX <- keep %>%
  group_by(dept) %>%
  summarise(n_dept=n()) %>%
  arrange(desc(n_dept), by.group=TRUE)
write.csv(deptconCX, "N:\\SharedFolder\\DataAnalysis\\CrossSectional\\CONCX_Dept.csv", row.names = TRUE)

########
#remove multiple rows by id
#pivot wider 
con_cx <- keep %>% select(id, case, age, sex, test, res, days_prior) %>% 
  pivot_wider(id_cols=c(id, sex, age, case), names_from=test, values_from=c(res, days_prior)) %>% ungroup() 
#rename columns - remove "res_" prefix 
colnames(con_cx) <-gsub("res_", "", colnames(con_cx))
#convert character vars to numeric
numnames<- c("Albumin", "ALk-Phos", "ALT", "Basophil", "Calcium", "Creatinine", "CRP", 
             "Eosinophil", "Haemoglobin", "Lymphocyte", "Mean Cell Volume", "Monocyte", "Neutrophil", 
             "Platelets", "WCC")
con_cx <-con_cx %>% mutate_at(vars(numnames), as.numeric) 
#rename variables with spaces in 
con_cx <-con_cx %>%  rename(ALP="ALk-Phos", MCV="Mean Cell Volume")
#any missing values across all columns
sapply(con_cx, anyNA)

#remove id any ids have all missing values 
colnames(con_cx)
btnames<-c("Albumin", "ALP", "ALT", "Basophil", "Calcium", "Creatinine", "CRP", "Eosinophil", 
           "Haemoglobin", "Lymphocyte", "MCV", "Monocyte", "Neutrophil", 
           "Platelets", "WCC") 

#exclude if all blood test results are missing
con_cx<- con_cx %>% 
  filter_at(vars(btnames), any_vars(! is.na(.)))

con_cx<- con_cx %>% select(id, age, sex, case, all_of(btnames))

save(con_cx, file="N:\\SharedFolder\\DataAnalysis\\CrossSectional\\ControlsCX.RData")

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##MGUS Cohort
mgus<- tc %>% filter (case=="MGUS")
mgus <-mgus %>%  filter((test %in% testkeep))
table(mgus$test)
#number of cases 
mgus %>% summarise(n_distinct(id))

#keep last measurement for each test (some have more than 1 test result on the same day)
keep <- mgus %>% 
  group_by(test, id) %>% 
  filter(days_prior==max(days_prior)) %>% 
  mutate(count=row_number(), tot=n())
table(keep$count)
table(keep$tot)
keep %>% group_by(case) %>% summarise(n_distinct(id))
#check those with more than 1 test on the same day
dup<-keep %>% filter(tot>=2)
table(dup$dept)

#only keep first test for duplicates 
keep <- keep %>% filter(count==1)

#need to check days prior remove those >2 years  from diagnosis 
hist(keep$days_prior)
summary(keep$days_prior)
#check those excluded as tests >2 years 
ty <- keep %>% filter(days_prior<=-730) 
keep <- keep %>% filter(days_prior>-730)
ggplot(data=keep, aes(x=days_prior)) +
  geom_histogram() +
  facet_wrap(~test) 
n_distinct(keep$id)
keep %>% group_by(case) %>% summarise(n_distinct(id))

#all data checks complete 
rm("dup", "ty")

#summarise dept for MGUS in cross sectional analysis
table(keep$dept)
deptMGUSCX <- keep %>%
  group_by(dept) %>%
  summarise(n_dept=n()) %>%
  arrange(desc(n_dept), by.group=TRUE)
write.csv(deptMGUSCX, "N:\\SharedFolder\\DataAnalysis\\CrossSectional\\MGUSCX_Dept.csv", row.names = TRUE)

#remove multiple rows by id
#pivot wider 
mgus_cx <- keep %>% select(id, case, age, sex, test, res, days_prior) %>% 
  pivot_wider(id_cols=c(id, sex, age, case), names_from=test, values_from=c(res, days_prior)) %>% ungroup() 
#rename columns - remove "res_" prefix 
colnames(mgus_cx) <-gsub("res_", "", colnames(mgus_cx))
#convert character vars to numeric
numnames<- c("Albumin", "ALk-Phos", "ALT", "Basophil", "Calcium", "Creatinine", "CRP", 
             "Eosinophil", "Haemoglobin", "Lymphocyte", "Mean Cell Volume", "Monocyte", "Neutrophil", 
             "Platelets", "WCC")
mgus_cx <-mgus_cx %>% mutate_at(vars(all_of(numnames)), as.numeric) 

#rename variabls with spaces in 
mgus_cx <-mgus_cx %>%  rename(ALP="ALk-Phos", MCV="Mean Cell Volume")
#any missing values across all columns
sapply(mgus_cx, anyNA)

btnames<-c("Albumin", "ALP", "ALT", "Basophil", "Calcium", "Creatinine", "CRP", "Eosinophil", 
           "Haemoglobin", "Lymphocyte", "MCV", "Monocyte", "Neutrophil", 
           "Platelets", "WCC") 
#exclude if all blood test results are missing
mgus_cx<- mgus_cx %>% 
  filter_at(vars(all_of(btnames)), any_vars(! is.na(.)))

mgus_cx<- mgus_cx %>% select(id, age, sex, case, all_of(btnames))

#save as dataframes to reload later 
save(mgus_cx, file="N:\\SharedFolder\\DataAnalysis\\CrossSectional\\MGUSCX.RData")

