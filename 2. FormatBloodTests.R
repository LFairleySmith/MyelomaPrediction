#Code to collate all test data: myeloma, MGUS and contols
#data in long format
#include age and sex for descriptive stats and plots
#each test saved separately 

library(tidyverse)
library(Hmisc)

#clear R environment
rm(list=ls())

#load saved myeloma, mgus and control data 
load("N:\\DataAnalysis\\FormattedData\\Myeloma.RData")
load("N:\\DataAnalysis\\FormattedData\\MGUS.RData")
load("N:\\DataAnalysis\\FormattedData\\Controls.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Albumin

#MYELOMA CASES
alb_my <-filter(pat_wbloods, test_name=="Albumin")
alb_my$res<-alb_my$result
alb_my$res[alb_my$res == ""] <- NA
alb_my$res <- as.numeric(alb_my$res)
alb_my<- alb_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (alb=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
alb_my <- within(alb_my, {
  id<-factor(id)
})
alb_my$sex<-factor(alb_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
alb_my$case<-1

#MGUS CASES
alb_mgus <-filter(mgus_pat_wbloods, test_name=="Albumin")
alb_mgus$res<-alb_mgus$result
alb_mgus$res[alb_mgus$res == ""] <- NA
alb_mgus$res <- as.numeric(alb_mgus$res)
alb_mgus <- alb_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (alb=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
alb_mgus <- within(alb_mgus, {
  id<-factor(id)
})
alb_mgus$sex<-factor(alb_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
alb_mgus$case<-2

#Controls
alb_con <-filter(con_wbloods, test_name=="Albumin")
alb_con %>% summarise(n_distinct(id))
alb_con$res<-alb_con$result
alb_con$res[alb_con$res == ""] <- NA
alb_con$res <- as.numeric(alb_con$res)
alb_con <- alb_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (alb=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
alb_con <- within(alb_con, {
  id<-factor(id)
})
alb_con$sex<-factor(alb_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
alb_con$case<-3

#append myeloma, MGUS and controls 
alb<-rbind(alb_my, alb_mgus, alb_con)
#case as factor variable
alb$case<-factor(alb$case, levels=c('1', '2', '3'), 
                 labels=c('Myeloma', 'MGUS', 'Control'))

table(alb$case)

#save as R file to use in further analysis
save(alb, file="N:\\DataAnalysis\\Trajectories\\Albumin\\Alb_long.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Alkaline Phosphatase (ALP)

#MYELOMA CASES
alk_my <-filter(pat_wbloods, test_name=="Alk phos"|test_name=="Alk Phos")
alk_my$res<-alk_my$result
alk_my$res[alk_my$res == ""] <- NA
alk_my$res <- as.numeric(alk_my$res)
alk_my <- alk_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (alk=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
alk_my <- within(alk_my, {
  id<-factor(id)
})
alk_my$sex<-factor(alk_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
alk_my$case<-1

#MGUS CASES
alk_mgus <-filter(mgus_pat_wbloods, test_name=="Alk phos"|test_name=="Alk Phos")
alk_mgus$res<-alk_mgus$result
alk_mgus$res[alk_mgus$res == ""] <- NA
alk_mgus$res <- as.numeric(alk_mgus$res)
alk_mgus <- alk_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (alk=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
alk_mgus <- within(alk_mgus, {
  id<-factor(id)
})
alk_mgus$sex<-factor(alk_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
alk_mgus$case<-2

#CONTROLS
alk_con <-filter(con_wbloods, test_name=="Alk phos"|test_name=="Alk Phos")
alk_con %>% summarise(n_distinct(id))
alk_con$res<-alk_con$result
alk_con$res[alk_con$res == ""] <- NA
alk_con$res <- as.numeric(alk_con$res)
alk_con <- alk_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (alk=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
alk_con <- within(alk_con, {
  id<-factor(id)
})
alk_con$sex<-factor(alk_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
alk_con$case<-3 

#append myeloma, MGUS and controls 
alkp<-rbind(alk_my, alk_mgus, alk_con)
#case as factor variable
alkp$case<-factor(alkp$case, levels=c('1', '2', '3'), 
                  labels=c('Myeloma', 'MGUS', 'Control'))
table(alkp$case)

#remove values >1000
alkp %>% subset(alk>1000) %>% count()
alkp %>% subset(alk>1000) %>% group_by(case) %>% count()
alkp %>% subset(alk>1000) %>% group_by(case) %>% distinct(id) %>% count()

alkp<- alkp %>% subset(alk<=1000)

#save as R file 
save(alkp, file="N:\\DataAnalysis\\Trajectories\\ALK-Phos\\ALKP_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Alanine transaminase (ALT)

#MYELOMA CASES
alt_my <-filter(pat_wbloods, test_name=="ALT")
table(alt_my$ref_range)
alt_my$flag<-ifelse(grepl("<", alt_my$result, ignore.case = T), 1, 0)
alt_my %>% tally(flag) #N=30
table(subset(alt_my, flag==1)$result)
alt_my$res<-alt_my$result
alt_my$res[alt_my$res == ""] <- NA
alt_my<-alt_my %>% drop_na(res)
alt_my$res <- sub("<", NA, alt_my$res)
alt_my$res <- as.numeric(alt_my$res)
# lowest rescorded ALT value is 3, for those with values recorded as <9, impute random value between 3 and 9
alt_my <- alt_my %>% 
  mutate(m=ifelse(is.na(res), 1, 0))
alt_my %>% tally(m)
#create columns of 1s to sum over
alt_my$c<-1
#ipute imissing values as random value between 0 and 5
alt_my <- alt_my %>% 
  mutate(alt_im =
           case_when(
             m==0 ~ res,
             m==1 ~ runif(sum(c), min=3, max=9)))
alt_my <- alt_my %>%
  select(id, alt_im, days_prior, age, sex, flag, dept) %>%
  rename (alt=alt_im,) %>% 
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
alt_my <- within(alt_my, {
  id<-factor(id)
})
alt_my$sex<-factor(alt_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
alt_my$case <- 1

#MGUS CASES
alt_mgus <-filter(mgus_pat_wbloods, test_name=="ALT")
alt_mgus$flag<-ifelse(grepl("<", alt_mgus$result, ignore.case = T), 1, 0)
alt_mgus %>% tally(flag)

#check values <
alt_mgus %>% filter(flag==1) %>% select(result) %>% table()
alt_mgus$res<-alt_mgus$result
alt_mgus$res[alt_mgus$res == ""] <- NA
alt_mgus<-alt_mgus %>% drop_na(res)
alt_mgus$res <- sub("<", NA, alt_mgus$res)
alt_mgus$res <- as.numeric(alt_mgus$res)
summary(alt_mgus$res, na.rm=TRUE)
#1 value recorded as <2 recode to 2
alt_mgus$res[alt_mgus$result == "<2"] <- 2
# lowest recorded ALT value is 2, for those with values recorded as <9, impute random value between 2 and 9
alt_mgus <- alt_mgus %>% 
  mutate(m=ifelse(is.na(res), 1, 0))
alt_mgus %>% tally(m)
#create columns of 1s to sum over
alt_mgus$c<-1
#ipute imissing values as random value between 2 and 9
alt_mgus <- alt_mgus %>% 
  mutate(alt_im =
           case_when(
             m==0 ~ res,
             m==1 ~ floor(runif(sum(c), min=2, max=9))))
alt_mgus<- alt_mgus %>%
  select(id, alt_im, days_prior, age, sex, flag, dept) %>%
  rename (alt=alt_im,) %>% 
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
alt_mgus <- within(alt_mgus, {
  id<-factor(id)
})
alt_mgus$sex<-factor(alt_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
alt_mgus$case <- 2

#CONTROLS
alt_con <-filter(con_wbloods, test_name=="ALT")
alt_con %>% summarise(n_distinct(id))

alt_con$flag<-ifelse(grepl("<", alt_con$result, ignore.case = T), 1, 0)
alt_con %>% tally(flag)
#check values <
alt_con %>% filter(flag==1) %>% select(result) %>% table()
alt_con$res<-alt_con$result
alt_con$res[alt_con$res == ""] <- NA
alt_con<-alt_con %>% drop_na(res)
alt_con$res <- sub("<", NA, alt_con$res)
alt_con$res <- as.numeric(alt_con$res)
summary(alt_con$res, na.rm=TRUE)
#5 value recorded as <2 recode to 2
alt_con$res[alt_con$result == "<2"] <- 2
# lowest recorded ALT value is 2, for those with values recorded as <9, impute random value between 2 and 9
alt_con <- alt_con %>% 
  mutate(m=ifelse(is.na(res), 1, 0))
alt_con %>% tally(m)
#create columns of 1s to sum over
alt_con$c<-1
#ipute imissing values as random value between 2 and 9
alt_con <- alt_con %>% 
  mutate(alt_im =
           case_when(
             m==0 ~ res,
             m==1 ~ floor(runif(sum(c), min=2, max=9))))
alt_con<- alt_con %>%
  select(id, alt_im, days_prior, age, sex, flag, dept) %>%
  rename (alt=alt_im,) %>% 
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
alt_con <- within(alt_con, {
  id<-factor(id)
})
alt_con$sex<-factor(alt_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
alt_con$case <- 3

#append myeloma and MGUS 
alt<-rbind(alt_my, alt_mgus, alt_con)
#case as factor variable
alt$case<-factor(alt$case, levels=c('1', '2', '3'), 
                 labels=c('Myeloma', 'MGUS', 'Control'))
table(alt$case)
alt %>% group_by(case) %>% summarise(count=n_distinct(id))

#remove values >300
alt %>% subset(alt>1000) %>% count()
alt %>% subset(alt>1000) %>% group_by(case) %>% count()
alt %>% subset(alt>1000) %>% group_by(case) %>% distinct(id) %>% count()

alt<- alt %>% subset(alt<=1000)

#save as R file 
save(alt, file="N:\\DataAnalysis\\Trajectories\\ALT\\Alt_long.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Activated Partial Thromboplastin Clotting Time (APTT)

#MYELOMA CASES
aptt_my <-filter(pat_wbloods, test_name=="APTT")
#create flag variables for < and > values
aptt_my$flag1<-ifelse(grepl("<", aptt_my$result, ignore.case = T), 1, 0)
aptt_my$flag2<-ifelse(grepl(">240", aptt_my$result, ignore.case = T), 1, 0)
aptt_my$flag3<-ifelse(grepl(">250", aptt_my$result, ignore.case = T), 1, 0)
aptt_my$flag<-ifelse(grepl("<|>", aptt_my$result, ignore.case = T), 1, 0)
aptt_my %>% tally(flag)
aptt_my %>%  tally(flag1) 
aptt_my %>% tally(flag2)
aptt_my %>% tally(flag3)
table(subset(aptt_my, flag==1)$res)
aptt_my$res<-aptt_my$result
aptt_my$res[aptt_my$res == ""] <- NA
aptt_my<-aptt_my %>% drop_na(res)
aptt_my$res <- sub("<", NA, aptt_my$res)
aptt_my$res <- sub(">", NA, aptt_my$res)
aptt_my$res <- as.numeric(aptt_my$res)
#replace < and > values with actual value
aptt_my <- aptt_my %>% mutate(aptt_imp=
                                case_when(
                                  (flag==0) ~ res, 
                                  (flag==1 & (result=="<20"|result=="<20.0")) ~ 20,
                                  (flag==1 & result==">240") ~ 240,
                                  (flag==1 & (result==">250"|result==">250.0") ~250)))
aptt_my$atpp_imp <- as.numeric(aptt_my$aptt_imp)
aptt_my<- aptt_my %>%
  select(id, aptt_imp, days_prior, age, sex, dept) %>%
  rename (aptt=aptt_imp,) %>% 
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
aptt_my <- within(aptt_my, {
  id<-factor(id)
})
aptt_my$sex<-factor(aptt_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
aptt_my$case <- 1

#MGUS CASES
aptt_mgus <-filter(mgus_pat_wbloods, test_name=="APTT")
aptt_mgus$res<-aptt_mgus$result
aptt_mgus$res[aptt_mgus$res == ""] <- NA
aptt_mgus<-aptt_mgus %>% drop_na(res)
aptt_mgus$res <- as.numeric(aptt_mgus$res)
aptt_mgus<- aptt_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  rename (aptt=res) %>% 
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
aptt_mgus <- within(aptt_mgus, {
  id<-factor(id)
})
aptt_mgus$sex<-factor(aptt_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
aptt_mgus$case <- 2

#CONTROLS
aptt_con <-filter(con_wbloods, test_name=="APTT")
aptt_con %>% summarise(n_distinct(id))
aptt_con$flag<-ifelse(grepl("<|>", aptt_con$result, ignore.case = T), 1, 0)
aptt_con %>% tally(flag)
table(subset(aptt_con, flag==1)$res)
aptt_con$res<-aptt_con$result
aptt_con$res[aptt_con$res == ""] <- NA
aptt_con<-aptt_con %>% drop_na(res)
aptt_con$res <- sub("<", NA, aptt_con$res)
aptt_con$res <- sub(">", NA, aptt_con$res)
aptt_con$res <- as.numeric(aptt_con$res)
#replace < and > values with actual value
aptt_con <- aptt_con %>% mutate(aptt_imp=
                                  case_when(
                                    (flag==0) ~ res, 
                                    (flag==1 & (result=="<20"|result=="<20.0")) ~ 20,
                                    (flag==1 & result==">240") ~ 240,
                                    (flag==1 & result==">245") ~245,
                                    (flag==1 & (result==">250"|result==">250.0") ~250)))
aptt_con$atpp_imp <- as.numeric(aptt_con$aptt_imp)
aptt_con<- aptt_con %>%
  select(id, aptt_imp, days_prior, age, sex, dept) %>%
  rename (aptt=aptt_imp,) %>% 
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
aptt_con <- within(aptt_con, {
  id<-factor(id)
})
aptt_con$sex<-factor(aptt_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
aptt_con$case <- 3

#append myeloma and MGUS 
aptt<-rbind(aptt_my, aptt_mgus, aptt_con)
#case as factor variable
aptt$case<-factor(aptt$case, levels=c('1', '2', '3'), 
                  labels=c('Myeloma', 'MGUS', 'Control'))
table(aptt$case)

#save as R file 
save(aptt, file="N:\\DataAnalysis\\Trajectories\\APTT\\APTT_long.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#BASOPHIL

#MYELOMA CASES
bas_my <-filter(pat_wbloods, test_name=="Basophil count")
bas_my$res<-bas_my$result
bas_my$res[bas_my$res == ""] <- NA
bas_my$res <- as.numeric(bas_my$res)
bas_my <- bas_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (bas=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
bas_my <- within(bas_my, {
  id<-factor(id)
})
bas_my$sex<-factor(bas_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
bas_my$case<-1

#MGUS CASES
bas_mgus <-filter(mgus_pat_wbloods, test_name=="Basophil count")
bas_mgus$res<-bas_mgus$result
bas_mgus$res[bas_mgus$res == ""] <- NA
bas_mgus$res <- as.numeric(bas_mgus$res)
bas_mgus <- bas_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (bas=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
bas_mgus <- within(bas_mgus, {
  id<-factor(id)
})
bas_mgus$sex<-factor(bas_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
bas_mgus$case<-2

#CONTROLS
bas_con <-filter(con_wbloods, test_name=="Basophil count")
bas_con$res<-bas_con$result
bas_con$res[bas_con$res == ""] <- NA
bas_con$res <- as.numeric(bas_con$res)
bas_con <- bas_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (bas=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
bas_con <- within(bas_con, {
  id<-factor(id)
})
bas_con$sex<-factor(bas_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
bas_con$case<-3

#append myeloma and MGUS 
bas<-rbind(bas_my, bas_mgus, bas_con)
#case as factor variable
bas$case<-factor(bas$case, levels=c('1', '2', '3'), 
                 labels=c('Myeloma', 'MGUS', 'Control'))
table(bas$case)

#save as R file
save(bas, file="N:\\DataAnalysis\\Trajectories\\Basophil\\Bas_long.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#BETA 2 MICROGLOBULIN

#MYELOMA CASES
b2m_my<-filter(pat_wbloods, test_name=="Beta-2-microglobulin")
b2m_my$res<-b2m_my$result
b2m_my$res[b2m_my$res == ""] <- NA
b2m_my$res <- as.numeric(b2m_my$res)
b2m_my <- b2m_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (b2m=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
b2m_my <- within(b2m_my, {
  id<-factor(id)
})
b2m_my$sex<-factor(b2m_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
b2m_my$case<-1

#MGUS CASES
b2m_mgus <-filter(mgus_pat_wbloods, test_name=="Beta-2-microglobulin")
b2m_mgus$res<-b2m_mgus$result
b2m_mgus$res[b2m_mgus$res == ""] <- NA
b2m_mgus$res <- as.numeric(b2m_mgus$res)
b2m_mgus <- b2m_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (b2m=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
b2m_mgus <- within(b2m_mgus, {
  id<-factor(id)
})
b2m_mgus$sex<-factor(b2m_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
b2m_mgus$case<-2

#CONTROLS
b2m_con <-filter(con_wbloods, test_name=="Beta-2-microglobulin")
b2m_con$res<-b2m_con$result
b2m_con$res[b2m_con$res == ""] <- NA
b2m_con$res <- as.numeric(b2m_con$res)
b2m_con <- b2m_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (b2m=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
b2m_con <- within(b2m_con, {
  id<-factor(id)
})
b2m_con$sex<-factor(b2m_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
b2m_con$case<-3

#append myeloma and MGUS 
b2m<-rbind(b2m_my, b2m_mgus, b2m_con)
#case as factor variable
b2m$case<-factor(b2m$case, levels=c('1', '2', '3'), 
                 labels=c('Myeloma', 'MGUS', 'Control'))
table(b2m$case)

#save as R file
save(b2m, file="N:\\DataAnalysis\\Trajectories\\Beta2Microglobulin\\B2M_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#CALCIUM

#MYELOMA CASES
cal_my <-filter(pat_wbloods, test_name=="Calcium-adjusted")
cal_my$res<-cal_my$result
cal_my$res[cal_my$res == ""] <- NA
cal_my$res <- as.numeric(cal_my$res)
cal_my <- cal_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (cal=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>% 
  
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
cal_my <- within(cal_my, {
  id<-factor(id)
})
cal_my$sex<-factor(cal_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
cal_my$case<-1

#MGUS CASES
cal_mgus <-filter(mgus_pat_wbloods, test_name=="Calcium-adjusted")
cal_mgus$res<-cal_mgus$result
cal_mgus$res[cal_mgus$res == ""] <- NA
cal_mgus$res <- as.numeric(cal_mgus$res)
cal_mgus <- cal_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (cal=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>% 
  
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
cal_mgus <- within(cal_mgus, {
  id<-factor(id)
})
cal_mgus$sex<-factor(cal_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
cal_mgus$case<-2

#Conrols
cal_con <-filter(con_wbloods, test_name=="Calcium-adjusted")
cal_con$res<-cal_con$result
cal_con$res[cal_con$res == ""] <- NA
cal_con$res <- as.numeric(cal_con$res)
cal_con <- cal_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (cal=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>% 
  
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
cal_con <- within(cal_con, {
  id<-factor(id)
})
cal_con$sex<-factor(cal_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
cal_con$case<-3

#append myeloma and MGUS 
cal<-rbind(cal_my, cal_mgus, cal_con)
#case as factor variable
cal$case<-factor(cal$case, levels=c('1', '2', '3'), 
                 labels=c('Myeloma', 'MGUS', 'Control'))
table(cal$case)

#save as R file 
save(cal, file="N:\\DataAnalysis\\Trajectories\\Calcium\\Cal_long.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#CREATININE

#MYELOMA CASES
creat_my <-filter(pat_wbloods, test_name=="Creatinine")
creat_my$res<-creat_my$result
creat_my$res[creat_my$res == ""] <- NA
creat_my$res <- as.numeric(creat_my$res)
creat_my<- creat_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (creat=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>% 
  
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
creat_my <- within(creat_my, {
  id<-factor(id)
})
creat_my$sex<-factor(creat_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
creat_my$case<-1

#MGUS CASES
creat_mgus <-filter(mgus_pat_wbloods, test_name=="Creatinine")
creat_mgus$res<-creat_mgus$result
creat_mgus$res[creat_mgus$res == ""] <- NA
creat_mgus$res <- as.numeric(creat_mgus$res)
creat_mgus<- creat_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (creat=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>% 
  
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
creat_mgus <- within(creat_mgus, {
  id<-factor(id)
})
creat_mgus$sex<-factor(creat_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
creat_mgus$case<-2

#Controls
creat_con <-filter(con_wbloods, test_name=="Creatinine")
creat_con$res<-creat_con$result
creat_con$res[creat_con$res == ""] <- NA
creat_con$res <- as.numeric(creat_con$res)
creat_con<- creat_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (creat=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>% 
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
creat_con <- within(creat_con, {
  id<-factor(id)
})
creat_con$sex<-factor(creat_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
creat_con$case<-3

#append myeloma and MGUS 
creat<-rbind(creat_my, creat_mgus, creat_con)
#case as factor variable
creat$case<-factor(creat$case, levels=c('1', '2', '3'), 
                   labels=c('Myeloma', 'MGUS', 'Control'))
table(creat$case)

#save as R file 
save(creat, file="N:\\DataAnalysis\\Trajectories\\Creatinine\\Creat_long.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#C-Reactive Protein (CRP)
#RECODE values <5 AS 5

#MYELOMA CASES
crp_my <-filter(pat_wbloods, test_name=="CRP")
crp_my$res<-crp_my$result
crp_my$res[crp_my$res == ""] <- NA
crp_my$res[crp_my$res == "<5.0"] <- 5
crp_my$res <- as.numeric(crp_my$res)
crp_my<-crp_my %>% drop_na(res)
crp_my<- crp_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  rename (crp=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
crp_my <- within(crp_my, {
  id<-factor(id)
})
crp_my$sex<-factor(crp_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
crp_my$case<-1

#MGUS CASES
crp_mgus <-filter(mgus_pat_wbloods, test_name=="CRP")
crp_mgus$res<-crp_mgus$result
crp_mgus$res[crp_mgus$res == ""] <- NA
crp_mgus$res[crp_mgus$res == "<5.0"] <- 5
crp_mgus$res <- as.numeric(crp_mgus$res)
crp_mgus<-crp_mgus %>% drop_na(res)
crp_mgus <- crp_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  rename (crp=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
crp_mgus <- within(crp_mgus, {
  id<-factor(id)
})
crp_mgus$sex<-factor(crp_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
crp_mgus$case <- 2

#CONRTOLS
crp_con <-filter(con_wbloods, test_name=="CRP")
crp_con$res<-crp_con$result
crp_con$res[crp_con$res == ""] <- NA
crp_con$res[crp_con$res == "<4.0"] <- 5
crp_con$res[crp_con$res == "<5.0"] <- 5
crp_con$res <- as.numeric(crp_con$res)
crp_con<-crp_con %>% drop_na(res)
crp_con <- crp_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  rename (crp=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
crp_con <- within(crp_con, {
  id<-factor(id)
})
crp_con$sex<-factor(crp_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
crp_con$case <- 3

#append myeloma, MGUS and controls 
crp<-rbind(crp_my, crp_mgus, crp_con)
#case as factor variable
crp$case<-factor(crp$case, levels=c('1', '2', '3'), 
                 labels=c('Myeloma', 'MGUS', 'Control'))
table(crp$case)
crp %>% group_by(case) %>% summarise(count=n_distinct(id))

#save as R file
save(crp, file="N:\\DataAnalysis\\Trajectories\\CRP\\CRP_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#D-Dimer

#MYELOMA CASES
dd_my <-filter(pat_wbloods, test_name=="DDimer level")
table(dd_my$units)
#some coded as <150 (n=5), <200 (n=5), >1000 (n=3) these included in full range of recorded values so replace with value
dd_my$flag<-ifelse(grepl("<|>", dd_my$result, ignore.case = T), 1, 0)
dd_my %>% tally(flag)
table(subset(dd_my, flag==1)$result)
dd_my$res<-dd_my$result
dd_my$res[dd_my$res == ""] <- NA
dd_my<-dd_my %>% drop_na(res)
dd_my$res <- sub("<", NA, dd_my$res)
dd_my$res <- sub(">", NA, dd_my$res)
dd_my$res <- as.numeric(dd_my$res)
#replace < and > values with actual value
dd_my <- dd_my %>% mutate(dd_imp=
                            case_when(
                              (flag==0) ~ res, 
                              (flag==1 & result=="<150") ~ 150,
                              (flag==1 & result=="<200") ~ 200,
                              (flag==1 & result==">1000") ~1000))
dd_my$dd_imp <- as.numeric(dd_my$dd_imp)
dd_my<- dd_my %>%
  select(id, dd_imp, days_prior, age, sex, dept) %>%
  rename (dd=dd_imp,) %>% 
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
dd_my <- within(dd_my, {
  id<-factor(id)
})
dd_my$sex<-factor(dd_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
dd_my$case <- 1

#MGUS
dd_mgus<-filter(mgus_pat_wbloods, test_name=="DDimer level")
table(dd_mgus$units)
describe(dd_mgus$result)
#some coded as <150 (n=6), <200 (n=1), >1000 (n=1), >5000 (n=1) 
#these included in full range of recorded values so replace with value
dd_mgus$flag<-ifelse(grepl("<|>", dd_mgus$result, ignore.case = T), 1, 0)
dd_mgus %>% tally(flag)
table(subset(dd_mgus, flag==1)$result)
dd_mgus$res<-dd_mgus$result
dd_mgus$res[dd_mgus$res == ""] <- NA
dd_mgus<-dd_mgus %>% drop_na(res)
dd_mgus$res <- sub("<", NA, dd_mgus$res)
dd_mgus$res <- sub(">", NA, dd_mgus$res)
dd_mgus$res <- as.numeric(dd_mgus$res)
#replace < and > values with actual value
dd_mgus <- dd_mgus %>% mutate(dd_imp=
                                case_when(
                                  (flag==0) ~ res, 
                                  (flag==1 & result=="<150") ~ 150,
                                  (flag==1 & result=="<200") ~ 200,
                                  (flag==1 & result==">1000") ~1000,
                                  (flag==1 & result==">5000") ~5000))
dd_mgus$dd_imp <- as.numeric(dd_mgus$dd_imp)
dd_mgus<- dd_mgus %>%
  select(id, dd_imp, days_prior, age, sex, dept) %>%
  rename (dd=dd_imp,) %>% 
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
dd_mgus <- within(dd_mgus, {
  id<-factor(id)
})
dd_mgus$sex<-factor(dd_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
dd_mgus$case <- 2

#Control
dd_con<-filter(con_wbloods, test_name=="DDimer level")
table(dd_con$units)
describe(dd_con$result)
#some coded as <150 (n=9), <200 (n=4), >1000 (n=4), >5000 (n=1) 
#these included in full range of recorded values so replace with value
dd_con$flag<-ifelse(grepl("<|>", dd_con$result, ignore.case = T), 1, 0)
dd_con %>% tally(flag)
table(subset(dd_con, flag==1)$result)
dd_con$res<-dd_con$result
dd_con$res[dd_con$res == ""] <- NA
#replace 2 with values of 0 as missing 
table(subset(dd_con, result=="0")$result)
dd_con$res[dd_con$res == "0"] <- NA
dd_con<-dd_con %>% drop_na(res)
dd_con$res <- sub("<", NA, dd_con$res)
dd_con$res <- sub(">", NA, dd_con$res)
dd_con$res <- as.numeric(dd_con$res)
#replace < and > values with actual value
dd_con <- dd_con %>% mutate(dd_imp=
                              case_when(
                                (flag==0) ~ res, 
                                (flag==1 & result=="<150") ~ 150,
                                (flag==1 & result=="<200") ~ 200,
                                (flag==1 & result==">1000") ~1000))
dd_con$dd_imp <- as.numeric(dd_con$dd_imp)
dd_con<- dd_con %>%
  select(id, dd_imp, days_prior, age, sex, dept) %>%
  rename (dd=dd_imp,) %>% 
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
dd_con <- within(dd_con, {
  id<-factor(id)
})
dd_con$sex<-factor(dd_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
dd_con$case <- 3

#append myeloma and MGUS 
dd<-rbind(dd_my, dd_mgus, dd_con)
#case as factor variable
dd$case<-factor(dd$case, levels=c('1', '2', '3'), 
                labels=c('Myeloma', 'MGUS', 'Control'))
table(dd$case)

#save as R file 
save(dd, file="N:\\DataAnalysis\\Trajectories\\DDimer\\DDimer_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Estimated Glomerular Filtration Rate (eGFR)
#use categories

#MYELOMA CASES
egfr_my <-filter(pat_wbloods, test_name=="eGFR")
egfr_my$flag_gt<-ifelse(grepl(">", egfr_my$result, ignore.case = T), 1, 0)
egfr_my$res<-egfr_my$result
egfr_my$res[egfr_my$res == ""] <- NA
egfr_my<-egfr_my %>% drop_na(res)
egfr_my$res <- sub(">", NA, egfr_my$res)
egfr_my$res <- as.numeric(egfr_my$res)
egfr_my<-egfr_my %>%
  mutate(res2=recode(result, ">90" = "91"))
egfr_my$res2<-as.numeric(paste(egfr_my$res2))
egfr_my<-egfr_my %>%
  mutate(egfr_cat=cut(res2, breaks=c(0, 15, 30, 60, 90, 100),
                      right=FALSE, label=c("<15", "15-29", "30-59", "60-89", "90+")
  ))
egfr_my<- egfr_my %>%
  select(id, egfr_cat, days_prior, age, sex, dept) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
egfr_my <- within(egfr_my, {
  id<-factor(id)
})
egfr_my$sex<-factor(egfr_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
egfr_my$case<-1

#MGUS CASES
egfr_mgus <-filter(mgus_pat_wbloods, test_name=="eGFR")
egfr_mgus$flag_gt<-ifelse(grepl(">", egfr_mgus$result, ignore.case = T), 1, 0)
egfr_mgus$res<-egfr_mgus$result
egfr_mgus$res[egfr_mgus$res == ""] <- NA
egfr_mgus<-egfr_mgus %>% drop_na(res)
egfr_mgus$res <- sub(">", NA, egfr_mgus$res)
egfr_mgus$res <- as.numeric(egfr_mgus$res)
egfr_mgus<-egfr_mgus %>%
  mutate(res2=recode(result, ">90" = "91"))
egfr_mgus$res2<-as.numeric(paste(egfr_mgus$res2))
egfr_mgus<-egfr_mgus %>%
  mutate(egfr_cat=cut(res2, breaks=c(0, 15, 30, 60, 90, 100),
                      right=FALSE, label=c("<15", "15-29", "30-59", "60-89", "90+")
  ))
egfr_mgus<- egfr_mgus %>%
  select(id, egfr_cat, days_prior, age, sex, dept) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
egfr_mgus <- within(egfr_mgus, {
  id<-factor(id)
})
egfr_mgus$sex<-factor(egfr_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
egfr_mgus$case<-2

#CONTROLS
egfr_con <-filter(con_wbloods, test_name=="eGFR")
egfr_con$flag_gt<-ifelse(grepl(">", egfr_con$result, ignore.case = T), 1, 0)
egfr_con$res<-egfr_con$result
egfr_con$res[egfr_con$res == ""] <- NA
egfr_con<-egfr_con %>% drop_na(res)
egfr_con$res <- sub(">", NA, egfr_con$res)
egfr_con$res <- as.numeric(egfr_con$res)
egfr_con<-egfr_con %>%
  mutate(res2=recode(result, ">90" = "91"))
egfr_con$res2<-as.numeric(paste(egfr_con$res2))
egfr_con<-egfr_con %>%
  mutate(egfr_cat=cut(res2, breaks=c(0, 15, 30, 60, 90, 100),
                      right=FALSE, label=c("<15", "15-29", "30-59", "60-89", "90+")
  ))
egfr_con<- egfr_con %>%
  select(id, egfr_cat, days_prior, age, sex, dept) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
egfr_con <- within(egfr_con, {
  id<-factor(id)
})
egfr_con$sex<-factor(egfr_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
egfr_con$case<-3

#append myeloma, MGUS and controls 
egfr<-rbind(egfr_my, egfr_mgus, egfr_con)
#case as factor variable
egfr$case<-factor(egfr$case, levels=c('1', '2', '3'), 
                  labels=c('Myeloma', 'MGUS', 'Control'))
table(egfr$case)

#save as R file
save(egfr, file="N:\\DataAnalysis\\Trajectories\\eGFR\\eGFR_long.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EOSINOPHIL
#MYELOMA CASES
eos_my <-filter(pat_wbloods, test_name=="Eosinophil count")
eos_my$res<-eos_my$result
eos_my$res[eos_my$res == ""] <- NA
eos_my$res <- as.numeric(eos_my$res)
eos_my <- eos_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (eos=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
eos_my <- within(eos_my, {
  id<-factor(id)
})
eos_my$sex<-factor(eos_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
eos_my$case<-1

#MGUS CASES
eos_mgus <-filter(mgus_pat_wbloods, test_name=="Eosinophil count")
eos_mgus$res<-eos_mgus$result
eos_mgus$res[eos_mgus$res == ""] <- NA
eos_mgus$res <- as.numeric(eos_mgus$res)
eos_mgus <- eos_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (eos=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
eos_mgus <- within(eos_mgus, {
  id<-factor(id)
})
eos_mgus$sex<-factor(eos_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
eos_mgus$case<-2

#CONTROLS
eos_con <-filter(con_wbloods, test_name=="Eosinophil count")
eos_con$res<-eos_con$result
eos_con$res[eos_con$res == ""] <- NA
eos_con$res <- as.numeric(eos_con$res)
eos_con <- eos_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (eos=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
eos_con <- within(eos_con, {
  id<-factor(id)
})
eos_con$sex<-factor(eos_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
eos_con$case<-3

#append myeloma, MGUS and controls 
eos<-rbind(eos_my, eos_mgus, eos_con)
#case as factor variable
eos$case<-factor(eos$case, levels=c('1', '2', '3'), 
                 labels=c('Myeloma', 'MGUS', 'Control'))
table(eos$case)

#save as R file 
save(eos, file="N:\\DataAnalysis\\Trajectories\\Eosinophil\\Eos_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Erythrocyte Sedimentation Rate (ESR) 

#MYELOMA CASES
esr_my <-filter(pat_wbloods, test_name=="Erythrocyte sed rate")
esr_my$res<-esr_my$result
esr_my$res[esr_my$res == ""] <- NA
esr_my$res <- as.numeric(esr_my$res)
esr_my <- esr_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (esr=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
esr_my <- within(esr_my, {
  id<-factor(id)
})
esr_my$sex<-factor(esr_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
esr_my$case<-1

#MGUS CASES
esr_mgus <-filter(mgus_pat_wbloods, test_name=="Erythrocyte sed rate")
esr_mgus$res<-esr_mgus$result
esr_mgus$res[esr_mgus$res == ""] <- NA
esr_mgus$res <- as.numeric(esr_mgus$res)
esr_mgus <- esr_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (esr=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
esr_mgus <- within(esr_mgus, {
  id<-factor(id)
})
esr_mgus$sex<-factor(esr_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
esr_mgus$ case <- 2

#CONTROL
esr_con <-filter(con_wbloods, test_name=="Erythrocyte sed rate")
esr_con$res<-esr_con$result
esr_con$res[esr_con$res == ""] <- NA
esr_con$res <- as.numeric(esr_con$res)
esr_con <- esr_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (esr=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
esr_con <- within(esr_con, {
  id<-factor(id)
})
esr_con$sex<-factor(esr_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
esr_con$ case <- 3

#append myeloma, MGUS and controls 
esr<-rbind(esr_my, esr_mgus, esr_con)
#case as factor variable
esr$case<-factor(esr$case, levels=c('1', '2', '3'), 
                 labels=c('Myeloma', 'MGUS', 'Control'))
table(esr$case)

#save as R file
save(esr, file="N:\\DataAnalysis\\Trajectories\\ESR\\ESR_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#FIBRINOGEN
# TWO GROUPS <4.5, >=4.5

#MYELOMA
df_my <-filter(pat_wbloods, test_name=="Derived Fibrinogen"|test_name=="Derived fibrinogen"|
                 test_name=="Fibrinogen Level"|test_name=="Fibrinogen level")
df_my$flag1<-ifelse(grepl(">", df_my$result, ignore.case = T), 1, 0)
df_my %>% tally(flag1)
table(subset(df_my, flag1==1)$result)
df_my$res<-df_my$result
df_my$res[df_my$res == ""] <- NA
df_my<-df_my %>% drop_na(res)
df_my$res <- sub(">", NA, df_my$res)
df_my$res <- as.numeric(df_my$res)
df_my<- df_my %>%
  mutate(res2=recode(result, ">4.5" = "4.55"))
df_my$res2<-as.numeric(paste(df_my$res2))
df_my<-df_my %>%
  mutate(df_cat=cut(res2, breaks=c(0, 4.5, 20),
                    right=FALSE, label=c("<4.5", ">4.5")
  ))
df_my<- df_my %>%
  select(id, df_cat, days_prior, age, sex, dept) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
df_my <- within(df_my, {
  id<-factor(id)
})
df_my$sex<-factor(df_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
df_my$case <- 1

#MGUS
df_mgus <-filter(mgus_pat_wbloods, test_name=="Derived Fibrinogen"|test_name=="Derived fibrinogen"|
                   test_name=="Fibrinogen Level"|test_name=="Fibrinogen level")
df_mgus$flag1<-ifelse(grepl(">", df_mgus$result, ignore.case = T), 1, 0)
df_mgus %>% tally(flag1)
table(subset(df_mgus, flag1==1)$result)
df_mgus$res<-df_mgus$result
df_mgus$res[df_mgus$res == ""] <- NA
df_mgus<-df_mgus %>% drop_na(res)
df_mgus$res <- sub(">", NA, df_mgus$res)
df_mgus$res <- as.numeric(df_mgus$res)
df_mgus<- df_mgus %>%
  mutate(res2=recode(result, ">4.5" = "4.55"))
df_mgus$res2<-as.numeric(paste(df_mgus$res2))
df_mgus<-df_mgus %>%
  mutate(df_cat=cut(res2, breaks=c(0, 4.5, 20),
                    right=FALSE, label=c("<4.5", ">4.5")
  ))
df_mgus<- df_mgus %>%
  select(id, df_cat, days_prior, age, sex, dept) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
df_mgus <- within(df_mgus, {
  id<-factor(id)
})
df_mgus$sex<-factor(df_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
df_mgus$case <- 2

#CONTROLS
df_con <-filter(con_wbloods, test_name=="Derived Fibrinogen"|test_name=="Derived fibrinogen"|
                  test_name=="Fibrinogen Level"|test_name=="Fibrinogen level")
df_con$flag1<-ifelse(grepl(">", df_con$result, ignore.case = T), 1, 0)
df_con %>% tally(flag1)
table(subset(df_con, flag1==1)$result)
df_con$res<-df_con$result
df_con$res[df_con$res == ""] <- NA
df_con<-df_con %>% drop_na(res)
df_con$res <- sub(">", NA, df_con$res)
df_con$res <- as.numeric(df_con$res)
df_con<- df_con %>%
  mutate(res2=recode(result, ">4.5" = "4.55"))
df_con$res2<-as.numeric(paste(df_con$res2))
df_con<-df_con %>%
  mutate(df_cat=cut(res2, breaks=c(0, 4.5, 20),
                    right=FALSE, label=c("<4.5", ">4.5")
  ))
df_con<- df_con %>%
  select(id, df_cat, days_prior, age, sex, dept) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
df_con <- within(df_con, {
  id<-factor(id)
})
df_con$sex<-factor(df_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
df_con$case <- 3

#append myeloma, MGUS and controls 
df<-rbind(df_my, df_mgus, df_con)
#case as factor variable
df$case<-factor(df$case, levels=c('1', '2', '3'), 
                labels=c('Myeloma', 'MGUS', 'Control'))
table(df$case)
table(df$case, df$df_cat)
#save as R file
save(df, file="N:\\DataAnalysis\\Trajectories\\Fibrinogen\\DF_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#HAEMOGLOBIN
#recode so all the same units

#MYELOMA CASES
haem_my<- filter(pat_wbloods, test_name=="Haemoglobin")
haem_my$res<-haem_my$result
haem_my$res[haem_my$res == ""] <- NA
haem_my$res <- as.numeric(haem_my$res)
haem_my<-haem_my %>% drop_na(res)
haem_my$flag_h <-ifelse((haem_my$units=="g/dl"|haem_my$units=="g/dL"), 1, 0)
haem_my<- haem_my %>%
  mutate(res=case_when(flag_h==0 ~  res, 
                       flag_h==1 ~ res*10,))
haem_my<- haem_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (haem=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>% 
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
haem_my <- within(haem_my, {
  id<-factor(id)
})
haem_my$sex<-factor(haem_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
haem_my$case<-1

#MGUS CASES
haem_mgus<- filter(mgus_pat_wbloods, test_name=="Haemoglobin")
haem_mgus$res<-haem_mgus$result
haem_mgus$res[haem_mgus$res == ""] <- NA
haem_mgus$res <- as.numeric(haem_mgus$res)
haem_mgus<-haem_mgus %>% drop_na(res)
haem_mgus$flag_h <-ifelse((haem_mgus$units=="g/dl"|haem_mgus$units=="g/dL"), 1, 0)
haem_mgus<- haem_mgus %>%
  mutate(res=case_when(flag_h==0 ~  res, 
                       flag_h==1 ~ res*10,))
haem_mgus<- haem_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (haem=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>% 
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
haem_mgus <- within(haem_mgus, {
  id<-factor(id)
})
haem_mgus$sex<-factor(haem_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
haem_mgus$case<-2

#CONTROLs
haem_con<- filter(con_wbloods, test_name=="Haemoglobin")
haem_con$res<-haem_con$result
haem_con$res[haem_con$res == ""] <- NA
haem_con$res <- as.numeric(haem_con$res)
haem_con<-haem_con %>% drop_na(res)
haem_con$flag_h <-ifelse((haem_con$units=="g/dl"|haem_con$units=="g/dL"), 1, 0)
haem_con<- haem_con %>%
  mutate(res=case_when(flag_h==0 ~  res, 
                       flag_h==1 ~ res*10,))
haem_con<- haem_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (haem=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>% 
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
haem_con <- within(haem_con, {
  id<-factor(id)
})
haem_con$sex<-factor(haem_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
haem_con$case<-3

#append myeloma, MGUS and controls 
haem<-rbind(haem_my, haem_mgus, haem_con)
#case as factor variable
haem$case<-factor(haem$case, levels=c('1', '2', '3'), 
                  labels=c('Myeloma', 'MGUS','Control'))
table(haem$case)

#save as R file to use in further analysis
save(haem, file="N:\\DataAnalysis\\Trajectories\\Haemoglobin\\Haem_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# % HYPOCHROMIC CELLS

#MYELOMA CASES
hc_my <-filter(pat_wbloods, test_name=="% Hypochromic cells")
hc_my$res<-hc_my$result
hc_my$res[hc_my$res == ""] <- NA
hc_my$res <- as.numeric(hc_my$res)
hc_my <- hc_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (hc=res)%>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
hc_my <- within(hc_my, {
  id<-factor(id)
})
hc_my$sex<-factor(hc_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
hc_my$case<-1

#MGUS CASES
hc_mgus <-filter(mgus_pat_wbloods, test_name=="% Hypochromic cells")
hc_mgus$res<-hc_mgus$result
hc_mgus$res[hc_mgus$res == ""] <- NA
hc_mgus$res <- as.numeric(hc_mgus$res)
hc_mgus <- hc_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (hc=res)%>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
hc_mgus <- within(hc_mgus, {
  id<-factor(id)
})
hc_mgus$sex<-factor(hc_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
hc_mgus$case <- 2

#CONTROLS
hc_con <-filter(con_wbloods, test_name=="% Hypochromic cells")
hc_con$res<-hc_con$result
hc_con$res[hc_con$res == ""] <- NA
hc_con$res <- as.numeric(hc_con$res)
hc_con <- hc_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (hc=res)%>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
hc_con <- within(hc_con, {
  id<-factor(id)
})
hc_con$sex<-factor(hc_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
hc_con$case <- 3

#append myeloma, MGUS and controls 
hc<-rbind(hc_my, hc_mgus, hc_con)
#case as factor variable
hc$case<-factor(hc$case, levels=c('1', '2', '3'),
                labels=c('Myeloma', 'MGUS', 'Control'))
table(hc$case)
#save as R 
save(hc, file="N:\\DataAnalysis\\Trajectories\\Hypochromic\\HC_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# International Normalized Ratio (INR)

#MYELOMA CASES
inr_my <-filter(pat_wbloods, test_name=="INR")
inr_my$res<-inr_my$result
inr_my$res[inr_my$res == ""] <- NA
inr_my$res <- as.numeric(inr_my$res)
inr_my <- inr_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (inr=res)  %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
inr_my <- within(inr_my, {
  id<-factor(id)
})
inr_my$sex<-factor(inr_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
inr_my$case<-1

#MGUS CASES
inr_mgus <-filter(mgus_pat_wbloods, test_name=="INR")
inr_mgus$res<-inr_mgus$result
inr_mgus$res[inr_mgus$res == ""] <- NA
inr_mgus$res <- as.numeric(inr_mgus$res)
inr_mgus <- inr_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (inr=res)  %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
inr_mgus <- within(inr_mgus, {
  id<-factor(id)
})
inr_mgus$sex<-factor(inr_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
inr_mgus$ case <- 2

#CONTROLS
inr_con <-filter(con_wbloods, test_name=="INR")
inr_con$res<-inr_con$result
inr_con$res[inr_con$res == ""] <- NA
inr_con$res <- as.numeric(inr_con$res)
inr_con <- inr_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (inr=res)  %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
inr_con <- within(inr_con, {
  id<-factor(id)
})
inr_con$sex<-factor(inr_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
inr_con$ case <- 3

#append myeloma, MGUS and controls 
inr<-rbind(inr_my, inr_mgus, inr_con)
#case as factor variable
inr$case<-factor(inr$case, levels=c('1', '2', '3'), 
                 labels=c('Myeloma', 'MGUS', 'Control'))
table(inr$case)
#save as R file
save(inr, file="N:\\DataAnalysis\\Trajectories\\INR\\INR_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#FREE KAPPA LIGHT CELLS

#MYELOMA CASES
fklc_my<-filter(pat_wbloods, test_name=="Free Kappa LC (Siem)")
fklc_my$res<-fklc_my$result
fklc_my$res[fklc_my$res == ""] <- NA
fklc_my$res <- as.numeric(fklc_my$res)
fklc_my <- fklc_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (fklc=res) %>% 
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))

fklc_my <- within(fklc_my, {
  id<-factor(id)
})
fklc_my$sex<-factor(fklc_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
fklc_my$case<-1

#MGUS CASES
fklc_mgus<-filter(mgus_pat_wbloods, test_name=="Free Kappa LC (Siem)")
fklc_mgus$res<-fklc_mgus$result
fklc_mgus$res[fklc_mgus$res == ""] <- NA
fklc_mgus$res <- as.numeric(fklc_mgus$res)
fklc_mgus <- fklc_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (fklc=res) %>% 
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))

fklc_mgus <- within(fklc_mgus, {
  id<-factor(id)
})
fklc_mgus$sex<-factor(fklc_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
fklc_mgus$case<-2

#CONTROLS
fklc_con<-filter(con_wbloods, test_name=="Free Kappa LC (Siem)")
fklc_con$res<-fklc_con$result
fklc_con$res[fklc_con$res == ""] <- NA
fklc_con$res <- as.numeric(fklc_con$res)
fklc_con <- fklc_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (fklc=res) %>% 
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))

fklc_con <- within(fklc_con, {
  id<-factor(id)
})
fklc_con$sex<-factor(fklc_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
fklc_con$case<-3

#append myeloma, MGUS and controls 
fklc<-rbind(fklc_my, fklc_mgus, fklc_con)
#case as factor variable
fklc$case<-factor(fklc$case, levels=c('1', '2', '3'), 
                  labels=c('Myeloma', 'MGUS', 'Control'))
table(fklc$case)
#save as R file 
save(fklc, file="N:\\DataAnalysis\\Trajectories\\KappaLC\\FKLC_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#FREE LAMBDA LIGHT CHAINS

#MYELOMA CASES
fllc_my<-filter(pat_wbloods, test_name=="Free Lambda LC(Siem)")
fllc_my$res<-fllc_my$result
fllc_my$res[fllc_my$res == ""] <- NA
fllc_my$res <- as.numeric(fllc_my$res)
fllc_my <- fllc_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (fllc=res) %>% 
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))

fllc_my <- within(fllc_my, {
  id<-factor(id)
})
fllc_my$sex<-factor(fllc_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
fllc_my$case<-1

#MGUS CASES
fllc_mgus<-filter(mgus_pat_wbloods, test_name=="Free Lambda LC(Siem)")
fllc_mgus$res<-fllc_mgus$result
fllc_mgus$res[fllc_mgus$res == ""] <- NA
fllc_mgus$res <- as.numeric(fllc_mgus$res)
fllc_mgus <- fllc_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (fllc=res) %>% 
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))

fllc_mgus <- within(fllc_mgus, {
  id<-factor(id)
})
fllc_mgus$sex<-factor(fllc_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
fllc_mgus$case<-2

#CONTROLS
fllc_con<-filter(con_wbloods, test_name=="Free Lambda LC(Siem)")
fllc_con$res<-fllc_con$result
fllc_con$res[fllc_con$res == ""] <- NA
fllc_con$res <- as.numeric(fllc_con$res)
fllc_con <- fllc_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (fllc=res) %>% 
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))

fllc_con <- within(fllc_con, {
  id<-factor(id)
})
fllc_con$sex<-factor(fllc_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
fllc_con$case<-3

#append myeloma, MGUS and controls 
fllc<-rbind(fllc_my, fllc_mgus, fllc_con)
#case as factor variable
fllc$case<-factor(fllc$case, levels=c('1', '2', '3'), 
                  labels=c('Myeloma', 'MGUS', 'Control'))
table(fllc$case)
#save as R file 
save(fllc, file="N:\\DataAnalysis\\Trajectories\\LambdaLC\\FLLC_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Lactate dehydrogenase (LDH) 

#MYELOMA CASES
ldh_my <-filter(pat_wbloods, test_name=="LDH")
ldh_my$res<-ldh_my$result
ldh_my$res[ldh_my$res == ""] <- NA
ldh_my$res <- as.numeric(ldh_my$res)
ldh_my <- ldh_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (ldh=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
ldh_my <- within(ldh_my, {
  id<-factor(id)
})
ldh_my$sex<-factor(ldh_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
ldh_my$case<-1

#MGUS CASES
ldh_mgus <-filter(mgus_pat_wbloods, test_name=="LDH")
ldh_mgus$res<-ldh_mgus$result
ldh_mgus$res[ldh_mgus$res == ""] <- NA
ldh_mgus$res <- as.numeric(ldh_mgus$res)
ldh_mgus <- ldh_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (ldh=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
ldh_mgus <- within(ldh_mgus, {
  id<-factor(id)
})
ldh_mgus$sex<-factor(ldh_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
ldh_mgus$case<-2

#CONTROLS
ldh_con <-filter(con_wbloods, test_name=="LDH")
ldh_con$res<-ldh_con$result
ldh_con$res[ldh_con$res == ""] <- NA
ldh_con$res <- as.numeric(ldh_con$res)
ldh_con <- ldh_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (ldh=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
ldh_con <- within(ldh_con, {
  id<-factor(id)
})
ldh_con$sex<-factor(ldh_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
ldh_con$case<-3

#append myeloma, MGUS and controls 
ldh<-rbind(ldh_my, ldh_mgus, ldh_con)
#case as factor variable
ldh$case<-factor(ldh$case, levels=c('1', '2', '3'),
                 labels=c('Myeloma', 'MGUS', 'Control'))
table(ldh$case)
#save as R file
save(ldh, file="N:\\DataAnalysis\\Trajectories\\LDH\\LDH_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LYMPHOCYTES

#MYELOMA CASES
lymph_my <-filter(pat_wbloods, test_name=="Lymphocyte count")
lymph_my$res<-lymph_my$result
lymph_my$res[lymph_my$res == ""] <- NA
lymph_my$res <- as.numeric(lymph_my$res)
lymph_my <- lymph_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (lymph=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
lymph_my <- within(lymph_my, {
  id<-factor(id)
})
lymph_my$sex<-factor(lymph_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
lymph_my$case<-1

#MGUS CASES
lymph_mgus <-filter(mgus_pat_wbloods, test_name=="Lymphocyte count")
lymph_mgus$res<-lymph_mgus$result
lymph_mgus$res[lymph_mgus$res == ""] <- NA
lymph_mgus$res <- as.numeric(lymph_mgus$res)
lymph_mgus <- lymph_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (lymph=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
lymph_mgus <- within(lymph_mgus, {
  id<-factor(id)
})
lymph_mgus$sex<-factor(lymph_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
lymph_mgus$case<-2

#CONTROLS
lymph_con <-filter(con_wbloods, test_name=="Lymphocyte count")
lymph_con$res<-lymph_con$result
lymph_con$res[lymph_con$res == ""] <- NA
lymph_con$res <- as.numeric(lymph_con$res)
lymph_con <- lymph_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (lymph=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
lymph_con <- within(lymph_con, {
  id<-factor(id)
})
lymph_con$sex<-factor(lymph_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
lymph_con$case<-3

#append myeloma, MGUS and controls 
lymph<-rbind(lymph_my, lymph_mgus, lymph_con)
#case as factor variable
lymph$case<-factor(lymph$case, levels=c('1', '2', '3'), 
                   labels=c('Myeloma', 'MGUS', 'Control'))
table(lymph$case)

#remove values >10
lymph %>% subset(lymph>10) %>% count()
lymph %>% subset(lymph>10) %>% group_by(case) %>% count()
lymph %>% subset(lymph>10) %>% group_by(case) %>% distinct(id) %>% count()

lymph<- lymph %>% subset(lymph<=10)

#save as R file 
save(lymph, file="N:\\DataAnalysis\\Trajectories\\Lymphocyte\\Lymph_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Mean Cell Volume (MCV)

#MYELOMA CASES
mcv_my <-filter(pat_wbloods, test_name=="Mean cell volume MCV")
mcv_my$res<-mcv_my$result
mcv_my$res[mcv_my$res == ""] <- NA
mcv_my$res <- as.numeric(mcv_my$res)
mcv_my <- mcv_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (mcv=res)%>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+")))%>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
mcv_my <- within(mcv_my, {
  id<-factor(id)
})
mcv_my$sex<-factor(mcv_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
mcv_my$case <- 1 

#MGUS CASES
mcv_mgus <-filter(mgus_pat_wbloods, test_name=="Mean cell volume MCV")
mcv_mgus$res<-mcv_mgus$result
mcv_mgus$res[mcv_mgus$res == ""] <- NA
mcv_mgus$res <- as.numeric(mcv_mgus$res)
mcv_mgus <- mcv_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (mcv=res)%>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+")))%>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
mcv_mgus <- within(mcv_mgus, {
  id<-factor(id)
})
mcv_mgus$sex<-factor(mcv_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
mcv_mgus$case <- 2

#CONTROLS
mcv_con <-filter(con_wbloods, test_name=="Mean cell volume MCV")
mcv_con$res<-mcv_con$result
mcv_con$res[mcv_con$res == ""] <- NA
mcv_con$res <- as.numeric(mcv_con$res)
mcv_con <- mcv_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (mcv=res)%>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+")))%>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
mcv_con <- within(mcv_con, {
  id<-factor(id)
})
mcv_con$sex<-factor(mcv_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
mcv_con$case <- 3

#append myeloma, MGUS and control 
mcv<-rbind(mcv_my, mcv_mgus, mcv_con)
#case as factor variable
mcv$case<-factor(mcv$case, levels=c('1', '2', '3'), 
                 labels=c('Myeloma', 'MGUS', 'Control'))
table(mcv$case)
#save as R file 
save(mcv, file="N:\\DataAnalysis\\Trajectories\\MCV\\MCV_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Mean CORPUSC HB 

#MYELOMA CASES
mch_my <-filter(pat_wbloods, test_name=="Mean corpusc Hb MCH")
mch_my$res<-mch_my$result
mch_my$res[mch_my$res == ""] <- NA
mch_my$res <- as.numeric(mch_my$res)
mch_my <- mch_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (mch=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+")))%>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
mch_my <- within(mch_my, {
  id<-factor(id)
})
mch_my$sex<-factor(mch_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
mch_my$case <- 1

#MGUS CASES
mch_mgus <-filter(mgus_pat_wbloods, test_name=="Mean corpusc Hb MCH")
mch_mgus$res<-mch_mgus$result
mch_mgus$res[mch_mgus$res == ""] <- NA
mch_mgus$res <- as.numeric(mch_mgus$res)
mch_mgus <- mch_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (mch=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+")))%>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
mch_mgus <- within(mch_mgus, {
  id<-factor(id)
})
mch_mgus$sex<-factor(mch_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
mch_mgus$case <- 2

#CONTROLS
mch_con <-filter(con_wbloods, test_name=="Mean corpusc Hb MCH")
mch_con$res<-mch_con$result
mch_con$res[mch_con$res == ""] <- NA
mch_con$res <- as.numeric(mch_con$res)
mch_con <- mch_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (mch=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+")))%>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
mch_con <- within(mch_con, {
  id<-factor(id)
})
mch_con$sex<-factor(mch_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
mch_con$case <- 3

#append myeloma and MGUS 
mch<-rbind(mch_my, mch_mgus, mch_con)
#case as factor variable
mch$case<-factor(mch$case, levels=c('1', '2', '3'), 
                 labels=c('Myeloma', 'MGUS', 'Control'))
table(mch$case)
#save as R file 
save(mch, file="N:\\DataAnalysis\\Trajectories\\MeanCorpusc\\MCH_long.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#MONOCYTES

#MYELOMA CASES
mono_my <-filter(pat_wbloods, test_name=="Monocyte count")
mono_my$res<-mono_my$result
mono_my$res[mono_my$res == ""] <- NA
mono_my$res <- as.numeric(mono_my$res)
mono_my <- mono_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (mono=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+")))%>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
mono_my <- within(mono_my, {
  id<-factor(id)
})
mono_my$sex<-factor(mono_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
mono_my$case <- 1

#MGUS CASES
mono_mgus <-filter(mgus_pat_wbloods, test_name=="Monocyte count")
mono_mgus$res<-mono_mgus$result
mono_mgus$res[mono_mgus$res == ""] <- NA
mono_mgus$res <- as.numeric(mono_mgus$res)
mono_mgus <- mono_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (mono=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+")))%>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
mono_mgus <- within(mono_mgus, {
  id<-factor(id)
})
mono_mgus$sex<-factor(mono_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
mono_mgus$case <- 2

#CONTROLS
mono_cons <-filter(con_wbloods, test_name=="Monocyte count")
mono_cons$res<-mono_cons$result
mono_cons$res[mono_cons$res == ""] <- NA
mono_cons$res <- as.numeric(mono_cons$res)
mono_cons <- mono_cons %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (mono=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+")))%>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
mono_cons <- within(mono_cons, {
  id<-factor(id)
})
mono_cons$sex<-factor(mono_cons$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
mono_cons$case <- 3

#append myeloma, MGUS and control 
mono<-rbind(mono_my, mono_mgus, mono_cons)
#case as factor variable
mono$case<-factor(mono$case, levels=c('1', '2', '3'), 
                  labels=c('Myeloma', 'MGUS', 'Control'))
table(mono$case)

#remove values >3
mono %>% subset(mono>3) %>% count()
mono %>% subset(mono>3) %>% group_by(case) %>% count()
mono %>% subset(mono>3) %>% group_by(case) %>% distinct(id) %>% count()

mono<- mono %>% subset(mono<=3)

#save as R file 
save(mono, file="N:\\DataAnalysis\\Trajectories\\Monocyte\\Mono_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#NEURTOPHILS

#MYELOMA CASES
neut_my <-filter(pat_wbloods, test_name=="Neutrophil count")
neut_my$res<-neut_my$result
neut_my$res[neut_my$res == ""] <- NA
neut_my$res <- as.numeric(neut_my$res)
neut_my <- neut_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (neut=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+")))%>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
neut_my <- within(neut_my, {
  id<-factor(id)
})
neut_my$sex<-factor(neut_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
neut_my$case <- 1

#MGUS CASES
neut_mgus <-filter(mgus_pat_wbloods, test_name=="Neutrophil count")
neut_mgus$res<-neut_mgus$result
neut_mgus$res[neut_mgus$res == ""] <- NA
neut_mgus$res <- as.numeric(neut_mgus$res)
neut_mgus <- neut_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (neut=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+")))%>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
neut_mgus <- within(neut_mgus, {
  id<-factor(id)
})
neut_mgus$sex<-factor(neut_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
neut_mgus$case <- 2

#CONTROLS
neut_con <-filter(con_wbloods, test_name=="Neutrophil count")
neut_con$res<-neut_con$result
neut_con$res[neut_con$res == ""] <- NA
neut_con$res <- as.numeric(neut_con$res)
neut_con <- neut_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (neut=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+")))%>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
neut_con <- within(neut_con, {
  id<-factor(id)
})
neut_con$sex<-factor(neut_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
neut_con$case <- 3

#append myeloma, MGUS and controsl 
neut<-rbind(neut_my, neut_mgus, neut_con)
#case as factor variable
neut$case<-factor(neut$case, levels=c('1', '2', '3'), 
                  labels=c('Myeloma', 'MGUS', 'Control'))
table(neut$case)

#remove values >30
neut %>% subset(neut>30) %>% count()
neut %>% subset(neut>30) %>% group_by(case) %>% count()
neut %>% subset(neut>30) %>% group_by(case) %>% distinct(id) %>% count()

neut<- neut %>% subset(neut<=30)

#save as R file
save(neut, file="N:\\DataAnalysis\\Trajectories\\Neutrophils\\Neut_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#PACKED CELL VOLUME

#MYELOMA CASES
pcv_my <-filter(pat_wbloods, test_name=="Packed cell volume")
pcv_my$res<-pcv_my$result
pcv_my$res[pcv_my$res == ""] <- NA
pcv_my$res <- as.numeric(pcv_my$res)
pcv_my <- pcv_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (pcv=res)%>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
pcv_my <- within(pcv_my, {
  id<-factor(id)
})
pcv_my$sex<-factor(pcv_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
pcv_my$case <- 1

#MGUS CASES
pcv_mgus <-filter(mgus_pat_wbloods, test_name=="Packed cell volume")
pcv_mgus$res<-pcv_mgus$result
pcv_mgus$res[pcv_mgus$res == ""] <- NA
pcv_mgus$res <- as.numeric(pcv_mgus$res)
pcv_mgus <- pcv_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (pcv=res)%>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
pcv_mgus <- within(pcv_mgus, {
  id<-factor(id)
})
pcv_mgus$sex<-factor(pcv_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
pcv_mgus$case <- 2

#CONTROLS 
pcv_con <-filter(con_wbloods, test_name=="Packed cell volume")
pcv_con$res<-pcv_con$result
pcv_con$res[pcv_con$res == ""] <- NA
pcv_con$res <- as.numeric(pcv_con$res)
pcv_con <- pcv_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (pcv=res)%>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
pcv_con <- within(pcv_con, {
  id<-factor(id)
})
pcv_con$sex<-factor(pcv_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
pcv_con$case <- 3

#append myeloma, MGUS and controls 
pcv<-rbind(pcv_my, pcv_mgus, pcv_con)
#case as factor variable
pcv$case<-factor(pcv$case, levels=c('1', '2', '3'), 
                 labels=c('Myeloma', 'MGUS', 'Control'))
table(pcv$case)
#save as R file
save(pcv, file="N:\\DataAnalysis\\Trajectories\\PackedCellVolume\\PCV_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#PARAPROTEINS

#MYELOMA CASES
para_my<-filter(pat_wbloods, test_name=="Paraprotein")
para_my$res<-para_my$result
para_my$res[para_my$res == ""] <- NA
para_my$res <- as.numeric(para_my$res)
para_my <- para_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (para=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
para_my <- within(para_my, {
  id<-factor(id)
})
para_my$sex<-factor(para_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
para_my$case <- 1

#MGUS CASES
para_mgus<-filter(mgus_pat_wbloods, test_name=="Paraprotein")
para_mgus$res<-para_mgus$result
para_mgus$res[para_mgus$res == ""] <- NA
para_mgus$res <- as.numeric(para_mgus$res)
para_mgus <- para_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (para=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
para_mgus <- within(para_mgus, {
  id<-factor(id)
})
para_mgus$sex<-factor(para_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
para_mgus$case <- 2

#CONTROLS
para_con<-filter(con_wbloods, test_name=="Paraprotein")
para_con$res<-para_con$result
para_con$res[para_con$res == ""] <- NA
para_con$res <- as.numeric(para_con$res)
para_con <- para_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (para=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
para_con <- within(para_con, {
  id<-factor(id)
})
para_con$sex<-factor(para_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
para_con$case <- 3

#append myeloma, MGUS and controls 
para<-rbind(para_my, para_mgus, para_con)
#case as factor variable
para$case<-factor(para$case, levels=c('1', '2', '3'),
                  labels=c('Myeloma', 'MGUS', 'Control'))
table(para$case)
#save as R file
save(para, file="N:\\DataAnalysis\\Trajectories\\Paraprotein\\Para_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#PLATELETS 

#MYELOMA CASES
plt_my <-filter(pat_wbloods, test_name=="Platelets")
plt_my$res<-plt_my$result
plt_my$res[plt_my$res == ""] <- NA
plt_my$res <- as.numeric(plt_my$res)
plt_my <- plt_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (plt=res)%>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+")))%>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
plt_my <- within(plt_my, {
  id<-factor(id)
})
plt_my$sex<-factor(plt_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
plt_my$case <- 1

#MGUS CASES
plt_mgus <-filter(mgus_pat_wbloods, test_name=="Platelets")
plt_mgus$res<-plt_mgus$result
plt_mgus$res[plt_mgus$res == ""] <- NA
plt_mgus$res <- as.numeric(plt_mgus$res)
plt_mgus <- plt_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (plt=res)%>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+")))%>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
plt_mgus <- within(plt_mgus, {
  id<-factor(id)
})
plt_mgus$sex<-factor(plt_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
plt_mgus$case <-2 

#CONTROLS
plt_con <-filter(con_wbloods, test_name=="Platelets")
plt_con$res<-plt_con$result
plt_con$res[plt_con$res == ""] <- NA
plt_con$res <- as.numeric(plt_con$res)
plt_con <- plt_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (plt=res)%>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+")))%>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
plt_con <- within(plt_con, {
  id<-factor(id)
})
plt_con$sex<-factor(plt_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
plt_con$case <-3 

#append myeloma, MGUS and controls 
plt<-rbind(plt_my, plt_mgus, plt_con)
#case as factor variable
plt$case<-factor(plt$case, levels=c('1', '2', '3'), 
                 labels=c('Myeloma', 'MGUS', 'Control'))
table(plt$case)

#remove values >1000
plt %>% subset(plt>1000) %>% count()
plt %>% subset(plt>1000) %>% group_by(case) %>% count()
plt %>% subset(plt>1000) %>% group_by(case) %>% distinct(id) %>% count()

plt<- plt %>% subset(plt<=1000)

#save as R file
save(plt, file="N:\\DataAnalysis\\Trajectories\\Platelets\\Plt_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#PROTHROMBIN TIME

#MYELOMA CASES
pt_my <-filter(pat_wbloods, test_name=="Prothrombin time")
pt_my$res<-pt_my$result
pt_my$res[pt_my$res == ""] <- NA
pt_my$res <- as.numeric(pt_my$res)
pt_my <- pt_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (pt=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
pt_my <- within(pt_my, {
  id<-factor(id)
})
pt_my$sex<-factor(pt_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
pt_my$case <- 1 

#MGUS CASES
pt_mgus <-filter(mgus_pat_wbloods, test_name=="Prothrombin time")
pt_mgus$res<-pt_mgus$result
pt_mgus$res[pt_mgus$res == ""] <- NA
pt_mgus$res <- as.numeric(pt_mgus$res)
pt_mgus <- pt_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (pt=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
pt_mgus <- within(pt_mgus, {
  id<-factor(id)
})
pt_mgus$sex<-factor(pt_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
pt_mgus$case <-2 

#CONTROLS
pt_con <-filter(con_wbloods, test_name=="Prothrombin time")
pt_con$res<-pt_con$result
pt_con$res[pt_con$res == ""] <- NA
pt_con$res <- as.numeric(pt_con$res)
pt_con <- pt_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (pt=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
pt_con <- within(pt_con, {
  id<-factor(id)
})
pt_con$sex<-factor(pt_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
pt_con$case <-3 

#append myeloma, MGUS and controls 
pt<-rbind(pt_my, pt_mgus, pt_con)
#case as factor variable
pt$case<-factor(pt$case, levels=c('1', '2', '3'), 
                labels=c('Myeloma', 'MGUS', 'Control'))
table(pt$case)
save(pt, file="N:\\DataAnalysis\\Trajectories\\ProthrombinTime\\PT_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#PLASMA VISCOSITY

#MYELOMA CASES
pv_my <-filter(pat_wbloods, test_name=="Plasma viscosity")
pv_my$res<-pv_my$result
pv_my$res[pv_my$res == ""] <- NA
pv_my$res <- as.numeric(pv_my$res)
pv_my <- pv_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (pv=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
pv_my <- within(pv_my, {
  id<-factor(id)
})
pv_my$sex<-factor(pv_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
pv_my$case <- 1

#MGUS CASES
pv_mgus <-filter(mgus_pat_wbloods, test_name=="Plasma viscosity")
pv_mgus$res<-pv_mgus$result
pv_mgus$res[pv_mgus$res == ""] <- NA
pv_mgus$res <- as.numeric(pv_mgus$res)
pv_mgus <- pv_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (pv=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
pv_mgus <- within(pv_mgus, {
  id<-factor(id)
})
pv_mgus$sex<-factor(pv_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
pv_mgus$case <-2 

#CONTROLS 
pv_con <-filter(con_wbloods, test_name=="Plasma viscosity")
pv_con$res<-pv_con$result
pv_con$res[pv_con$res == ""] <- NA
pv_con$res <- as.numeric(pv_con$res)
pv_con <- pv_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (pv=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
pv_con <- within(pv_con, {
  id<-factor(id)
})
pv_con$sex<-factor(pv_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
pv_con$case <-3 

#append myeloma, MGUS and controls and save
pv<-rbind(pv_my, pv_mgus, pv_con)
#case as factor variable
pv$case<-factor(pv$case, levels=c('1', '2', '3'), 
                labels=c('Myeloma', 'MGUS', 'Control'))
table(pv$case)
save(pv, file="N:\\DataAnalysis\\Trajectories\\PV\\PV_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#RED CELL WIDTH 

#MYELOMA CASES
rcw_my <-filter(pat_wbloods, test_name=="RBC dist. width")
rcw_my$res<-rcw_my$result
rcw_my$res[rcw_my$res == ""] <- NA
rcw_my$res <- as.numeric(rcw_my$res)
rcw_my <- rcw_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (rcw=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
rcw_my <- within(rcw_my, {
  id<-factor(id)
})
rcw_my$sex<-factor(rcw_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
rcw_my$case <- 1

#MGUS CASES
rcw_mgus <-filter(mgus_pat_wbloods, test_name=="RBC dist. width")
rcw_mgus$res<-rcw_mgus$result
rcw_mgus$res[rcw_mgus$res == ""] <- NA
rcw_mgus$res <- as.numeric(rcw_mgus$res)
rcw_mgus <- rcw_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (rcw=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
rcw_mgus <- within(rcw_mgus, {
  id<-factor(id)
})
rcw_mgus$sex<-factor(rcw_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
rcw_mgus$case <-2 

#CONTROLS 
rcw_con <-filter(con_wbloods, test_name=="RBC dist. width")
rcw_con$res<-rcw_con$result
rcw_con$res[rcw_con$res == ""] <- NA
rcw_con$res <- as.numeric(rcw_con$res)
rcw_con <- rcw_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (rcw=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
rcw_con <- within(rcw_con, {
  id<-factor(id)
})
rcw_con$sex<-factor(rcw_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
rcw_con$case <-3 

#append myeloma, MGUS controls and save 
rcw<-rbind(rcw_my, rcw_mgus, rcw_con)
#case as factor variable
rcw$case<-factor(rcw$case, levels=c('1', '2', '3'), 
                 labels=c('Myeloma', 'MGUS', 'Control'))
table(rcw$case)
save(rcw, file="N:\\DataAnalysis\\Trajectories\\RBCDist\\RCW_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#RED CELL COUNT 

#MYELOMA CASES
rbc_my <-filter(pat_wbloods, test_name=="Red cell count RBC")
rbc_my$res<-rbc_my$result
rbc_my$res[rbc_my$res == ""] <- NA
rbc_my$res <- as.numeric(rbc_my$res)
rbc_my <- rbc_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (rbc=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+")))%>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
rbc_my <- within(rbc_my, {
  id<-factor(id)
})
rbc_my$sex<-factor(rbc_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
rbc_my$case <-1

#MGUS CASES
rbc_mgus <-filter(mgus_pat_wbloods, test_name=="Red cell count RBC")
rbc_mgus$res<-rbc_mgus$result
rbc_mgus$res[rbc_mgus$res == ""] <- NA
rbc_mgus$res <- as.numeric(rbc_mgus$res)
rbc_mgus <- rbc_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (rbc=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+")))%>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
rbc_mgus <- within(rbc_mgus, {
  id<-factor(id)
})
rbc_mgus$sex<-factor(rbc_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
rbc_mgus$case <-2

#CONTROLS 
rbc_con <-filter(con_wbloods, test_name=="Red cell count RBC")
rbc_con$res<-rbc_con$result
rbc_con$res[rbc_con$res == ""] <- NA
rbc_con$res <- as.numeric(rbc_con$res)
rbc_con <- rbc_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (rbc=res) %>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+")))%>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
rbc_con <- within(rbc_con, {
  id<-factor(id)
})
rbc_con$sex<-factor(rbc_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
rbc_con$case <-3

#append myeloma, MGUS and controls and save 
rbc<-rbind(rbc_my, rbc_mgus, rbc_con)
#case as factor variable
rbc$case<-factor(rbc$case, levels=c('1', '2', '3'), 
                 labels=c('Myeloma', 'MGUS', 'Control'))
table(rbc$case)
save(rbc, file="N:\\DataAnalysis\\Trajectories\\RedCellCount\\RBC_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#SERUM FREE LIGHT CHAINS 

#MYELOMA CASES
sflc_my<-filter(pat_wbloods, test_name=="sFLC Ratio (Siem)")
sflc_my$res<-sflc_my$result
sflc_my$res[sflc_my$res == ""] <- NA
sflc_my$res <- as.numeric(sflc_my$res)
sflc_my <- sflc_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (sflc=res) %>% 
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))

sflc_my <- within(sflc_my, {
  id<-factor(id)
})
sflc_my$sex<-factor(sflc_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
sflc_my$case<-1

#MGUS CASES
sflc_mgus<-filter(mgus_pat_wbloods, test_name=="sFLC Ratio (Siem)")
sflc_mgus$res<-sflc_mgus$result
sflc_mgus$res[sflc_mgus$res == ""] <- NA
sflc_mgus$res <- as.numeric(sflc_mgus$res)
sflc_mgus <- sflc_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (sflc=res) %>% 
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))

sflc_mgus <- within(sflc_mgus, {
  id<-factor(id)
})
sflc_mgus$sex<-factor(sflc_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
sflc_mgus$case<-2

#CONTROLS 
sflc_con<-filter(con_wbloods, test_name=="sFLC Ratio (Siem)")
sflc_con$res<-sflc_con$result
sflc_con$res[sflc_con$res == ""] <- NA
sflc_con$res <- as.numeric(sflc_con$res)
sflc_con <- sflc_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (sflc=res) %>% 
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))

sflc_con <- within(sflc_con, {
  id<-factor(id)
})
sflc_con$sex<-factor(sflc_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
sflc_con$case<-3

#append myeloma, MGUS, controls and save 
sflc<-rbind(sflc_my, sflc_mgus, sflc_con)
#case as factor variable
sflc$case<-factor(sflc$case, levels=c('1', '2', '3'), 
                  labels=c('Myeloma', 'MGUS', 'Control'))
table(sflc$case)
save(sflc, file="N:\\DataAnalysis\\Trajectories\\SFLC\\SFLC_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#UREA

#MYELOMA CASES
urea_my <-filter(pat_wbloods, test_name=="Urea")
urea_my$res<-urea_my$result
urea_my$res[urea_my$res == ""] <- NA
urea_my$res <- as.numeric(urea_my$res)
urea_my<- urea_my %>%
  select(id, res, days_prior, sex, age, dept) %>%
  na.omit() %>%
  rename (urea=res)%>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+")))%>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
urea_my <- within(urea_my, {
  id<-factor(id)
})
urea_my$sex<-factor(urea_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
urea_my$case <- 1

#MGUS CASES
urea_mgus <-filter(mgus_pat_wbloods, test_name=="Urea")
urea_mgus$res<-urea_mgus$result
urea_mgus$res[urea_mgus$res == ""] <- NA
urea_mgus$res <- as.numeric(urea_mgus$res)
urea_mgus <- urea_mgus %>%
  select(id, res, days_prior, sex, age, dept) %>%
  na.omit() %>%
  rename (urea=res)%>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+")))%>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
urea_mgus <- within(urea_mgus, {
  id<-factor(id)
})
urea_mgus$sex<-factor(urea_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
urea_mgus$case <- 2

#CONTROLS 
urea_con <-filter(con_wbloods, test_name=="Urea")
urea_con$res<-urea_con$result
urea_con$res[urea_con$res == ""] <- NA
urea_con$res <- as.numeric(urea_con$res)
urea_con <- urea_con %>%
  select(id, res, days_prior, sex, age, dept) %>%
  na.omit() %>%
  rename (urea=res)%>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+")))%>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
urea_con <- within(urea_con, {
  id<-factor(id)
})
urea_con$sex<-factor(urea_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
urea_con$case <- 3

#append myeloma, MGUS, controls and save 
urea<-rbind(urea_my, urea_mgus, urea_con)
#case as factor variable
urea$case<-factor(urea$case, levels=c('1', '2', '3'), 
                  labels=c('Myeloma', 'MGUS', 'Control'))
table(urea$case)

#remove values >50
urea %>% subset(urea>50) %>% count()
urea %>% subset(urea>50) %>% group_by(case) %>% count()
urea %>% subset(urea>50) %>% group_by(case) %>% distinct(id) %>% count()

urea<- urea %>% subset(urea<=50)

save(urea, file="N:\\DataAnalysis\\Trajectories\\Urea\\Urea_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# URINE CREATININE

#MYELOMA CASES
uc_my <-filter(pat_wbloods, test_name=="Urine Creatinine")
uc_my$flag1<-ifelse(grepl("<0.1", uc_my$result, ignore.case = T), 1, 0)
uc_my$flag2<-ifelse(grepl("<1.5", uc_my$result, ignore.case = T), 1, 0)
uc_my$flag<-ifelse(grepl("<|>", uc_my$result, ignore.case = T), 1, 0)
table(subset(uc_my, flag==1)$result)
uc_my$res<-uc_my$result
uc_my$res[uc_my$res == ""] <- NA
uc_my<-uc_my %>% drop_na(res)
uc_my$res <- sub("<", NA, uc_my$res)
uc_my$res <- as.numeric(uc_my$res)
#replace < and > values with actual value
uc_my <- uc_my %>% mutate(uc_imp=
                            case_when(
                              (flag==0) ~ res, 
                              (flag==1 & result=="<0.1") ~ 0.1,
                              (flag==1 & result=="<1.5") ~ 1.5))
uc_my$uc_imp <- as.numeric(uc_my$uc_imp)
uc_my<- uc_my %>%
  select(id, uc_imp, days_prior, age, sex, dept) %>%
  rename (uc=uc_imp,) %>% 
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
uc_my <- within(uc_my, {
  id<-factor(id)
})
uc_my$sex<-factor(uc_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
uc_my$case <- 1

#MGUS CASES
uc_mgus <-filter(mgus_pat_wbloods, test_name=="Urine Creatinine")
uc_mgus$flag1<-ifelse(grepl("<0.1", uc_mgus$result, ignore.case = T), 1, 0)
uc_mgus$flag2<-ifelse(grepl("<1.5", uc_mgus$result, ignore.case = T), 1, 0)
uc_mgus$flag<-ifelse(grepl("<|>", uc_mgus$result, ignore.case = T), 1, 0)
table(subset(uc_mgus, flag==1)$result)
uc_mgus$res<-uc_mgus$result
uc_mgus$res[uc_mgus$res == ""] <- NA
uc_mgus<-uc_mgus %>% drop_na(res)
uc_mgus$res <- sub("<", NA, uc_mgus$res)
uc_mgus$res <- as.numeric(uc_mgus$res)
#replace < and > values with actual value
uc_mgus <- uc_mgus %>% mutate(uc_imp=
                                case_when(
                                  (flag==0) ~ res, 
                                  (flag==1 & result=="<0.1") ~ 0.1,
                                  (flag==1 & result=="<1.5") ~ 1.5))
uc_mgus$uc_imp <- as.numeric(uc_mgus$uc_imp)
uc_mgus<- uc_mgus %>%
  select(id, uc_imp, days_prior, age, sex, dept) %>%
  rename (uc=uc_imp,) %>% 
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
uc_mgus <- within(uc_mgus, {
  id<-factor(id)
})
uc_mgus$sex<-factor(uc_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
uc_mgus$case <- 2

#CONTROLS 
uc_con <-filter(con_wbloods, test_name=="Urine Creatinine")
uc_con$flag1<-ifelse(grepl("<0.1", uc_con$result, ignore.case = T), 1, 0)
uc_con$flag2<-ifelse(grepl("<1.5", uc_con$result, ignore.case = T), 1, 0)
uc_con$flag<-ifelse(grepl("<|>", uc_con$result, ignore.case = T), 1, 0)
table(subset(uc_con, flag==1)$result)
uc_con$res<-uc_con$result
uc_con$res[uc_con$res == ""] <- NA
uc_con<-uc_con %>% drop_na(res)
uc_con$res <- sub("<", NA, uc_con$res)
uc_con$res <- as.numeric(uc_con$res)
#replace < and > values with actual value
uc_con <- uc_con %>% mutate(uc_imp=
                              case_when(
                                (flag==0) ~ res, 
                                (flag==1 & result=="<0.1") ~ 0.1,
                                (flag==1 & result=="<1.5") ~ 1.5))
uc_con$uc_imp <- as.numeric(uc_con$uc_imp)
uc_con<- uc_con %>%
  select(id, uc_imp, days_prior, age, sex, dept) %>%
  rename (uc=uc_imp,) %>% 
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
uc_con <- within(uc_con, {
  id<-factor(id)
})
uc_con$sex<-factor(uc_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
uc_con$case <- 3

#append myeloma and MGUS 
uc<-rbind(uc_my, uc_mgus, uc_con)
#case as factor variable
uc$case<-factor(uc$case, levels=c('1', '2', '3'), 
                labels=c('Myeloma', 'MGUS', 'Control'))
table(uc$case)
save(uc, file="N:\\DataAnalysis\\Trajectories\\UrineCreatinine\\UC_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Urine Prot/Cre Ratio 

#MYELOMA CASES
upcr_my <-filter(pat_wbloods, test_name=="Urine Prot/Cre Ratio")
upcr_my$res<-upcr_my$result
upcr_my$res[upcr_my$res == ""] <- NA
upcr_my$res <- as.numeric(upcr_my$res)
upcr_my <- upcr_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (upcr=res)%>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
upcr_my <- within(upcr_my, {
  id<-factor(id)
})
upcr_my$sex<-factor(upcr_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
upcr_my$case <- 1

#MGUS CASES
upcr_mgus <-filter(mgus_pat_wbloods, test_name=="Urine Prot/Cre Ratio")
upcr_mgus$res<-upcr_mgus$result
upcr_mgus$res[upcr_mgus$res == ""] <- NA
upcr_mgus$res <- as.numeric(upcr_mgus$res)
upcr_mgus <- upcr_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (upcr=res)%>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
upcr_mgus <- within(upcr_mgus, {
  id<-factor(id)
})
upcr_mgus$sex<-factor(upcr_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
upcr_mgus$case <- 2

#CONTROLS 
upcr_con <-filter(con_wbloods, test_name=="Urine Prot/Cre Ratio")
upcr_con$res<-upcr_con$result
upcr_con$res[upcr_con$res == ""] <- NA
upcr_con$res <- as.numeric(upcr_con$res)
upcr_con <- upcr_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (upcr=res)%>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
upcr_con <- within(upcr_con, {
  id<-factor(id)
})
upcr_con$sex<-factor(upcr_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
upcr_con$case <- 3

#append myeloma, MGUS and controls and save 
upcr<-rbind(upcr_my, upcr_mgus, upcr_con)
#case as factor variable
upcr$case<-factor(upcr$case, levels=c('1', '2', '3'),
                  labels=c('Myeloma', 'MGUS', 'Control'))
table(upcr$case)
save(upcr, file="N:\\DataAnalysis\\Trajectories\\UrineProtCreRatio\\UPCR_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#VITAMIN D

#MYELOMA CASES
vitd_my<-filter(pat_wbloods, test_name=="25 OH Vit D (Total)"| test_name=="Total 25OH Vitamin D")
table(vitd_my$units)
describe(vitd_my$result)
#some coded as <20.20.0 (n=14) replace with 20 (lowest recorded value is 11)
vitd_my$flag<-ifelse(grepl("<", vitd_my$result, ignore.case = T), 1, 0)
vitd_my %>% tally(flag)
table(subset(vitd_my, flag==1)$result)
vitd_my$res<-vitd_my$result
vitd_my$res[vitd_my$res == ""] <- NA
vitd_my<-vitd_my %>% drop_na(res)
vitd_my$res <- sub("<", NA, vitd_my$res)
vitd_my$res <- as.numeric(vitd_my$res)
summary(vitd_my$res)
#replace < and > values with actual value
vitd_my <- vitd_my %>% mutate(vitd_imp=
                                case_when(
                                  (flag==0) ~ res, 
                                  (flag==1) ~ 20))
vitd_my$vitd <- as.numeric(vitd_my$vitd_imp)
vitd_my<- vitd_my %>%
  select(id, vitd_imp, days_prior, age, sex, dept) %>%
  rename (vitd=vitd_imp,) %>% 
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
vitd_my <- within(vitd_my, {
  id<-factor(id)
})
vitd_my$sex<-factor(vitd_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
vitd_my$case <- 1

#MGUS
vitd_mgus<-filter(mgus_pat_wbloods, test_name=="25 OH Vit D (Total)"| test_name=="Total 25OH Vitamin D")
table(vitd_mgus$units)
describe(vitd_mgus$result)
#some coded as <20.20.0 (n=17) replace with 20 (lowest recorded value is 11)
vitd_mgus$flag<-ifelse(grepl("<", vitd_mgus$result, ignore.case = T), 1, 0)
vitd_mgus %>% tally(flag)
table(subset(vitd_mgus, flag==1)$result)
vitd_mgus$res<-vitd_mgus$result
vitd_mgus$res[vitd_mgus$res == ""] <- NA
vitd_mgus<-vitd_mgus %>% drop_na(res)
vitd_mgus$res <- sub("<", NA, vitd_mgus$res)
vitd_mgus$res <- as.numeric(vitd_mgus$res)
summary(vitd_mgus$res)
#replace < and > values with actual value
vitd_mgus <- vitd_mgus %>% mutate(vitd_imp=
                                    case_when(
                                      (flag==0) ~ res, 
                                      (flag==1) ~ 20))
vitd_mgus$vitd <- as.numeric(vitd_mgus$vitd_imp)
vitd_mgus<- vitd_mgus %>%
  select(id, vitd_imp, days_prior, age, sex, dept) %>%
  rename (vitd=vitd_imp,) %>% 
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
vitd_mgus <- within(vitd_mgus, {
  id<-factor(id)
})
vitd_mgus$sex<-factor(vitd_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
vitd_mgus$case <- 2

#Control
vitd_con<-filter(con_wbloods, test_name=="25 OH Vit D (Total)"| test_name=="Total 25OH Vitamin D")
table(vitd_con$units)
describe(vitd_con$result)
#some coded as <20.20.0 (n=25) replace with 20 (lowest recorded value is 11)
vitd_con$flag<-ifelse(grepl("<", vitd_con$result, ignore.case = T), 1, 0)
vitd_con %>% tally(flag)
table(subset(vitd_con, flag==1)$result)
vitd_con$res<-vitd_con$result
vitd_con$res[vitd_con$res == ""] <- NA
vitd_con<-vitd_con %>% drop_na(res)
vitd_con$res <- sub("<", NA, vitd_con$res)
vitd_con$res <- as.numeric(vitd_con$res)
summary(vitd_con$res)
#replace < and > values with actual value
vitd_con <- vitd_con %>% mutate(vitd_imp=
                                  case_when(
                                    (flag==0) ~ res, 
                                    (flag==1) ~ 20))
vitd_con$vitd <- as.numeric(vitd_con$vitd_imp)
vitd_con<- vitd_con %>%
  select(id, vitd_imp, days_prior, age, sex, dept) %>%
  rename (vitd=vitd_imp,) %>% 
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+"))) %>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
vitd_con <- within(vitd_con, {
  id<-factor(id)
})
vitd_con$sex<-factor(vitd_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
vitd_con$case <- 3

#append myeloma and MGUS 
vitd<-rbind(vitd_my, vitd_mgus, vitd_con)
#case as factor variable
vitd$case<-factor(vitd$case, levels=c('1', '2', '3'), 
                  labels=c('Myeloma', 'MGUS', 'Control'))
table(vitd$case)

#remove values >200
vitd %>% subset(vitd>200) %>% count()
vitd %>% subset(vitd>200) %>% group_by(case) %>% count()
vitd %>% subset(vitd>200) %>% group_by(case) %>% distinct(id) %>% count()
vitd<- vitd %>% subset(vitd<=200)

#save as R file 
save(vitd, file="N:\\DataAnalysis\\Trajectories\\VitaminD\\VitaminD_long.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#WHITE CELL COUNT

#MYELOMA CASES
wcc_my <-filter(pat_wbloods, test_name=="White cell count")
wcc_my$res<-wcc_my$result
wcc_my$res[wcc_my$res == ""] <- NA
wcc_my$res <- as.numeric(wcc_my$res)
wcc_my <- wcc_my %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (wcc=res)%>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+")))%>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
wcc_my <- within(wcc_my, {
  id<-factor(id)
})
wcc_my$sex<-factor(wcc_my$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
wcc_my$case <- 1 

#MGUS CASES
wcc_mgus <-filter(mgus_pat_wbloods, test_name=="White cell count")
wcc_mgus$res<-wcc_mgus$result
wcc_mgus$res[wcc_mgus$res == ""] <- NA
wcc_mgus$res <- as.numeric(wcc_mgus$res)
wcc_mgus <- wcc_mgus %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (wcc=res)%>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+")))%>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
wcc_mgus <- within(wcc_mgus, {
  id<-factor(id)
})
wcc_mgus$sex<-factor(wcc_mgus$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
wcc_mgus$case <- 2

#CONTROLS 
wcc_con <-filter(con_wbloods, test_name=="White cell count")
wcc_con$res<-wcc_con$result
wcc_con$res[wcc_con$res == ""] <- NA
wcc_con$res <- as.numeric(wcc_con$res)
wcc_con <- wcc_con %>%
  select(id, res, days_prior, age, sex, dept) %>%
  na.omit() %>%
  rename (wcc=res)%>%
  mutate(age_cat=cut(age, breaks=c(30, 50, 60, 70, 80, 100), 
                     right=FALSE, 
                     label=c("<50", "50-59", "60-69", "70-79", "80+")))%>%
  filter(days_prior>-1800) %>%
  mutate(month_prior=cut(days_prior, breaks=seq(0, -1800, by=-30),
                         right=FALSE,label=-60:-1))
wcc_con <- within(wcc_con, {
  id<-factor(id)
})
wcc_con$sex<-factor(wcc_con$sex, levels=c('1', '2'), labels=c('Male', 'Female'))
wcc_con$case <- 3

#append myeloma, MGUS and controls and save 
wcc<-rbind(wcc_my, wcc_mgus, wcc_con)
#case as factor variable
wcc$case<-factor(wcc$case, levels=c('1', '2', '3'),
                 labels=c('Myeloma', 'MGUS', 'Control'))
table(wcc$case)

#remove values >50
wcc %>% subset(wcc>50) %>% count()
wcc %>% subset(wcc>50) %>% group_by(case) %>% count()
wcc %>% subset(wcc>50) %>% group_by(case) %>% distinct(id) %>% count()

wcc<- wcc %>% subset(wcc<=50)

save(wcc, file="N:\\DataAnalysis\\Trajectories\\Wcc\\Wcc_long.RData")

#########################################################################
#load all cleaned blood test data to append all 
load("N:\\DataAnalysis\\Trajectories\\Albumin\\Alb_long.RData")
tc1<- alb %>% select(id, alb, case, days_prior, age, sex, dept) %>% 
  rename(res=alb) %>% mutate(test="Albumin")
rm(alb)

load("N:\\DataAnalysis\\Trajectories\\ALK-Phos\\ALKP_long.RData")
tc2<- alkp %>% select(id, alk, case, days_prior, age, sex, dept) %>% 
  rename(res=alk) %>% mutate(test="ALk-Phos")
rm(alkp)

load("N:\\DataAnalysis\\Trajectories\\ALT\\Alt_long.RData")
tc3<- alt %>% select(id, alt, case, days_prior, age, sex, dept) %>% 
  rename(res=alt) %>% mutate(test="ALT")
rm(alt)

load("N:\\DataAnalysis\\Trajectories\\APTT\\APTT_long.RData")
tc4<- aptt %>% select(id, aptt, case, days_prior, age, sex, dept) %>% 
  rename(res=aptt) %>% mutate(test="APTT")
rm(aptt)

load("N:\\DataAnalysis\\Trajectories\\Basophil\\Bas_long.RData")
tc5<- bas %>% select(id, bas, case, days_prior, age, sex, dept) %>% 
  rename(res=bas) %>% mutate(test="Basophil")
rm(bas)

load("N:\\DataAnalysis\\Trajectories\\Beta2Microglobulin\\B2M_long.RData")
tc6<- b2m %>% select(id, b2m, case, days_prior, age, sex, dept) %>% 
  rename(res=b2m) %>% mutate(test="Beta2microglobulin")
rm(b2m)

load("N:\\DataAnalysis\\Trajectories\\Calcium\\Cal_long.RData")
tc7<- cal %>% select(id, cal, case, days_prior, age, sex, dept) %>% 
  rename(res=cal) %>% mutate(test="Calcium")
rm(cal)

load("N:\\DataAnalysis\\Trajectories\\Creatinine\\Creat_long.RData")
tc8<- creat %>% select(id, creat, case, days_prior, age, sex, dept) %>% 
  rename(res=creat) %>% mutate(test="Creatinine")
rm(creat)

load("N:\\DataAnalysis\\Trajectories\\CRP\\CRP_long.RData")
tc9<- crp %>% select(id, crp, case, days_prior, age, sex, dept) %>% 
  rename(res=crp) %>% mutate(test="CRP")
rm(crp)

load("N:\\DataAnalysis\\Trajectories\\DDimer\\DDimer_long.RData")
tc10<- dd %>% select(id, dd, case, days_prior, age, sex, dept) %>% 
  rename(res=dd) %>% mutate(test="DDimer")
rm(dd)

load("N:\\DataAnalysis\\Trajectories\\eGFR\\eGFR_long.RData")
tc11<- egfr %>% select(id, egfr_cat, case, days_prior, age, sex, dept) %>% 
  rename(res=egfr_cat) %>% mutate(test="eGFR")
rm(egfr)

load("N:\\DataAnalysis\\Trajectories\\Eosinophil\\Eos_long.RData")
tc12<- eos %>% select(id, eos, case, days_prior, age, sex, dept) %>% 
  rename(res=eos) %>% mutate(test="Eosinophil")
rm(eos)

load("N:\\DataAnalysis\\Trajectories\\ESR\\ESR_long.RData")
tc13 <- esr %>% select(id, esr, case, days_prior, age, sex, dept) %>% 
  rename(res=esr) %>% mutate(test="ESR")
rm(esr)

load("N:\\DataAnalysis\\Trajectories\\Fibrinogen\\DF_long.RData")
tc14 <- df %>% select(id, df_cat, case, days_prior, age, sex, dept) %>% 
  rename(res=df_cat) %>% mutate(test="Derived Fibrinogen")
rm(df)

load("N:\\DataAnalysis\\Trajectories\\Haemoglobin\\Haem_long.RData")
tc15<- haem %>% select(id, haem, case, days_prior, age, sex, dept) %>% 
  rename(res=haem) %>% mutate(test="Haemoglobin")
rm(haem)

load("N:\\DataAnalysis\\Trajectories\\Hypochromic\\HC_long.RData")
tc16<- hc %>% select(id, hc, case, days_prior, age, sex, dept) %>% 
  rename(res=hc) %>% mutate(test="Hypochromic")
rm(hc)

load("N:\\DataAnalysis\\Trajectories\\INR\\INR_long.RData")
tc17<- inr %>% select(id, inr, case, days_prior, age, sex, dept) %>% 
  rename(res=inr) %>% mutate(test="INR")
rm(inr)

load("N:\\DataAnalysis\\Trajectories\\KappaLC\\FKLC_long.RData")
tc18<- fklc %>% select(id, fklc, case, days_prior, age, sex, dept) %>% 
  rename(res=fklc) %>% mutate(test="Free Kappa LC")
rm(fklc)

load("N:\\DataAnalysis\\Trajectories\\LambdaLC\\FLLC_long.RData")
tc19<- fllc %>% select(id, fllc, case, days_prior, age, sex, dept) %>% 
  rename(res=fllc) %>% mutate(test="Free Lambda LC")
rm(fllc)

load("N:\\DataAnalysis\\Trajectories\\LDH\\LDH_long.RData")
tc20<- ldh %>% select(id, ldh, case, days_prior, age, sex, dept) %>% 
  rename(res=ldh) %>% mutate(test="LDH")
rm(ldh)

load("N:\\DataAnalysis\\Trajectories\\Lymphocyte\\Lymph_long.RData")
tc21<- lymph %>% select(id, lymph, case, days_prior, age, sex, dept) %>% 
  rename(res=lymph) %>% mutate(test="Lymphocyte")
rm(lymph)

load("N:\\DataAnalysis\\Trajectories\\MCV\\MCV_long.RData")
tc22<- mcv %>% select(id, mcv, case, days_prior, age, sex, dept) %>% 
  rename(res=mcv) %>% mutate(test="Mean Cell Volume")
rm(mcv)

load("N:\\DataAnalysis\\Trajectories\\MeanCorpusc\\MCH_long.RData")
tc23<- mch %>% select(id, mch, case, days_prior, age, sex, dept) %>% 
  rename(res=mch) %>% mutate(test="Mean Corpusc")
rm(mch)

load("N:\\DataAnalysis\\Trajectories\\Monocyte\\Mono_long.RData")
tc24<- mono %>% select(id, mono, case, days_prior, age, sex, dept) %>% 
  rename(res=mono) %>% mutate(test="Monocyte")
rm(mono)

load("N:\\DataAnalysis\\Trajectories\\Neutrophils\\Neut_long.RData")
tc25<- neut %>% select(id, neut, case, days_prior, age, sex, dept) %>% 
  rename(res=neut) %>% mutate(test="Neutrophil")
rm(neut)

load("N:\\DataAnalysis\\Trajectories\\PackedCellVolume\\PCV_long.RData")
tc26<- pcv %>% select(id, pcv, case, days_prior, age, sex, dept) %>% 
  rename(res=pcv) %>% mutate(test="Packed Cell Volume")
rm(pcv)

load("N:\\DataAnalysis\\Trajectories\\Platelets\\Plt_long.RData")
tc27<- plt %>% select(id, plt, case, days_prior, age, sex, dept) %>% 
  rename(res=plt) %>% mutate(test="Platelets")
rm(plt)

load("N:\\DataAnalysis\\Trajectories\\ProthrombinTime\\PT_long.RData")
tc28<- pt %>% select(id, pt, case, days_prior, age, sex, dept) %>% 
  rename(res=pt) %>% mutate(test="Prothrombin Time")
rm(pt)

load("N:\\DataAnalysis\\Trajectories\\PV\\PV_long.RData")
tc29<- pv %>% select(id, pv, case, days_prior, age, sex, dept) %>% 
  rename(res=pv) %>% mutate(test="PV")
rm(pv)

load("N:\\DataAnalysis\\Trajectories\\RBCDist\\RCW_long.RData")
tc30<- rcw %>% select(id, rcw, case, days_prior, age, sex, dept) %>% 
  rename(res=rcw) %>% mutate(test="RBC Dist")
rm(rcw)

load("N:\\DataAnalysis\\Trajectories\\RedCellCount\\RBC_long.RData")
tc31<- rbc %>% select(id, rbc, case, days_prior, age, sex, dept) %>% 
  rename(res=rbc) %>% mutate(test="Red Cell Count")
rm(rbc)

load("N:\\DataAnalysis\\Trajectories\\SFLC\\SFLC_long.RData")
tc32<- sflc %>% select(id, sflc, case, days_prior, age, sex, dept) %>% 
  rename(res=sflc) %>% mutate(test="Serum Free LC")
rm(sflc)

load("N:\\DataAnalysis\\Trajectories\\Urea\\Urea_long.RData")
tc33<- urea %>% select(id, urea, case, days_prior, age, sex, dept) %>% 
  rename(res=urea) %>% mutate(test="Urea")
rm(urea)

load("N:\\DataAnalysis\\Trajectories\\UrineCreatinine\\UC_long.RData")
tc34<- uc %>% select(id, uc, case, days_prior, age, sex, dept) %>% 
  rename(res=uc) %>% mutate(test="Urine Creatinine")
rm(uc)

load("N:\\DataAnalysis\\Trajectories\\UrineProtCreRatio\\UPCR_long.RData")
tc35<- upcr %>% select(id, upcr, case, days_prior, age, sex, dept) %>% 
  rename(res=upcr) %>% mutate(test="Urine Prot Creat Ratio")
rm(upcr)

load("N:\\DataAnalysis\\Trajectories\\VitaminD\\VitaminD_long.RData")
tc36<- vitd %>% select(id, vitd, case, days_prior, age, sex, dept) %>% 
  rename(res=vitd) %>% mutate(test="Vit D")
rm(vitd)

load("N:\\DataAnalysis\\Trajectories\\WCC\\Wcc_long.RData")
tc37<- wcc %>% select(id, wcc, case, days_prior, age, sex, dept) %>% 
  rename(res=wcc) %>% mutate(test="WCC")
rm(wcc)

#append all data 
tc<-rbind(tc1, tc2, tc3, tc4, tc5, tc6, tc7, tc8, tc9, tc10,
          tc11, tc12, tc13, tc14, tc15, tc16, tc17, tc18, tc19, tc20,
          tc21, tc22, tc23, tc24, tc25, tc26, tc27, tc28, tc29, tc30,
          tc31, tc32, tc33, tc34, tc35, tc36, tc37)

#save 
save(tc, file="N:\\DataAnalysis\\FormattedData\\AllTestsCombined.RData")

tc %>% group_by(case) %>% summarise(n_distinct(id))
tc %>% group_by(case) %>% summarise(n())
table(tc$test)
table(tc$dept)

