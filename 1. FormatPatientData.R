########################################
#Development and internal validation of a risk prediction model to identify 
#myeloma based on routine blood tests: a case-control study


#Format myeloma, MGUS and control patient level data including comorbidities 


library(tidyverse)
library(Hmisc)

#clear R environment
rm(list=ls())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#MYELOMA CASES

#read in patient level data and rename columns
pat<- read.delim("N:\\ImportedData\\Incoming_2021-07-15\\LTH20102_Myeloma_data\\Myeloma_data\\LTH20102_Myeloma_Patients_Demographics.txt", header=TRUE, sep="|")
#rename columns
colnames(pat)<- c("id", "sex", "ethnicity", "diag_month", "diag_year", "age", "diagcode", 
                  "performance", "performance_lab", "imd2015", "imd2019")
summarise(pat, n_distinct(id)) 

#recode missing ethnicity
pat$ethnicity[pat$ethnicity=="Missing"]<-NA

#read in myeloma bloods and rename columns 
bloods<- read.delim("N:\\ImportedData\\Incoming_2021-07-15\\LTH20102_Myeloma_data\\Myeloma_data\\LTH20102_MyelomaPats_BloodResults.txt", header=TRUE, sep="|")
colnames(bloods)<- c("id", "t_code", "test_name", "result", "units", "ref_range", "days", "dept")
bloods$days_prior<-bloods$days*-1
summarise(bloods, n_distinct(id)) 

#summarise department
dept <- bloods %>%
  group_by(dept) %>%
  summarise(n_dept=n()) %>%
  arrange(desc(n_dept), by.group=TRUE)

#check no ids with ICU tests
icuids <- bloods %>% filter(dept=="Intensive Care") %>% select(id) %>% distinct(id) 
#remove blood tests in ICU
bloods <-bloods %>%  filter(dept!="Intensive Care")

#merge pat and bloods for those with bloods only 
pat_wbloods<-left_join(bloods, pat)
summarise(pat_wbloods, n_distinct(id)) 

#Keep distinct Ids of those with bloods only
ids_bl<-bloods %>%
  distinct(id, .keep_all=TRUE) %>%
  select(id)

#only keep patient file for those with bloods
pat <- ids_bl %>% 
  left_join(pat)

#save myeloma data
save(pat, pat_wbloods, file="N:\\DataAnalysis\\FormattedData\\Myeloma.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#MGUS patients
#read in patient data and rename columns
mgus_pat<- read.delim("N:\\ImportedData\\Incoming_2021-07-15\\LTH20102_MGUS_data\\MGUS_data\\LTH20102_MGUS_Patients_Demographics.txt",  header=TRUE, sep="|")
colnames(mgus_pat)<- c("id", "sex", "ethnicity", "diag_month", "diag_year", "age", "diagcode", 
                       "imd2015", "imd2019", "myeloma_flag")
summarise(mgus_pat, n_distinct(id)) 
#read in MGUS bloods and rename columns
mgus_bloods<- read_tsv("N:\\ImportedData\\Incoming_2021-07-15\\LTH20102_MGUS_data\\MGUS_data\\LTH20102_MGUSPats_BloodResults.txt",col_names=TRUE)
colnames(mgus_bloods)<- c("id", "t_code", "test_name", "result", "units", "ref_range", "days", "dept")
#recode days prior as negative 
mgus_bloods$days_prior<-mgus_bloods$days*-1
summarise(mgus_bloods, n_distinct(id)) 

#summarise department
dept <- mgus_bloods %>%
  group_by(dept) %>%
  summarise(n_dept=n()) %>%
  arrange(desc(n_dept), by.group=TRUE)

#check no ids with ICU tests
icuids <- mgus_bloods %>% filter(dept=="Intensive Care") %>% select(id) %>% distinct(id) 
#remove blood tests in ICU
mgus_bloods <-mgus_bloods %>%  filter(dept!="Intensive Care")
summarise(mgus_bloods, n_distinct(id)) 
summarise(mgus_bloods, n()) 

#merge pat and bloods for those with bloods only 
mgus_pat_wbloods<-left_join(mgus_bloods, mgus_pat)
summarise(mgus_pat_wbloods, n_distinct(id)) 

#Keep distinct Ids of those with bloods only
mgus_ids_bl<-mgus_bloods %>%
  distinct(id, .keep_all=TRUE) %>%
  select(id)

#only keep patient file for those with bloods
mgus_pat <- mgus_ids_bl %>% 
  left_join(mgus_pat)

#save MGUS
save(mgus_pat, mgus_pat_wbloods, file="N:\\DataAnalysis\\FormattedData\\MGUS.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Control data
#exclude those with previous cancer
#exclude those with tests from ICU only and lymphoctye >10
#exclude tests in ICU and Palliative care

#read in patient level data and rename columns
c_pat<- read.delim("N:\\ImportedData\\Incoming_2021-11-24\\LTH20102_Control_Demographics.txt", header=TRUE, sep="|")
#rename columns
colnames(c_pat)<- c("id", "sex", "ethnicity", "diag_month", "diag_year", "age", "imd2015", "imd2019")
summarise(c_pat, n_distinct(id)) 

#recode missing ethnicity
c_pat$ethnicity[c_pat$ethnicity=="Missing"]<-NA

#read in control bloods and rename columns
c_bloods<- read.delim("N:\\ImportedData\\Incoming_2021-11-24\\LTH20102_Control_BloodResults.txt", header=TRUE, sep="|")
colnames(c_bloods)<- c("id", "t_code", "test_name", "result", "units", "ref_range", "days", "dept")
c_bloods$days_prior<-c_bloods$days*-1
summarise(c_bloods, n_distinct(id)) 

#updated on 8/2/22 to exclude those with previous cancer based on comorbidity data
cm<- read.delim("N:\\ImportedData\\Incoming_2021-11-24\\LTH20102_Control_Prev_Diagnosis.txt", header=TRUE, sep="|")

#idenitfy thise with previous cancer
prevcan <- cm %>% 
  mutate(prevc= ifelse(
    grepl("C",cm$Diag_Code, ignore.case = T), 1, 0)) %>% 
  filter(prevc==1) %>% 
  select(PseudoStudy_ID) %>% distinct(PseudoStudy_ID)
omit<-as.vector(prevcan$PseudoStudy_ID)

#remove from both c_pat and c_bloods
c_pat <-c_pat %>%  filter(!(id %in% omit))
c_bloods <-c_bloods %>%  filter(!(id %in% omit))

#remove controls who only had blood tests in ICU 
icuids <- c_bloods %>% filter(dept=="Intensive Care") %>% select(id) %>% distinct(id) 
nonicuids <- c_bloods %>% filter(dept!="Intensive Care") %>% select(id) %>% distinct(id)
omiticu <-setdiff(icuids, nonicuids) %>% pull(id)

#remove from both c_pat and c_bloods
c_pat <-c_pat %>%  filter(!(id %in% omiticu))
c_bloods <-c_bloods %>%  filter(!(id %in% omiticu))

#remove controls who only had blood tests in palliarive care 
pcids <- c_bloods %>% filter(dept=="Palliative Care")%>% select(id) %>% distinct(id) 
nonpcids <- c_bloods %>% filter(dept!="Palliative Care") %>% select(id) %>% distinct(id)
omitpc <-setdiff(pcids, nonpcids) %>% pull(id)

#remove from both c_pat and c_bloods
c_pat <-c_pat %>%  filter(!(id %in% omitpc))
c_bloods <-c_bloods %>%  filter(!(id %in% omitpc))

#remove ids with lymphocyte value >10 
omitlymph<- c_bloods %>% filter(test_name=="Lymphocyte count") %>% 
  mutate(res=as.numeric(result)) %>% subset(res>10) %>% 
  select(id) %>% distinct(id) %>% pull(id)

c_pat <-c_pat %>%  filter(!(id %in% omitlymph))
c_bloods <-c_bloods %>%  filter(!(id %in% omitlymph))

###remove bloods from ICU
c_bloods <-c_bloods %>%  filter(dept!="Intensive Care")
summarise(c_bloods, n_distinct(id)) 
summarise(c_bloods, n()) 
###remove vlodds from Pallative care
c_bloods <-c_bloods %>%  filter(dept!="Palliative Care")
summarise(c_bloods, n_distinct(id)) 
summarise(c_bloods, n()) 

#merge pat and bloods for those with bloods only 
con_wbloods<-left_join(c_bloods, c_pat)
summarise(con_wbloods, n_distinct(id)) 

##rename dataframes to save
con<- c_pat

#save Control data
save(con, con_wbloods, file="N:\\DataAnalysis\\FormattedData\\Controls.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Format co-morbidity data for cases, MGUS and controls

#based on categories in Fowler et at, BMC Cancer 2020
#https://bmccancer.biomedcentral.com/articles/10.1186/s12885-019-6472-9

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Myeloma comorbidity data 
cm<- read.delim("N:\\ImportedData\\Incoming_2021-07-15\\LTH20102_Myeloma_data\\Myeloma_data\\LTH20102_Myeloma_Prev_Diagnosis.txt", header=TRUE, sep="|")
names(cm)

#create indicator for each co morbidity
#Diabetes
cm$diab<- ifelse(
  grepl("E10.0|E10.1|E10.2|E10.3|E10.4|E10.5|E10.6|E10.7|E10.8|E10.9|
              E11.0|E11.1|E11.2|E11.3|E11.4|E11.5|E11.6|E11.7|E11.8|E11.9|
              E12.0|E12.1|E12.2|E12.3|E12.4|E12.5|E12.6|E12.7|E12.8|E12.9|
              E13.0|E13.1|E13.2|E13.3|E13.4|E13.5|E13.6|E13.7|E13.8|E13.9|
              E14.0|E14.1|E14.2|E14.3|E14.4|E14.5|E14.6|E14.7|E14.8|E14.9",
        cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(diab)

#Liver disease
cm$liver<- ifelse(
  grepl("B18.0|B18.1|B18.2|B18.8|B18.9|
              I85.0|I85.9|I86.4|I98.2|
              K70|K70.0|K70.1|K70.2|K70.3|K70.4|K70.9|
              K71|K71.1|K71.3|K71.4|K71.5|K71.7|
              K72|K72.1|K72.9|K73|K74|K76|
              K76.0|K76.2|K76.3|K76.4|K76.5|K76.6|K76.7|K76.8|K76.9|Z94.4",
        cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(liver)

#Previous Cancer 
cm$can <- ifelse(
  grepl("C",cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(can)

#Metastatic cancer
cm$metscan <- ifelse(
  grepl("C77|C78|C79|C80",cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(metscan)

#Obesity
cm$obese<-ifelse(grepl("E66", cm$Diag_Code, ignore.case=T), 1, 0)
cm %>% tally(obese)

#Dementia
cm$dem <- ifelse(
  grepl("F00|F01|F02|F03|F05.1|G30|G31.1"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(dem)

#Hemiplegia or paraplegia
cm$para <- ifelse(
  grepl("G04.1|G11.4|G80.1|G80.2|G81|G82|G83.0|G83.1|G83.2|G83.4|G83.9"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(para)

#Stroke
cm$stroke<- ifelse(
  grepl("G45|G46|H34.0|I60|I61|I62|I63|I64|I65|I66|I67|I68|I69"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(stroke)

#hypertension
cm$hyp<- ifelse(
  grepl("I10|I11.9|I12.9|I13.9"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(hyp)

#renal disease
cm$renal<- ifelse(
  grepl("I12.0|I13.1|N03.2|N03.4|N03.5|N03.6|N03.7|
        N05.2|N05.3|N05.4|N05.5|N05.6|N05.7|N18|N19|N25.0|
        Z49.0|Z49.1|Z49.2|Z94.0|Z99.2"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(renal)

#Myocardial infarction
cm$mi<- ifelse(
  grepl("I21|I22|I25.2"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(mi)

#COPD
cm$copd<- ifelse(
  grepl("I27.8|I27.9|
        J40|J41|J42|J43|J4|J45|J46|J47|
        J60|J61|J62|J63|J64|J65|J66|J67|
        J68.4|J70.1|J70.3"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(copd)

#Congestive heart failure
cm$chf<- ifelse(
  grepl("I43|I50|I09.0|I11.0|I13.0|I13.2|I25.5|
        I42.0|I42.5|I42.6|I42.7|I42.8|I42.9|P29.0"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(chf)

#peripheral vascular disease
cm$pvd<- ifelse(
  grepl("I70|I71|I73.1|I73.8|I73.9|I77.1|I79.0|
        K55.1|K55.8|K55.9|Z95.8|Z95.9"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(pvd)

#Rheumatological conditioins
cm$rh<- ifelse(
  grepl("M05|M06|M31.5|M32|M33|M34|M35.1|M35.3|M36.0"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(rh)

#aggregate by study Id to get totals for each id
cm<- cm %>%
  group_by(PseudoStudy_ID) %>%
  summarise(c1=sum(diab), c2=sum(liver), c3=sum(can), c4=sum(metscan), 
            c5=sum(obese), c6=sum(dem), c7=sum(para), c8=sum(stroke), 
            c9=sum(hyp), c10=sum(renal), c11=sum(mi), c12=sum(copd), 
            c13=sum(chf), c14=sum(pvd), c15=sum(rh)) 
#recode any vales>1 as 1 to get binary indicator for each ID
cm$c1[cm$c1>1] <-1
cm$c2[cm$c2>1] <-1
cm$c3[cm$c3>1] <-1
cm$c4[cm$c4>1] <-1
cm$c5[cm$c5>1] <-1
cm$c6[cm$c6>1] <-1
cm$c7[cm$c7>1] <-1
cm$c8[cm$c8>1] <-1
cm$c9[cm$c9>1] <-1
cm$c10[cm$c10>1] <-1
cm$c11[cm$c11>1] <-1
cm$c12[cm$c12>1] <-1
cm$c13[cm$c13>1] <-1
cm$c14[cm$c14>1] <-1
cm$c15[cm$c15>1] <-1
#calcualte total comorbidities
cm$totcm<-cm$c1 + cm$c2  + cm$c3  + cm$c4  + cm$c5  + cm$c6  + cm$c7  + cm$c8 
+ cm$c9 + cm$c10 + cm$c11 + cm$c12 + cm$c13 + cm$c14 + cm$c15
head(cm)
summary(cm$totcm)
table(cm$totcm)

#rename columns in comrobidity data frame
colnames(cm)<- c("id", "diab", "liver", "can", "metscan", "obese", "dem", "para", 
                 "stroke", "hyp", "renal", "mi", "copd", "chf", "pvd", "rh", "totcm")

cm_my<-cm
save(cm_my, file="N:\\DataAnalysis\\FormattedData\\ComorbidMyeloma.Rdata")

#merge with patients with bloods - load myeloma saved files
load("N:\\DataAnalysis\\FormattedData\\Myeloma.RData")

#merge with patient data
pat_cm<-left_join(pat, cm)
summarise(pat_cm, n_distinct(id)) 

pat_cm <- pat_cm %>%
  mutate_at(vars(diab, liver, can, metscan, obese, dem, para, 
                 stroke, hyp, renal, mi, copd, chf, pvd, rh, totcm), 
            ~ replace_na(., 0))
pat_cm %>% tally(diab)
#table of total no. comorbidities
cmtab<-table(pat_cm$totcm)
t2<-prop.table(cmtab)
cbind(cmtab, t2)

#select comorbidity data only to create sums and %s
c1<- pat_cm %>% 
  select("diab", "liver", "can", "metscan", "obese", "dem", "para", 
         "stroke", "hyp", "renal", "mi", "copd", "chf", "pvd", "rh")
#total with each condition
apply(c1, 2, sum)
#% with each condition
apply(c1, 2, mean)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MGUS cohort 
rm(list=ls())
#comorbidity data for MGUS 
cm<- read.delim("N:\\ImportedData\\Incoming_2021-07-15\\LTH20102_MGUS_data\\MGUS_data\\LTH20102_MGUS_Prev_Diagnosis.txt", header=TRUE, sep="|")
names(cm)

#create indicator for each comorbidity
#Diabetes
cm$diab<- ifelse(
  grepl("E10.0|E10.1|E10.2|E10.3|E10.4|E10.5|E10.6|E10.7|E10.8|E10.9|
              E11.0|E11.1|E11.2|E11.3|E11.4|E11.5|E11.6|E11.7|E11.8|E11.9|
              E12.0|E12.1|E12.2|E12.3|E12.4|E12.5|E12.6|E12.7|E12.8|E12.9|
              E13.0|E13.1|E13.2|E13.3|E13.4|E13.5|E13.6|E13.7|E13.8|E13.9|
              E14.0|E14.1|E14.2|E14.3|E14.4|E14.5|E14.6|E14.7|E14.8|E14.9",
        cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(diab)

#Liver disease
cm$liver<- ifelse(
  grepl("B18.0|B18.1|B18.2|B18.8|B18.9|
              I85.0|I85.9|I86.4|I98.2|
              K70|K70.0|K70.1|K70.2|K70.3|K70.4|K70.9|
              K71|K71.1|K71.3|K71.4|K71.5|K71.7|
              K72|K72.1|K72.9|K73|K74|K76|
              K76.0|K76.2|K76.3|K76.4|K76.5|K76.6|K76.7|K76.8|K76.9|Z94.4",
        cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(liver)

#Previous Cancer
cm$can <- ifelse(
  grepl("C",cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(can)

#metastatic cancer
cm$metscan <- ifelse(
  grepl("C77|C78|C79|C80",cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(metscan)

#Obesity
cm$obese<-ifelse(grepl("E66", cm$Diag_Code, ignore.case=T), 1, 0)
cm %>% tally(obese)

#Dementia
cm$dem <- ifelse(
  grepl("F00|F01|F02|F03|F05.1|G30|G31.1"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(dem)

#Hemiplegia or paraplegia
cm$para <- ifelse(
  grepl("G04.1|G11.4|G80.1|G80.2|G81|G82|G83.0|G83.1|G83.2|G83.4|G83.9"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(para)

#Stroke
cm$stroke<- ifelse(
  grepl("G45|G46|H34.0|I60|I61|I62|I63|I64|I65|I66|I67|I68|I69"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(stroke)

#hypertension
cm$hyp<- ifelse(
  grepl("I10|I11.9|I12.9|I13.9"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(hyp)

#renal disease
cm$renal<- ifelse(
  grepl("I12.0|I13.1|N03.2|N03.4|N03.5|N03.6|N03.7|
        N05.2|N05.3|N05.4|N05.5|N05.6|N05.7|N18|N19|N25.0|
        Z49.0|Z49.1|Z49.2|Z94.0|Z99.2"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(renal)

#Myocardial infarction
cm$mi<- ifelse(
  grepl("I21|I22|I25.2"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(mi)

#COPD
cm$copd<- ifelse(
  grepl("I27.8|I27.9|
        J40|J41|J42|J43|J4|J45|J46|J47|
        J60|J61|J62|J63|J64|J65|J66|J67|
        J68.4|J70.1|J70.3"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(copd)

#Congestive heart failure
cm$chf<- ifelse(
  grepl("I43|I50|I09.0|I11.0|I13.0|I13.2|I25.5|
        I42.0|I42.5|I42.6|I42.7|I42.8|I42.9|P29.0"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(chf)

#peripheral vascular disease
cm$pvd<- ifelse(
  grepl("I70|I71|I73.1|I73.8|I73.9|I77.1|I79.0|
        K55.1|K55.8|K55.9|Z95.8|Z95.9"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(pvd)

#Rheumatological conditioins
cm$rh<- ifelse(
  grepl("M05|M06|M31.5|M32|M33|M34|M35.1|M35.3|M36.0"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(rh)

#aggregate by study Id to get totals for each id
cm<- cm %>%
  group_by(PseudoStudy_ID) %>%
  summarise(c1=sum(diab), c2=sum(liver), c3=sum(can), c4=sum(metscan), 
            c5=sum(obese), c6=sum(dem), c7=sum(para), c8=sum(stroke), 
            c9=sum(hyp), c10=sum(renal), c11=sum(mi), c12=sum(copd), 
            c13=sum(chf), c14=sum(pvd), c15=sum(rh)) 
#recode any vales>1 as 1 to get binary indicator for each ID
cm$c1[cm$c1>1] <-1
cm$c2[cm$c2>1] <-1
cm$c3[cm$c3>1] <-1
cm$c4[cm$c4>1] <-1
cm$c5[cm$c5>1] <-1
cm$c6[cm$c6>1] <-1
cm$c7[cm$c7>1] <-1
cm$c8[cm$c8>1] <-1
cm$c9[cm$c9>1] <-1
cm$c10[cm$c10>1] <-1
cm$c11[cm$c11>1] <-1
cm$c12[cm$c12>1] <-1
cm$c13[cm$c13>1] <-1
cm$c14[cm$c14>1] <-1
cm$c15[cm$c15>1] <-1
#calcualte total comorbidities
cm$totcm<-cm$c1 + cm$c2  + cm$c3  + cm$c4  + cm$c5  + cm$c6  + cm$c7  + cm$c8 
+ cm$c9 + cm$c10 + cm$c11 + cm$c12 + cm$c13 + cm$c14 + cm$c15
head(cm)
summary(cm$totcm)
table(cm$totcm)

#rename columns in comrobidity data frame
colnames(cm)<- c("id", "diab", "liver", "can", "metscan", "obese", "dem", "para", 
                 "stroke", "hyp", "renal", "mi", "copd", "chf", "pvd", "rh", "totcm")
cm_mgus<-cm
save(cm_mgus, file="N:\\DataAnalysis\\FormattedData\\ComorbidMGUS.Rdata")


#merge with patient data
mgus_cm<-left_join(mgus_pat, cm)
summarise(mgus_cm, n_distinct(id)) 

mgus_cm <- mgus_cm %>%
  mutate_at(vars(diab, liver, can, metscan, obese, dem, para, 
                 stroke, hyp, renal, mi, copd, chf, pvd, rh, totcm), 
            ~ replace_na(., 0))
mgus_cm %>% tally(diab)
#table of total no. comorbidities
cmtab<-table(mgus_cm$totcm)
t2<-prop.table(cmtab)
cbind(cmtab, t2)

#select comorbidity data only to create sums and %s
c1<- mgus_cm %>% 
  select("diab", "liver", "can", "metscan", "obese", "dem", "para", 
         "stroke", "hyp", "renal", "mi", "copd", "chf", "pvd", "rh")
#total with each condition
apply(c1, 2, sum)
#% with each condition
apply(c1, 2, mean)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Control cohort 
#Cancers have been excluded as a comorbidity in control group
rm(list=ls())

#comorbidity data for controls 
cm<- read.delim("N:\\ImportedData\\Incoming_2021-11-24\\LTH20102_Control_Prev_Diagnosis.txt", header=TRUE, sep="|")
names(cm)

#create indicator for each comorbidity
#Diabetes
cm$diab<- ifelse(
  grepl("E10.0|E10.1|E10.2|E10.3|E10.4|E10.5|E10.6|E10.7|E10.8|E10.9|
              E11.0|E11.1|E11.2|E11.3|E11.4|E11.5|E11.6|E11.7|E11.8|E11.9|
              E12.0|E12.1|E12.2|E12.3|E12.4|E12.5|E12.6|E12.7|E12.8|E12.9|
              E13.0|E13.1|E13.2|E13.3|E13.4|E13.5|E13.6|E13.7|E13.8|E13.9|
              E14.0|E14.1|E14.2|E14.3|E14.4|E14.5|E14.6|E14.7|E14.8|E14.9",
        cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(diab)

#Liver disease
cm$liver<- ifelse(
  grepl("B18.0|B18.1|B18.2|B18.8|B18.9|
              I85.0|I85.9|I86.4|I98.2|
              K70|K70.0|K70.1|K70.2|K70.3|K70.4|K70.9|
              K71|K71.1|K71.3|K71.4|K71.5|K71.7|
              K72|K72.1|K72.9|K73|K74|K76|
              K76.0|K76.2|K76.3|K76.4|K76.5|K76.6|K76.7|K76.8|K76.9|Z94.4",
        cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(liver)

#Previous Cancer (check if need to exclude C44 (14 counts))
#cm$can <- ifelse(
#          grepl("C",cm$Diag_Code, ignore.case = T), 1, 0)
#cm %>% tally(can)

#metastatic cancer
#cm$metscan <- ifelse(
#  grepl("C77|C78|C79|C80",cm$Diag_Code, ignore.case = T), 1, 0)
#cm %>% tally(metscan)

#Obesity
cm$obese<-ifelse(grepl("E66", cm$Diag_Code, ignore.case=T), 1, 0)
cm %>% tally(obese)

#Dementia
cm$dem <- ifelse(
  grepl("F00|F01|F02|F03|F05.1|G30|G31.1"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(dem)

#Hemiplegia or paraplegia
cm$para <- ifelse(
  grepl("G04.1|G11.4|G80.1|G80.2|G81|G82|G83.0|G83.1|G83.2|G83.4|G83.9"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(para)

#Stroke
cm$stroke<- ifelse(
  grepl("G45|G46|H34.0|I60|I61|I62|I63|I64|I65|I66|I67|I68|I69"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(stroke)

#hypertension
cm$hyp<- ifelse(
  grepl("I10|I11.9|I12.9|I13.9"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(hyp)

#renal disease
cm$renal<- ifelse(
  grepl("I12.0|I13.1|N03.2|N03.4|N03.5|N03.6|N03.7|
        N05.2|N05.3|N05.4|N05.5|N05.6|N05.7|N18|N19|N25.0|
        Z49.0|Z49.1|Z49.2|Z94.0|Z99.2"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(renal)

#Myocardial infarction
cm$mi<- ifelse(
  grepl("I21|I22|I25.2"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(mi)

#COPD
cm$copd<- ifelse(
  grepl("I27.8|I27.9|
        J40|J41|J42|J43|J4|J45|J46|J47|
        J60|J61|J62|J63|J64|J65|J66|J67|
        J68.4|J70.1|J70.3"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(copd)

#Congestive heart failure
cm$chf<- ifelse(
  grepl("I43|I50|I09.0|I11.0|I13.0|I13.2|I25.5|
        I42.0|I42.5|I42.6|I42.7|I42.8|I42.9|P29.0"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(chf)

#peripheral vascular disease
cm$pvd<- ifelse(
  grepl("I70|I71|I73.1|I73.8|I73.9|I77.1|I79.0|
        K55.1|K55.8|K55.9|Z95.8|Z95.9"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(pvd)

#Rheumatological conditioins
cm$rh<- ifelse(
  grepl("M05|M06|M31.5|M32|M33|M34|M35.1|M35.3|M36.0"
        ,cm$Diag_Code, ignore.case = T), 1, 0)
cm %>% tally(rh)

#aggregate by study Id to get totals for each id
cm<- cm %>%
  group_by(PseudoStudy_ID) %>%
  summarise(c1=sum(diab), c2=sum(liver), 
            c5=sum(obese), c6=sum(dem), c7=sum(para), c8=sum(stroke), 
            c9=sum(hyp), c10=sum(renal), c11=sum(mi), c12=sum(copd), 
            c13=sum(chf), c14=sum(pvd), c15=sum(rh)) 
#recode any vales>1 as 1 to get binary indicator for each ID
cm$c1[cm$c1>1] <-1
cm$c2[cm$c2>1] <-1
cm$c5[cm$c5>1] <-1
cm$c6[cm$c6>1] <-1
cm$c7[cm$c7>1] <-1
cm$c8[cm$c8>1] <-1
cm$c9[cm$c9>1] <-1
cm$c10[cm$c10>1] <-1
cm$c11[cm$c11>1] <-1
cm$c12[cm$c12>1] <-1
cm$c13[cm$c13>1] <-1
cm$c14[cm$c14>1] <-1
cm$c15[cm$c15>1] <-1
#calcualte total comorbidities
cm$totcm<-cm$c1 + cm$c2  + cm$c5  + cm$c6  + cm$c7  + cm$c8 
+ cm$c9 + cm$c10 + cm$c11 + cm$c12 + cm$c13 + cm$c14 + cm$c15
head(cm)
summary(cm$totcm)
table(cm$totcm)

#rename columns in comrobidity data frame
colnames(cm)<- c("id", "diab", "liver", "obese", "dem", "para", 
                 "stroke", "hyp", "renal", "mi", "copd", "chf", "pvd", "rh", "totcm")

cm_con<-cm
save(cm_con, file="N:\\DataAnalysis\\FormattedData\\ComorbidControl.Rdata")

#merge with patients with bloods - load myeloma saved files
load("N:\\DataAnalysis\\FormattedData\\Controls.RData")

#this needs to be updated as have removed cancers 10/2/22
#merge with patient data
con_cm<-left_join(c_pat, cm)
summarise(con_cm, n_distinct(id)) 

con_cm <- con_cm %>%
  mutate_at(vars(diab, liver, can, metscan, obese, dem, para, 
                 stroke, hyp, renal, mi, copd, chf, pvd, rh, totcm), 
            ~ replace_na(., 0))
con_cm %>% tally(diab)
#table of total no. comorbidities
cmtab<-table(con_cm$totcm)
t2<-prop.table(cmtab)
cbind(cmtab, t2)

#select comorbidity data only to create sums and %s
c1<- con_cm %>% 
  select("diab", "liver", "can", "metscan", "obese", "dem", "para", 
         "stroke", "hyp", "renal", "mi", "copd", "chf", "pvd", "rh")
#total with each condition
apply(c1, 2, sum)
#% with each condition
apply(c1, 2, mean)

