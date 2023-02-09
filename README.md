# MyelomaPrediction
Supporting code for the manuscript

**Development and internal validation of a risk prediction model to identify myeloma based on routine blood tests: a case-control study**
Lesley Smith, Jonathan Carmichael, Gordon Cook, Bethany Shinkins, Richard D Neal  
Cancers 2023, 15(3), 975; https://doi.org/10.3390/cancers15030975

Code written by Lesley Smith 

The study used de-identified data from Leeds Teaching Hospitals Trust (LTHT). The R code used for data formatting, cleaning and analysis is incuding in the following R scripts: 

1.	FormatPatientData.R Patient level data cleaning and formatting including comorbidities
2.	FormatBloodTest.R  Blood test data cleaning and formatting all tests in long form
3.	FormatCXDara.R Extracting data to include in cross sectional modelling
4.	PatientDescriptiveStats.R Patient descriptive statistic
5.	BloodTestSummaryStats.R descriptive statistics on number and timing of blood tests
6.	BloodTestsPlots_90datIntervals.R Blood test results trajectory plots based on means in 90 day intervals
7.	DescStatsUnivariableModels.R Descriptive stats and univariable models
8.	ImputationModel.R Main model based on imputed data
9.	ImputationModelwithMGUSasControls Imputation model with MGUS cohort included as controls
10.	SensSpecThresholds.R Diagnostic accuracy statistics with different prevalence estimates 
