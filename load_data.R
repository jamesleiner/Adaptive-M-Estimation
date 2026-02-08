###################################################################################################
# Compile data from Osteoarthritis Initiative (OAI) for semi-synthetic experiments
#
# Download at https://www.niams.nih.gov/grants-funding/funded-research/osteoarthritis-initiative
#
###################################################################################################




library(haven)
library(dplyr)
library(mice)
library(ggplot2)
library(caret)
library(pROC)
library(nnet)
library(clustermq)
library(sandwich)
library(reshape2)
library(nnet)
library(dplyr)
library(MASS)
library(randomForest)



#Variables to include: demographic variables, changes in KL grade, etc
VARS_TO_INCLUDE = c("ID","P02HISP","P02RACE","P02SEX","V00AGE","P01BMI","V00XRKL","V00XRJSM","V00KOOSQOL","V06KOOSQOL","V06XRJSM","V06XRKL") 
VARS_RIGHT = c("V00RKFHDEG","P01PMRKRCV","V00WOMADLR","V00WOMKPR","P01RSXKOA","V06RKFHDEG","V06WOMADLR","V06WOMKPR","V99ERKRPCF")
VARS_LEFT = c("V00LKFHDEG","P01PMLKRCV","V00WOMADLL","V00WOMKPL","P01LSXKOA","V06LKFHDEG","V06WOMADLL","V06WOMKPL","V99ELKRPCF")


#Load SAS datasets, downloadable from OAI
enrollees <- read_sas("data/enrollees.sas7bdat")
kxr  <- read_sas("data/kxr_sq_bu00.sas7bdat")
kxr06  <- read_sas("data/kxr_sq_bu06.sas7bdat")
allclinical <- read_sas("data/allclinical00.sas7bdat")
allclinical06 <- read_sas("data/allclinical06.sas7bdat")
outcomes99 <- read_sas("data/outcomes99.sas7bdat")



#Filter, merge and save datasets
enrollees <-  enrollees[enrollees$V00COHORT %in% c(1,2),names(enrollees) %in% VARS_TO_INCLUDE]
kxr <- kxr[kxr$READPRJ ==15,]
kxr06 <- kxr06[kxr06$READPRJ ==15,]

kxr_right <- kxr[kxr$SIDE ==1,names(kxr)%in% c(VARS_TO_INCLUDE,VARS_RIGHT)]
kxr_left <- kxr[kxr$SIDE ==2,names(kxr)%in% c(VARS_TO_INCLUDE,VARS_LEFT)]

kxr06_right <- kxr06[kxr06$SIDE ==1,names(kxr06)%in% c(VARS_TO_INCLUDE,VARS_RIGHT)]
kxr06_left <- kxr06[kxr06$SIDE ==2,names(kxr06)%in% c(VARS_TO_INCLUDE,VARS_LEFT)]


allclinical_right <- allclinical[,names(allclinical)%in% c(VARS_TO_INCLUDE,VARS_RIGHT)]
allclinical_left <- allclinical[,names(allclinical)%in% c(VARS_TO_INCLUDE,VARS_LEFT)]


allclinical06_right <- allclinical06[,names(allclinical06)%in% c(VARS_TO_INCLUDE,VARS_RIGHT)]
allclinical06_left <- allclinical06[,names(allclinical06)%in% c(VARS_TO_INCLUDE,VARS_LEFT)]


outcomes_right <- outcomes99[,names(outcomes99)%in% c(VARS_TO_INCLUDE,VARS_RIGHT)]
outcomes_left <- outcomes99[,names(outcomes99)%in% c(VARS_TO_INCLUDE,VARS_LEFT)]



dataset_right <- merge(enrollees,kxr_right,by="ID")
dataset_right <- merge(dataset_right,allclinical_right,by="ID")
dataset_right <- merge(dataset_right,kxr06_right,by="ID")
dataset_right <- merge(dataset_right,allclinical06_right,by="ID")



dataset_left <- merge(enrollees,kxr_right,by="ID")
dataset_left <- merge(dataset_left,allclinical_left,by="ID")
dataset_left <- merge(dataset_left,kxr06_left,by="ID")
dataset_left <- merge(dataset_left,allclinical06_left,by="ID")
save(dataset_right, file="data/dataset_right.Rdata")
save(dataset_left, file="data/dataset_left.Rdata")
