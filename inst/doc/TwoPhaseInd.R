### R code from vignette source 'TwoPhaseInd.Rnw'

###################################################
### code chunk number 1: loadLibrary
###################################################

library(TwoPhaseInd)


###################################################
### code chunk number 2: whiBioMarker
###################################################
data(whiBioMarker)

dim(whiBioMarker)
str(whiBioMarker)



###################################################
### code chunk number 3: spmleNonIndNoExtra
###################################################
spmleNonIndNoExtra <- spmle(data=whiBioMarker,  ## dataset
                 response="stroke",	## response variable
                 treatment="hrtdisp", ## treatment variable
                 BaselineMarker="papbl",	## environment variable
                 extra=NULL,
                 phase="phase",	## phase indicator (1 and 2)
                 ind=FALSE	## independent or non-indepentent
)

spmleNonIndNoExtra


###################################################
### code chunk number 4: spmleIndNoExtra
###################################################
spmleIndNoExtra <- spmle(data=whiBioMarker,	## dataset
              response="stroke", ## response variable
              treatment="hrtdisp",	## treatment variable
              BaselineMarker="papbl",		## environment variable
              extra=NULL,
              phase="phase",	## phase indicator
              ind=TRUE	## independent or non-indepentent
)

spmleIndNoExtra


###################################################
### code chunk number 5: spmleNonIndExtra
###################################################
spmleNonIndExtra <- spmle(data=whiBioMarker,	## dataset
               response="stroke",	## response variable
               treatment="hrtdisp",	## treatment variable
               BaselineMarker="papbl",	## environment variable
               extra=c(
                       "age" 		## age
                       ## physical activity levels
                        , "dias" 	## diabetes
                        , "hyp" 	## hypertension
                        , "syst" 	## systolic
                        , "diabtrt"	## diastolic BP
                        , "lmsepi" ## waist:hip ratio
                            ),	## extra variable(s)
               phase="phase",	## phase indicator
               ind=FALSE	## independent or non-indepentent
)

spmleNonIndExtra


###################################################
### code chunk number 6: spmleIndExtra
###################################################
spmleIndExtra <- spmle(data=whiBioMarker,	## dataset
            response="stroke",	## response variable
            treatment="hrtdisp",	## treatment variable
            BaselineMarker="papbl",	## environment variable
            extra=c(
               "age" 	## age
                		## physical activity levels
              , "dias"	## diabetes
              , "hyp" ## hypertension
              , "syst" ## systolic
              , "diabtrt"	## diastolic BP
              , "lmsepi" ## waist:hip ratio
                 ),	## extra variable(s)
            phase="phase", ## phase indicator
            ind=TRUE ## independent or non-indepentent
)

spmleIndExtra


###################################################
### code chunk number 7: melIndNoExtra
###################################################
melIndNoExtra <- mele(data=whiBioMarker,	## dataset
            response="stroke", ## response variable
            treatment="hrtdisp",	## treatment variable
            BaselineMarker="papbl",	## environment variable
            extra=NULL,
            phase="phase",	## variable for phase indicator
            ind=TRUE ## independent or non-indepentent
)
melIndNoExtra


###################################################
### code chunk number 8: melNoIndNoExtra
###################################################
melNoIndNoExtra <- mele(data=whiBioMarker,	## dataset
              response="stroke",	## response variable
              treatment="hrtdisp",	## treatment variable
              BaselineMarker="papbl",	## environment variable
              extra=NULL,
              phase="phase", ## phase indicator
              ind=FALSE	## independent or non-indepentent
)
melNoIndNoExtra


###################################################
### code chunk number 9: melIndExtra
###################################################
melIndExtra <- mele(data=whiBioMarker,	## dataset
          response="stroke",	## response variable
          treatment="hrtdisp",		## treatment variable
          BaselineMarker="papbl",		## environment variable
          extra=c(
             "age" 	## age
                		## physical activity levels
              , "dias" 	## diabetes
              , "hyp" ## hypertension
              , "syst" 	## systolic
              , "diabtrt"	## diastolic BP
              , "lmsepi" ## waist:hip ratio
              ),	## extra variable(s)
          phase="phase",	## phase indicator
          ind=TRUE	## independent or non-indepentent
)
melIndExtra


###################################################
### code chunk number 10: melNoIndExtra
###################################################
melNoIndExtra <- mele(data=whiBioMarker,	## dataset
            response="stroke",	## response variable
            treatment="hrtdisp",	## treatment variable
            BaselineMarker="papbl",	## environment variable
            extra=c(
                "age"	## age
                		## physical activity levels
                , "dias" 	## diabetes
                , "hyp" 	## hypertension
                , "syst" 	## systolic
                , "diabtrt"	## diastolic BP
                , "lmsepi" ## waist:hip ratio
                ),	## extra variable(s)
            phase="phase",	## phase indicator
            ind=FALSE	## independent or non-indepentent
)
melNoIndExtra


###################################################
### code chunk number 11: data
###################################################
data(acodata)

dim(acodata)
str(acodata)



###################################################
### code chunk number 12: rfit0
###################################################
rfit0 <- acoarm(data=acodata,
                 svtime="vacc1_evinf",
                 event="f_evinf",
                 treatment="f_treat",
                 BaselineMarker="fcgr2a.3",
                 id="ptid",
                 subcohort="subcoh",
                 esttype=1,
                 augment=0,
                 extra=c("f_agele30","f_hsv_2","f_ad5gt18","f_crcm","any_drug","num_male_part_cat","uias","uras")) 
rfit0


###################################################
### code chunk number 13: rfit1
###################################################
rfit1 <- acoarm(data=acodata,
                 svtime="vacc1_evinf",
                 event="f_evinf",
                 treatment="f_treat",
                 BaselineMarker="fcgr2a.3",
                 id="ptid",
                 subcohort="subcoh",
                 esttype=1,
                 augment=1,
                 extra=c("f_agele30","f_hsv_2","f_ad5gt18","f_crcm","any_drug","num_male_part_cat","uias","uras")) 
rfit1


###################################################
### code chunk number 14: rfit2
###################################################
rfit2 <- acoarm(data=acodata,
                 svtime="vacc1_evinf",
                 event="f_evinf",
                 treatment="f_treat",
                 BaselineMarker="fcgr2a.3",
                 id="ptid",
                 subcohort="subcoh",
                 esttype=1,
                 augment=2,
                 extra=c("f_agele30","f_hsv_2","f_ad5gt18","f_crcm","any_drug","num_male_part_cat","uias","uras")) 
rfit2


###################################################
### code chunk number 15: sessionInfo
###################################################
sessionInfo()


