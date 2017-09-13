## load data
patient    <- readRDS(paste0(save.path.data, "patientVC.Rds"))
daily      <- readRDS(paste0(save.path.data, "mergedAndCleanedDataVC.Rds"))
icu        <- readRDS(paste0(save.path.data, "ICUvc.Rds"))
nutri      <- readRDS(paste0(save.path.data, "dataWithNutrition.Rds"))
nutri.hosp <- readRDS(paste0(save.path.data.hosp, "dataWithNutritionHosp.Rds"))

## change reference level of admission category (otherwise problems in
# coefficient plots when subgroups are displayed).
nutri$AdmCatID      <- relevel(nutri$AdmCatID, ref="Surgical Elective")
nutri.hosp$AdmCatID <- relevel(nutri.hosp$AdmCatID, ref="Surgical Elective")

### create by dummy for random effect (can be set to zero for prediction)
nutri$icuByDummy      <- 1
nutri.hosp$icuByDummy <- 1

## for description (Oral Intake as one vector variable)
daily$OI <- factor(daily$OralIntake,
  labels = c("No Oral Intake", "Oral Intake"))
anyOI <- aggregate(OralIntake ~ CombinedID, data = daily, function(z)
            any(z == "Yes"))
daily$OI[daily$CombinedID %in% anyOI$CombinedID[anyOI$OralIntake]] <- "Oral Intake"

## create Status variable in daily data set (for labels in graphics)
daily$Status <- factor(daily$PatientDied, labels=c("Survived", "Deceased"))
## define Study_Day 1 (calendar day 1 at ICU) as day 0, days 2-12 as days 1-11
daily$Study_Day <- daily$Study_Day - 1

## create data set that contains last observation per patient only
# usefull for descriptive analyses
data.id <- nutri[get_last(nutri$CombinedID), ]
data.id$intMax <- with(nutri, tapply(intmid, CombinedID, max))
# create integer variable to be used in lasagna plots -> cleanNutri.Rnw ->
# descripts.Rnw -> plotFunctionsNutrition.R@plot_lasagna
data.id$AdequacyCalsTotF <-  with(data.id, AdequacyCalsTot0to30 +
  2*(AdequacyCalsTot30To70) + 3*(AdequacyCalsTotAbove70))
# same for last day imputed vars
data.id$AdequacyCals2TotF <-  with(data.id, AdequacyCals2Tot0to30 +
  2*(AdequacyCals2Tot30To70) + 3*(AdequacyCals2TotAbove70))

# protocols had been imputed if patient was discharged from ICU before end of
# logged calendar days.
# We imputed these days, but for description of actually observed calories we
# set these imputed observations back to NA
maxdays.nutri <- ncol(data.id$DaysMat)
setNA <- matrix(1:maxdays.nutri, ncol=maxdays.nutri, nrow=nrow(data.id),
  byrow=TRUE) > (data.id[,'dayMax'] - 1)

data.id$AdequacyCalsTotF[setNA]   <- NA
data.id$AdequacyCals2TotF[setNA]  <- NA

# adjusted caloriesPercentage
data.id$calPerc                        <- data.id$caloriesPercentage[, -1]
data.id$calPerc[data.id$calPerc > 200] <- 200
data.id$calPerc[setNA]                 <- NA
data.id$calIntake                      <- data.id$caloriesIntake[, -1]
data.id$calIntake[setNA]               <- NA
data.id$calPerKG                       <- with(data.id, calIntake/Weight)

data.id$Status <- factor(data.id$PatientDied, labels=c("Survived", "Deceased"))

## check dimensions
assert_data_frame(data.id,    nrows = 9717,   ncols = 90)
assert_data_frame(patient,    nrows = 9717,   ncols = 34)
assert_data_frame(daily,      nrows = 101145, ncols = 58)
assert_data_frame(icu,        nrows = 695,    ncols = 3)
assert_data_frame(nutri,      nrows = 220891, ncols = 83)
assert_data_frame(nutri.hosp, nrows = 157936, ncols = 83)


saveRDS(data.id,    paste0(save.path.data, "dataID.Rds"))
saveRDS(patient,    paste0(save.path.data, "patientVC.Rds"))
saveRDS(daily,      paste0(save.path.data, "mergedAndCleanedDataVC.Rds"))
saveRDS(icu,        paste0(save.path.data, "ICUvc.Rds"))
saveRDS(nutri,      paste0(save.path.data, "dataWithNutrition.Rds"))
saveRDS(nutri.hosp, paste0(save.path.data.hosp, "dataWithNutritionHosp.Rds"))
