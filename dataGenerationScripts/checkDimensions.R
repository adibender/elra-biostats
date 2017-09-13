## check dimensions
# patient
patient <- readRDS("../dataCurrent/patientVC.Rds")
assert_data_frame(patient, nrows=9661, ncols=34)
# ICU
icu <- readRDS("../dataCurrent/ICUvc.Rds")
assert_data_frame(icu, nrows=695, ncols=3)
# daily (time-dependent covariates)
daily <- readRDS("../dataCurrent/mergedAndCleanedDataVC.Rds")
assert_data_frame(daily, nrows=100559, ncols=58)
# transformed data with nutrition
ped <- readRDS("../dataCurrent/dataWithNutrition.Rds")
assert_data_frame(ped, nrows=219625, ncols=83)
# last observation ped data set
did <- readRDS("../dataCurrent/dataID.Rds")
assert_data_frame(did, nrows=9661, ncols=90)

assert_set_equal(unique(ped$CombinedID), unique(patient$CombinedID), ordered=TRUE)
assert_set_equal(did$CombinedID, patient$CombinedID, ordered=TRUE)

## same for data where hospital released was considered censoring event
ped.hosp <- readRDS("../dataCurrentHosp/dataWithNutritionHosp.Rds")
assert_data_frame(ped.hosp, nrows=157936, ncols=83)
