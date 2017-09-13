## get raw data
data.long <- readRDS("../data/mergedAndCleanedData.Rds")
patient   <- readRDS("../data/patient.Rds")
ICU       <- readRDS("../data/ICU.Rds")

# Define time-varying variables (exclude ID and Study_Day)
timevarying <- c("EN", "PN", "OralIntake", "Propofol", "PropofolCal",
	"caloriesIntake", "proteinIntakeAdjusted", "caloriesPercentage",
	"caloriesPercentage2", "proteinAdjustedPercentage", "proteinGproKG", "inMV")

data.pem <- transform_to_ped(data.long, timevarying, patient, ICU, brks,
	survvar="Survdays", calSurvVar="calendarSurvdays", max.follow=max.follow)

data.long <- subset(data.long, !is.na(data.long[["survhosp"]]))
data.pem.hosp <- transform_to_ped(data.long, timevarying, patient, ICU, brks,
	survvar="survhosp", calSurvVar="calendarSurvhosp", max.follow=max.follow)

## check dimensions
assert_data_frame(data.pem,      nrows = 328513, ncols = 64)
assert_data_frame(data.pem.hosp, nrows = 232713, ncols = 64)

saveRDS(data.pem,      file = paste0(save.path.data, "dataPEM.Rds"))
saveRDS(data.pem.hosp, file = paste0(save.path.data.hosp, "dataPEMhosp.Rds"))
