## data to adjust
patient <- readRDS("../data/patient.Rds")
# daily   <- readRDS(paste0(save.path.data, "daily.Rds"))
mc      <- readRDS("../data/mergedAndCleanedData.Rds")
icu     <- readRDS("../data/ICU.Rds")

# data containing valid cases
valid.cases <- readRDS(paste0(save.path.data, "validCases.Rds"))
id.valid    <- unique(valid.cases$CombinedID)
icu.valid   <- unique(valid.cases$CombinedicuID)

# save adjusted data
saveRDS(subset(patient, CombinedID %in% id.valid),
  paste0(save.path.data, "patientVC.Rds"))
# saveRDS(subset(daily, CombinedID %in% id.valid),
#   paste0(save.path.data, "dailyVC.Rds"))
saveRDS(subset(mc, CombinedID %in% id.valid),
  paste0(save.path.data, "mergedAndCleanedDataVC.Rds"))
saveRDS(subset(icu, CombinedicuID %in% icu.valid),
  paste0(save.path.data, "ICUvc.Rds"))
