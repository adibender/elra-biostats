nutri      <- readRDS(paste0(save.path.data, "dataWithNutrition.Rds"))
nutri.hosp <- readRDS(paste0(save.path.data.hosp, "dataWithNutritionHosp.Rds"))

## define variables that should be kept in the data set for modeling purposes
needed.vars <- c(
  "CombinedID", "event", "intmid", "survtime", "tstart", "tend",
  "ApacheIIScore", "Age", "BMI", "Year", "DiagID2", "AdmCatID", "Gender",
  "freqPNd2to4", "freqPFd2to4", "freqOId2to4", "freqInMVd2to4",
  "intMat", "DaysMat", "incompleteProtocol", "calendarDaysInICU",
  "AdequacyCalsTot0to30", "AdequacyCalsTot30To70", "AdequacyCalsTotAbove70",
  "AdequacyCals2Tot0to30", "AdequacyCals2Tot30To70", "AdequacyCals2TotAbove70",
  "v3AdequacyCalsTot0to30", "v3AdequacyCalsTot30To70", "v3AdequacyCalsTotAbove70",
  "LHartlDynf", "LSimple460f",
  "CombinedicuID", "icuByDummy",
  "offset")

rmValidCases      <- readRDS(paste0(save.path.data, "rmValidCases.Rds"))
nutri.small       <- nutri[, needed.vars]
rmValidCases$rmNA <- length(unique(nutri.small$CombinedID))

daily   <- readRDS(paste0(save.path.data, "mergedAndCleanedDataVC.Rds"))
data.id <- readRDS(paste0(save.path.data, "dataID.Rds"))

## remove not existent factors
nutri.small <- droplevels(nutri.small)
nutri.small <- na.omit(nutri.small)
saveRDS(nutri.small, file = paste0(save.path.data, "nutriOrigSmall.Rds"))

rmValidCases$rmNA <- rmValidCases$rmNA - length(unique(nutri.small$CombinedID))
saveRDS(rmValidCases, "../dataCurrent/rmValidCases.Rds")


## subset other data sets to match patients used for modeling
patient <- readRDS(paste0(save.path.data, "patientVC.Rds"))
daily   <- readRDS(paste0(save.path.data, "mergedAndCleanedDataVC.Rds"))
icu     <- readRDS(paste0(save.path.data, "ICUvc.Rds"))

nutri   <- subset(nutri, CombinedID %in% nutri.small$CombinedID)
nutri   <- droplevels(nutri)
patient <- subset(patient, CombinedID %in% nutri.small$CombinedID)
patient <- droplevels(patient)
daily   <- subset(daily, daily$CombinedID %in% nutri.small$CombinedID)
daily   <- droplevels(daily)
icu     <- subset(icu, CombinedicuID %in% nutri.small$CombinedicuID)
icu     <- droplevels(icu)
data.id <- subset(data.id, CombinedID %in% nutri.small$CombinedID)
data.id <- droplevels(data.id)


saveRDS(data.id, paste0(save.path.data, "dataID.Rds"))
saveRDS(patient, paste0(save.path.data, "patientVC.Rds"))
saveRDS(daily,   paste0(save.path.data, "mergedAndCleanedDataVC.Rds"))
saveRDS(icu,     paste0(save.path.data, "ICUvc.Rds"))
saveRDS(nutri,   paste0(save.path.data, "dataWithNutrition.Rds"))


data.list  <- list(nutri.small)
names(data.list)  <- c("nutriOrigSmall")

nutri.hosp <- nutri.hosp[, needed.vars]
for(i in seq_along(data.list)) {

  tmp.i <- subset(nutri.hosp, CombinedID %in% unique(data.list[[i]]$CombinedID))
  tmp.i <- droplevels(subset(tmp.i, intmid %in% data.list[[i]]$intmid))
  saveRDS(tmp.i, paste0(save.path.data.hosp, names(data.list)[i], ".Rds"))

}
