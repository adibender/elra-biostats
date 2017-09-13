
valid.cases <- readRDS(paste0(save.path.data, "validCases.Rds"))
vc.hosp     <- readRDS(paste0(save.path.data.hosp, "validCasesHosp.Rds"))


# variables for which first entry/column will be removed
# (as calendar day of admission is not a full 24h period and protocol may be unreliable)
vars.ignore.first <- c(
  "LHartlDynf", "LSimple460f", "L", "Lf",
  "OralIntakeMat", "inMVmat", "DaysMat", "intMat",
  "AdequacyCalsTot0to30", "AdequacyCalsTot30To70", "AdequacyCalsTotAbove70",
  "AdequacyCals2Tot0to30", "AdequacyCals2Tot30To70", "AdequacyCals2TotAbove70",
  "v3AdequacyCalsTot0to30", "v3AdequacyCalsTot30To70", "v3AdequacyCalsTotAbove70")

nutri <- createNutritionVars(
  valid.cases,
  brks         = brks,
  ignore.first = TRUE,
  vars.ignore  = vars.ignore.first)

nutri.hosp <- createNutritionVars(
  vc.hosp,
  brks         = brks,
  ignore.first = TRUE,
  vars.ignore  = vars.ignore.first)

## check dimensions
assert_data_frame(nutri,      nrows = 220891, ncols = 82)
assert_data_frame(nutri.hosp, nrows = 157936, ncols = 82)

saveRDS(nutri,      file = paste0(save.path.data, "dataWithNutrition.Rds"))
saveRDS(nutri.hosp, file = paste0(save.path.data.hosp, "dataWithNutritionHosp.Rds"))
