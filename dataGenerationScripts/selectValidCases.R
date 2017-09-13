valid.cases <- readRDS(paste0(save.path.data, "dataPEM.Rds"))
vc.hosp     <- readRDS(paste0(save.path.data.hosp, "dataPEMhosp.Rds"))
daily       <- readRDS("../data/mergedAndCleanedData.Rds")

## set parameters for exclusion
# parameters regarding time scale of survival (days after ICU admission)
mindays.surv  = 5 # minimal number of days after ICU admission survived
mindays.icu   = 4 # minimal number of days spent on ICU after ICU admission
# parameters regarding time scale of nutrition (calendar days/logged protocol days)
mindays.nutri = 4L # minimal number of calendar/protocol days with certain nutrition
skip.days     = 1L # days to skip when calculating certain nutrition information
maxdays.nutri = ncol(valid.cases$DaysMat) # maximal number of logged
# calendar/protocol days


## create list that keeps track of excluded patient numbers
rm.df <- data.frame(
  BeforeSelection = NA, LifeTooShort   = NA, StayIsTooShort = NA,
  NoENandPNbutOI  = NA, OIonlyAfterExt = NA, NoMV1to4       = NA,
  AgeTooYoung     = NA, BMITooSmall    = NA)
# number of patients before application of exclusion criteria
rm.df$BeforeSelection <- length(unique((valid.cases$CombinedID)))

## apply following exclusion criteria
# - survival time shorter than 4 days or discharged alive from ICU within 4 days
# - nutrition protocol before the end of calendar day 12 while patient still on ICU
# (we consider this to be missing data)
# - neither enteral nor parenteral nutrition during the first 4 days of protocol,
# but additional oral intake during first 4 days
# - no mechenical ventilation during first 4 days
# - patients that have not reached age of at least 18

## since our analysis assumes a lag of 4 protocol days we are interested only in
## patients that were alive in interval (4, 5]
# sort valid.cases by ID and interval
valid.cases <- valid.cases[order(valid.cases$CombinedID, valid.cases$intmid), ]
## extract patients with Survival < 5 days after ICU admission
# length of observation periods each patient survived
id.rle <- rle(valid.cases$CombinedID)
# IDs of patients with survival < 5 days
life.is.too.short  <- id.rle$values[id.rle$lengths < mindays.surv]
rm.df$LifeTooShort <- length(life.is.too.short)
valid.cases        <- subset(valid.cases, !(CombinedID %in% life.is.too.short))

## similarly, patients that have not been at ICU for at least 4 days, will
# be excluded, as nutrition information was only logged on ICU
stay.is.too.short    <- unique(subset(valid.cases, DaysInICU <= mindays.icu)$CombinedID)
rm.df$StayIsTooShort <- length(stay.is.too.short)
valid.cases          <- subset(valid.cases, !(CombinedID %in% stay.is.too.short))

# remove if patient on ICU for more days than have been logged
# maxdays.protocol gives maximum number of logged days (for all patients),
# dayMax the number of logged days (per patient)
icu.after.protocol.end    <- unique(subset(valid.cases,
  pmin(calendarDaysInICU, maxdays.nutri) > dayMax)$CombinedID)
rm.df$ICUafterProtocolEnd <- length(icu.after.protocol.end)

valid.cases <- subset(valid.cases, !(CombinedID %in% icu.after.protocol.end))


## apply exclusion criteria based on protocol days
# first exclude allready excluded patients from daily data set
daily <- filter(daily, CombinedID %in% unique(valid.cases$CombinedID))
daily <- mutate(daily, lastDay = Study_Day==maxday)
# data frame containing booleans for exclusion criteria
# - no en or pn on calendar days 1-4 at ICU , but oral intake on at least one of
# them -
# no mechanical ventilation on calendar days 1-4 at ICU
exclude.d1to4 <- daily %>%
  group_by(CombinedID) %>%
  filter(Study_Day %in% seq_len(mindays.nutri)) %>%
  summarize(
    exclude.noenpn.oi = all(EN=="No") & all(PN=="No") & any(OralIntake=="Yes"),
    exclude.nomv4 = sum(inMV)==0 | all(DaysMechVent==0))

rm.df$NoENandPNbutOI <- sum(exclude.d1to4$exclude.noenpn.oi)
rm.df$NoMV1to4 <- sum(exclude.d1to4$exclude.nomv4[!exclude.d1to4$exclude.noenpn.oi])
# remove from data set
daily <- anti_join(daily, filter(exclude.d1to4, exclude.noenpn.oi | exclude.nomv4))

valid.cases <- subset(valid.cases, CombinedID %in% unique(daily$CombinedID))

## create nutriton information that will be used in models
# exclude day1 (i.e. day of ICU admission), as it is usually not a full 24h
# period
freq.d2to4 <- daily %>% group_by(CombinedID) %>%
  filter(Study_Day %in% (1+skip.days):mindays.nutri) %>%
  summarize(
    freqENd2to4   = sum(EN=="Yes"),
    freqPNd2to4   = sum(PN=="Yes"),
    freqPFd2to4   = sum(Propofol=="Yes"),
    freqOId2to4   = sum(OralIntake=="Yes"),
    freqInMVd2to4 = sum(inMV))

valid.cases <- merge(valid.cases, freq.d2to4, by="CombinedID")
vc.hosp     <- merge(vc.hosp, freq.d2to4, by="CombinedID")
valid.cases <- valid.cases[order(valid.cases$CombinedID, valid.cases$intmid), ]
vc.hosp     <- vc.hosp[order(vc.hosp$CombinedID, vc.hosp$intmid), ]

## general exclusion criteria
# age > 18
# BMI > 13
age.too.young     <- unique(valid.cases$CombinedID[valid.cases$Age < 18])
valid.cases       <- subset(valid.cases, !(CombinedID %in% age.too.young))
rm.df$AgeTooYoung <- sum(!is.na(age.too.young))

# remove patients with extremely low BMI
bmi.too.small     <- unique(valid.cases$CombinedID[valid.cases$BMI < 13])
valid.cases       <- subset(valid.cases, !(CombinedID %in% bmi.too.small))
rm.df$BMITooSmall <- sum(!is.na(bmi.too.small))

## since we are only interested in modeling from day 5 on,
## remove all observations for day 1-4
valid.cases <- subset(valid.cases, survtime > (mindays.surv - 1))
vc.hosp     <- subset(vc.hosp, CombinedID %in% unique(valid.cases$CombinedID))
vc.hosp     <- subset(vc.hosp, int %in% unique(valid.cases$int))

## check dimensions
assert_data_frame(valid.cases, nrows = 220891, ncols = 69)
assert_data_frame(vc.hosp,     nrows = 157936, ncols = 69)

## save processed data
saveRDS(valid.cases, file = paste0(save.path.data, "validCases.Rds"))
saveRDS(vc.hosp,     file = paste0(save.path.data.hosp, "validCasesHosp.Rds"))
saveRDS(rm.df,       file = paste0(save.path.data, "rmValidCases.Rds"))
