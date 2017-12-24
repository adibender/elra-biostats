# This code was added for the second revision
library(dplyr)
library(purrr)
library(mgcv)

## read data
nutri <- readRDS("../dataCurrent/nutriOrigSmall.Rds")
seq_fun <- function(n) {
  if(n > 11) {
    c(1:11, rep(11, n - 11))
  } else {
    seq_len(n)
  }
}
adeq_vars <- c("AdequacyCalsTot0to30", "AdequacyCalsTot30To70", "AdequacyCalsTotAbove70")
nutri_id <- nutri[!duplicated(nutri$CombinedID), adeq_vars]
# number of intervals per subject (id)
n_id <- nutri %>%
  select(CombinedID) %>%
  group_by(CombinedID) %>%
  summarize(n = n())
# running index for each subject
seq_id <- purrr:::map(n_id$n, seq_fun)

# create sequence index, that repeats last observation until the last
# observation of each subject
seq_nutri <- function(adeq_mat, seq_list) {
  n_cat <- map(seq_len(nrow(adeq_mat)), ~cumsum(adeq_mat[.,]))
  map2(n_cat, seq_list, ~.x[.y]) %>% unlist()
}

## create "simple" TD nutrition variable for sensitivity analysis suggested
# by reviewer
CI_seq <- seq_nutri(nutri_id$AdequacyCalsTot0to30, seq_id)
CII_seq <- seq_nutri(nutri_id$AdequacyCalsTot30To70, seq_id)
CIII_seq <- seq_nutri(nutri_id$AdequacyCalsTotAbove70, seq_id)
# note: use lag of 1 as data starts at (4,5]. To be comparable to
# the lag-lead window of the main analysis, however, nutrtion must not start
# affect hazard before (5,6]
# We cannot set the effect of nutrtion in (4,5] to zero using the simple model
# however, thus we are forced to set the covariate to 0 (default = 0)
# NOTE: There is a general problem with this definition of the nutrition variable
# We do not know in advance if a person survived 11 days, thus for those who
# did, when calculating the % of days in category CI nutrtion, we use future
# information, for those who didn't, it is not clear what 2 of 11 days with CI
# nutrition means
nutri$adeqTot0to30   <- lag(CI_seq, default = 0)
nutri$adeqTot30To70  <- lag(CII_seq, default = 0)
nutri$adeqTotAbove70 <- lag(CIII_seq, default = 0)


## alternative sensitivity analysis. Here we use # of days with CI, CII and CIII
# within the dynamic time-window used in the main analysis, respectively
nadeq   <- rowSums(nutri$LHartlDynf)
CI_ll   <- rowSums(nutri$AdequacyCalsTot0to30*nutri$LHartlDynf)
CII_ll  <- rowSums(nutri$AdequacyCalsTot30To70*nutri$LHartlDynf)
CIII_ll <- rowSums(nutri$AdequacyCalsTotAbove70*nutri$LHartlDynf)
head(cbind(CI_ll, CII_ll, CIII_ll, nadeq))

nutri$LLadequTot0to30 <- CI_ll/nadeq
nutri$LLadequTot30To70<- CII_ll/nadeq
nutri$LLadequTotAbove70 <- CIII_ll/nadeq
nutri_LLadequ <- nutri[, c("LLadequTot0to30", "LLadequTot30To70", "LLadequTotAbove70")]
summary(nutri_LLadequ)
head(nutri_LLadequ)
## for first interval (4,5], where effects are set to 0 due to the definition
# of the lag-lead window, there are 0 eligible nutrition contributions,
# therefore, in the previous step, deviding by 0 created NAs
# For this simple model, we set them to 0 to be able to calculate the model
# however, it will be treated as a nutrition contribution of 0 by the model
nutri$LLadequTot0to30[is.na(nutri$LLadequTot0to30)] <- 0
nutri$LLadequTot30To70[is.na(nutri$LLadequTot30To70)] <- 0
nutri$LLadequTotAbove70[is.na(nutri$LLadequTotAbove70)] <- 0
saveRDS(nutri, "../dataCurrent/nutri2.Rds")

## model formula for models suggested by reviewer
# (only lag time considered and definition of TDC has some caveats,
# see comments above)
form <- event ~ s(intmid, bs = "ps") +
  ApacheIIScore + ApacheIIScore:intmid +
  s(Age, by = intmid, bs = "ps") +
  s(BMI, by = intmid, bs = "ps") +
  Year + DiagID2 + AdmCatID + Gender +
  freqInMVd2to4 + freqPFd2to4 +  freqOId2to4 + freqPNd2to4 +
  s(CombinedicuID, bs = "re", by=icuByDummy) +
  adeqTot0to30 +
  adeqTot30To70 +
  adeqTotAbove70

## Alternative "simple" model (only use % of CI, CII, CIII nutrition within
# dynamic time window from main application)
# Same caveats as above :
# 1) covariate value has to be set to 0 for interval (4,5])
# 2) starting with interval (16,17] covariate info remains constant and thus
# hazard ratio between two patients with different nutrition profiles remains
# constant throughout the rest of the follow up (c.p)
form2 <- event ~ s(intmid, bs = "ps") +
  ApacheIIScore + ApacheIIScore:intmid +
  s(Age, by = intmid, bs = "ps") +
  s(BMI, by = intmid, bs = "ps") +
  Year + DiagID2 + AdmCatID + Gender +
  freqInMVd2to4 + freqPFd2to4 +  freqOId2to4 + freqPNd2to4 +
  s(CombinedicuID, bs = "re", by=icuByDummy) +
  LLadequTot30To70 +
  LLadequTotAbove70


mod_simple <- gam(form, data = nutri, family = poisson(), offset = offset,
  method = "REML", control=list(trace=TRUE))
saveRDS(mod_simple, "simple_tdc_mod.Rds")


mod_simple2 <- gam(form2, data = nutri, family = poisson(), offset = offset,
  method = "REML", control=list(trace=TRUE))
saveRDS(mod_simple2, "simple_tdc_modLLadequ.Rds")

mod_no_nutri <- update(mod_simple, .~. - adeqTot0to30 - adeqTot30To70 -
 adeqTotAbove70)
saveRDS(mod_no_nutri, "mod_no_nutrition.Rds")
