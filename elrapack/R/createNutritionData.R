##' Creates nutrition adequacy variables and Lag/Lead matrices

#' @param valid.cases Data set in piece-wise exponential data format.
#' @param brks Breaks used to cut the time line for piece-wise exponential format.
#' @param ignore.first Indicates whether first observation of the TDC
# covariates should be ignored.
#' @param vars.ignore List of variable names specifying TDC for which first
#' column should be removed.
#' @import checkmate
#' @keywords internal
#' @export
createNutritionVars <- function(
  valid.cases,
  brks         = c(0:12, seq(15, 55, by = 5), 61.1),
  ignore.first = TRUE,
  vars.ignore  = NULL) {

  # check inputs
  assert_data_frame(valid.cases)
  assert_numeric(brks, lower=0, upper=ceiling(max(brks)), finite=TRUE,
    any.missing=FALSE, unique=TRUE)
  assert_flag(ignore.first)
  if(!is.null(vars.ignore)) {
    assert_character(vars.ignore, any.missing=FALSE, unique=TRUE)
  }

  # reassure correct ordering of data set by ID and time
  valid.cases <- valid.cases[order(valid.cases$CombinedID, valid.cases$intmid), ]
  ## rle.id contains to elements, ID and the number of observations/intervals per ID
  rle.id <- rle(valid.cases$CombinedID)
  ## seq.id contains interval index per patient
  seq.id <- unlist(lapply(rle.id$lengths, seq_len))
  stopifnot(length(seq.id)==nrow(valid.cases))

  ## get number of maximally logged calendar/prtocol days
  maxdays.nutri = ncol(valid.cases$DaysMat)

  ## Relevant Oral Intake = OralIntake and NO mech. vent. on the same calendar day
  ## Here beginMechVent and EndMechVent refer to nutrition level i.e. calendar days
  ## Define days from calendar day of begin MV (if < 12) until calendar day of MV (or 12)
  ## as inMV (per patient)
  inMV <- t(
    apply(
      X = cbind(
        begin = valid.cases$calendarBeginMechVent,
        end   = valid.cases$calendarEndMechVent),
      MARGIN = 1,
      function(x) {
        ret <- rep(0, maxdays.nutri)
        if (x["begin"] < maxdays.nutri) {
          ret[max(1, x["begin"]):min(maxdays.nutri, x["end"])] <- 1
        }
        ret
      }))

  stopifnot(nrow(inMV)==nrow(valid.cases))

  ## Define Oral Intake matrix. Needed later b/c days with MV will be
  # set to the next higher category of nutrition,
  # e.g. caloric adequacy < 30% + Oral Intake -> caloric adequacy 30-70 %
  # However, we assume next higher category of caloric intake if additional oral
  # intake occured only when there was no mechanical ventilation on the same day
  OralIntakeMat <- (valid.cases$OralIntake == "Yes" & !inMV) * 1
  NoOralIntakeMat <- (!OralIntakeMat) * 1

  stopifnot(nrow(OralIntakeMat) == nrow(valid.cases))

  # add newly created variables to data set
  valid.cases$OralIntakeMat <- OralIntakeMat
  valid.cases$inMVmat       <- inMV


  ############################ CALORIES ######################################

  # As we consider 3 categories of caloric adequacy <30%, 30-70%, >=70%
  # We create new variables for each of these categories
  # For later use these need to be Matrices with
  # ncol = number of days with maximal number of logged calendar/protocol days
  # i.e. the maximal number of days for which nutrition information had been
  # recorded (here 12)
  # nrow = equal number of rows in valid.cases, for each patient the same protocol
  # (from days 1-12) repeats for each row of the patient, i.e. for each
  # interval a patient was alive
  # Days with incomplete protocols, for which we assumed nutrition in category
  # >70 % are adjusted further below
  valid.cases$AdequacyCalsTot0to30 <- with(valid.cases, {
    (caloriesPercentage <= 30) * NoOralIntakeMat
  })
  # 0 to 30 and oral intake OR 30 to 70 and no oral intake
  valid.cases$AdequacyCalsTot30To70 <-with(valid.cases, {
    ((caloriesPercentage <= 70) * (caloriesPercentage > 30) * NoOralIntakeMat) |
    ((caloriesPercentage <= 30) * OralIntakeMat)
  })
  # 30 to 70 and oral intake OR 70+ OR oral and no EN, no PN
  valid.cases$AdequacyCalsTotAbove70 <-  with(valid.cases, {
    ((caloriesPercentage <= 70) * (caloriesPercentage > 30) * OralIntakeMat)  |
    (caloriesPercentage > 70)
  })
  # Check that on each day only one category
  stopifnot(all(rowSums(with(valid.cases,
    AdequacyCalsTot0to30 + AdequacyCalsTot30To70 + AdequacyCalsTotAbove70))==maxdays.nutri))

  ## Same for sensitivity analysis where last day of ICU stay was imputed
  # (if last day of ICU stay <= 12)
  # The motivation is that most patients didn't spend the whole day at ICU on the
  # day of release, therefore couldn't have had potentially reached higher
  # levels of caloric adequacy.
  # caloriesPercentage2 contains those imputed values (see importHypocaloric.R)
  # The rest is equivalent to the procedures above
  valid.cases$AdequacyCals2Tot0to30 <- with(valid.cases, {
    (caloriesPercentage2 <= 30) * NoOralIntakeMat
  })
  # 0 to 30 and oral intake OR 30 to 70 and no oral intake
  valid.cases$AdequacyCals2Tot30To70 <-with(valid.cases, {
    ((caloriesPercentage2 <= 70) * (caloriesPercentage2 > 30) * NoOralIntakeMat) |
    ((caloriesPercentage2 <= 30) * OralIntakeMat)
  })
  # 30 to 70 and oral intake OR 70+ OR oral and no EN, no PN
  valid.cases$AdequacyCals2TotAbove70 <-  with(valid.cases, {
    ((caloriesPercentage2 <= 70) * (caloriesPercentage2 > 30) * OralIntakeMat)  |
    (caloriesPercentage2 > 70)
  })
  # Check that on each day only one category
  stopifnot(all(rowSums(with(valid.cases,
    AdequacyCals2Tot0to30 + AdequacyCals2Tot30To70 + AdequacyCals2TotAbove70))==maxdays.nutri))

  ## v3: here Oral Intake only improves category of caloric adequacy, when
  # patient was not in MV and received no calories through EN, PN, PF,
  # i.e. caloriesPercentage=0 and (OI and !inMV)=OralIntakeMat
  valid.cases$v3AdequacyCalsTot0to30 <- with(valid.cases, {
    (caloriesPercentage <= 30)  & !((caloriesPercentage==0)*OralIntakeMat)
  })
  # 0 to 30 and oral intake OR 30 to 70 and no oral intake
  valid.cases$v3AdequacyCalsTot30To70 <- with(valid.cases, {
    ((caloriesPercentage <= 70) * (caloriesPercentage > 30)) |
    ((caloriesPercentage == 0) * OralIntakeMat)
  })
  # 30 to 70 and oral intake OR 70+ OR oral and no EN, no PN
  valid.cases$v3AdequacyCalsTotAbove70 <-  with(valid.cases, {
    (caloriesPercentage > 70)
  })
   # Check that on each day only one category
  stopifnot(all(rowSums(with(valid.cases,
    v3AdequacyCalsTot0to30 + v3AdequacyCalsTot30To70 + v3AdequacyCalsTotAbove70))==maxdays.nutri))



  ## if no protocol is available, patient has been released from ICU in which case
  # we assume full nutrition
  # protocol Variable is created in importHypocaloric.R and adjusted in
  # transformToPEM.R
  # protocol == "No" check for days 1-12 if protocol is available, this produces
  # a matrix
  # calendarDaysInICU <= valid.cases
  no.protocol <- valid.cases$protocol=="No"
  valid.cases$AdequacyCalsTot0to30[no.protocol]   <- 0
  valid.cases$AdequacyCalsTot30To70[no.protocol]  <- 0
  valid.cases$AdequacyCalsTotAbove70[no.protocol] <- 1
  stopifnot(all(rowSums(with(valid.cases,
    AdequacyCalsTot0to30 + AdequacyCalsTot30To70 + AdequacyCalsTotAbove70))==maxdays.nutri))


  ## adjust variables for which last day of protocol imputed
  valid.cases$AdequacyCals2Tot0to30[no.protocol]   <- 0
  valid.cases$AdequacyCals2Tot30To70[no.protocol]  <- 0
  valid.cases$AdequacyCals2TotAbove70[no.protocol] <- 1
   stopifnot(all(rowSums(with(valid.cases,
    AdequacyCals2Tot0to30 + AdequacyCals2Tot30To70 + AdequacyCals2TotAbove70))==maxdays.nutri))

  ## adjust v3 version (OI only -> Cat II, else remains)
  valid.cases$v3AdequacyCalsTot0to30[no.protocol]   <- 0
  valid.cases$v3AdequacyCalsTot30To70[no.protocol]  <- 0
  valid.cases$v3AdequacyCalsTotAbove70[no.protocol] <- 1
   stopifnot(all(rowSums(with(valid.cases,
    v3AdequacyCalsTot0to30 + v3AdequacyCalsTot30To70 + v3AdequacyCalsTotAbove70))==maxdays.nutri))



  ############################################################################
  ## Create lag lead matrices

  ## simple lag scheme: 4 day lead, lag until the end
  # should be renamed to LSimple460f
  LSimple460f <- createLmatDyn(lead = 4, lag.c = NULL, n.days = maxdays.nutri,
    brks = brks, time.shift=4)
  valid.cases$LSimple460f <- LSimple460f[seq.id, ]


  ####### hartl's lag scheme  : longer lags if nutrition period was longer ###
  # first interval:  (te+ t.lag)
  # last interval: (te + t.lag) + 4 + te*2
  # --> adjust L accordingly
  # NOTE here we create the matrix for days 1:12, however, in our analyses we
  # only use days 2:12 (as day 1 is first day on ICU and usually not a 24h period)
  # Later we denote this first day on ICU as day 0 and days 2:12 accordingly as
  # days 1:11
  # The first colum of LHartlDyn and other time-dependent covariates is removed
  # below.

  # create L matrix for all intervals
  LHartlDyn <- create_Lmat(te=1:12, te.add.lag=0, t.lag=4, te.add.lead=0,
    lead.const=4, lead.factor=2, interval.seq=brks)
  # subset LHartlDyn based on number of intervals a patient survived
  # see supplementary appendix for description of the Lag/Lead matrices
  valid.cases$LHartlDynf <- Lmat_data(data=valid.cases, Lmat=LHartlDyn)


  ## remove first column from time-dependent covariates specified in vars.ignore
  if(ignore.first) {
    for(i in vars.ignore) {
      valid.cases[[i]] <- valid.cases[[i]][, -1]
    }
  }

  return(valid.cases)

}
