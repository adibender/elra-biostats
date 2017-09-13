#### INPUT:
## data.long: data in long format (merge of patient, daily and ICU)
## patient: patient data (1 row, 1 patient)
## daily: daily data (1 patient, number of rows depends on protocolled days)
## ICU: data for the different ICUs (accessed by CombinedicuID)
## timevar: timevar used in reshape operation
## idvar: idvar used in reshape operations
## timevarying: used as v.names in reshape


#### OUTPUT:
## data in wide format

#' Transform data in long format to wide format
#'
#' @importFrom stats reshape
#' @import checkmate dplyr
#' @keywords internal
#' @export
preprocess_for_ped <- function(
  data.long,
  patient,
  ICU,
  timevarying,
  timevar    = "Study_Day",
  idvar      = "CombinedID",
  calSurvVar = "calendarSurvdays") {

  # check inputs
  assert_data_frame(data.long)
  assert_data_frame(patient, all.missing=FALSE)
  assert_data_frame(ICU, all.missing=FALSE)
  assert_character(timevarying, min.chars=1, any.missing=FALSE, min.len=1,
    unique=TRUE)
  assert_subset(timevarying, colnames(data.long))
  assert_string(timevar)
  assert_subset(timevar, colnames(data.long))
  assert_string(idvar)
  assert_subset(idvar, colnames(data.long))
  assert_string(calSurvVar)
  assert_subset(calSurvVar, colnames(data.long))

  ## order by id var, such that rle function works
  data.long <- arrange_(data.long, .dots=idvar)

  ## create wide data to store variables in wide format in
  data.wide <- data.long %>% select_(.dots=idvar) %>% unique %>% as.data.frame

  ## reshape every timevarying variable to wide format,
  # store in data.wide as matrix, b/c mgcv needs one variable which is a
  # matrix for linear functionals
  # need to transform to data.frame, as tbl_df objects would not work here
  data.long <- as.data.frame(data.long)
  for(i in timevarying) {

    tmp <- reshape(data.long,
      timevar = timevar, idvar = idvar, v.names = i,
      direction = "wide",
      drop = names(data.long)[!(names(data.long) %in%
        c(idvar, timevar, i))])

    tmp <- apply(tmp[, !(colnames(tmp) %in% idvar)], 2, unclass)
    dimnames(tmp) <- NULL
    data.wide[[i]] <- tmp

  }

  ## merge patient data with dataWide (both 1 row -> one patient)
  datap.wide <- left_join(patient, data.wide, by = idvar)

  ## compute number of days with nutrition protocol/length ICU stay
  study.days <- do.call(cbind.data.frame, rle(data.long[, idvar]))
  colnames(study.days) <- c("dayMax", idvar)
  maxdays.nutri <- max(study.days$dayMax)

  ## merge patient/daily data with study days (basically add new column,
  # but saver in case of different orderings)
  datap.wide <- merge(datap.wide, study.days, by = idvar)

  cal.surv.var <- datap.wide[[calSurvVar]]

  ## add info about incomplete protocols
  datap.wide$incompleteProtocol <- with(datap.wide,
    (DiedBeforeOrAtMax & (dayMax < cal.surv.var) & (calendarDaysInICU <= dayMax)) |
    (!DiedBeforeOrAtMax & (dayMax < maxdays.nutri) & (calendarDaysInICU <= dayMax)))

  id.died.inc <- with(datap.wide, DiedBeforeOrAtMax & (dayMax < cal.surv.var))
  id.censored.inc <- with(datap.wide, !DiedBeforeOrAtMax & (dayMax < maxdays.nutri))
  # adjust protocol variable (dayMax + 1, b/c dayMax itself has information)
  datap.wide$protocol <- matrix("Yes", nrow = nrow(datap.wide), ncol = maxdays.nutri)
  for(i in which(id.died.inc)) {
    # print(i)
    datap.wide$protocol[i,
    (datap.wide$dayMax[i] + 1):datap.wide[[calSurvVar]][i]] <- "No"
  }# NOTE: b/c intilized with "Yes", its now "Yes" for days beyond survival
   # doesn't matter b/c of Lag-Lead-Mat set to 0 later
  for(i in which(id.censored.inc)) {
    datap.wide$protocol[i, (datap.wide$dayMax[i] + 1):maxdays.nutri] <- "No"
  }

  ## merge patient/daily data with ICU data
  data.all.wide <- merge(datap.wide, ICU, by = c("CombinedicuID", "Year"))

  ## add matrix to patient/daily/icu data (still in wide format)
  # basically a 1:11 vector per row
  data.all.wide$DaysMat <- matrix(1:maxdays.nutri, byrow = TRUE,
    ncol = maxdays.nutri, nrow = nrow(data.all.wide))

  ## order by ID and return
  data.all.wide[order(data.all.wide[, idvar]), ]

}


#' transform data in wide format back to long/ped format
#'
#' @inherit preprocess_for_ped
#' @keywords internal
#' @export
wide_to_ped <- function(
  data.wide,
  timevarying,
  brks    = c(0:12, seq(15, 55, by = 5), 61.1),
  idvar   = "CombinedID",
  timevar = "Study_Day",
  survvar = "Survdays",
  max.follow = 30) {

  assert_data_frame(data.wide)
  assert_character(timevarying, min.chars=1, any.missing=FALSE, min.len=1,
    unique=TRUE)
  assert_numeric(brks, lower=0, finite=TRUE, any.missing=FALSE, min.len=2,
    unique=TRUE)
  assert_subset(timevarying, colnames(data.wide))
  assert_string(timevar)
  assert_subset(timevar, colnames(data.long))
  assert_string(idvar)
  assert_subset(idvar, colnames(data.long))


  ## create intervals factor,
  int.names <- paste0("(", paste0(brks[-length(brks)], ",", brks[-1]), "]")
  intervals <- factor(int.names, levels = int.names)
  ##compute interval lengths
  diffs <- diff(brks)

  maxdays.nutri <- ncol(data.wide$DaysMat)
  ## length of mechvent
  mech.diff <- data.wide$endMechVent - data.wide$beginMechVent
  calendar.mechdiff <- with(data.wide, calendarEndMechVent - calendarBeginMechVent)

  ## number of intervals survived
  # NOTE: don't use calendar survival days here, b/c we interested in actual survival
  # not nutrition-/calendar-days
  ints.vec <- pmin(pmax(colSums(sapply(data.wide[[survvar]], ">", y = brks)), 1), max.follow)
  use.vec  <- pmin(data.wide$dayMax, ints.vec)

  ## index for which days nutrition protocol available
  use.vec <- pmin(use.vec, maxdays.nutri)

  ## prep functional covariates:
  T12 <- rep(TRUE, maxdays.nutri)
  not.use <- lapply(use.vec, function(z) {
    T12[seq_len(z)] <- FALSE
    T12
  })
  # matrix with one row for each patient, TRUE indicates that no nutrition
  # protocol was available at this day of nutrition protocol
  not.use.mat <- do.call(rbind, not.use)

  ## set all values after dayMax/death to 0 (NA trouble otherwise)
  ## define integration weights L (0 after dayMax)
  ## set unused values in L to 0
  df.list <- lapply(seq_along(data.wide[[idvar]]),  function(z) {

    d <- data.wide[z, ]
    surv <- min(d[[survvar]], max.follow)
    # what's the last interval border the patient survives
    ints <- ints.vec[z]

    # offset in interval k := max(0, min(length(interval k),
    # time- left side of interval k)
    offset <- if ( ints == 1 ) {
      log(surv)
    } else {
      tmp <- diff(brks[1:(ints+1)])
      tmp[ints] <- surv - brks[ints]
      log(tmp)
    }
    # It's OK to use Patient Died here, b/c PatientDied was set to zero
    # in importHypocaloric.R when death occured after maximal follow up time
    event    <- c(rep(0, ints-1), d$PatientDied)
    survtime <- cumsum(exp(offset))
    int      <- intervals[1:length(survtime)]
    intmid   <- (brks[2:(ints + 1)] + brks[1:ints]) / 2

    ## make time varying covariate: proportion of each interval
    # spent in ICU/MechVent
    ## edit: round because proportions make no sense if intervals
    # don't have same length
    inICU <- round(pmin(1,
      pmax(0, (d$beginHosp - brks[1:ints]) / diffs[1:(ints)])))

    #inMechVent==1 if there was _any_ intubation in that interval
    inMechVent <- (d$beginMechVent <= brks[2:(ints+1)]) *
      (d$endMechVent >= brks[1:ints]) *
      ceiling(pmin(1, mech.diff[z]/diffs[1:(ints)]))

    ## define matrix covariate for current interval:
    intMat <- matrix(intmid, nrow = ints, ncol = maxdays.nutri)

    # setup lag/lead matrix
    L <- matrix(1, nrow=ints, ncol=maxdays.nutri) *
      lower.tri(matrix(1, nrow=ints, ncol=maxdays.nutri), diag=TRUE) *
      matrix(!not.use.mat[z,], nrow=ints, ncol=maxdays.nutri, byrow=TRUE)

    ## create new df
    list(
      offset     = offset,     event  = event,  survtime = survtime,
      int        = int,        intmid = intmid, inICU    = inICU,
      inMechVent = inMechVent, intMat = intMat, L=L)

  })

  for(j in timevarying) {
    if(is.numeric(data.wide[[j]])) {
      data.wide[[j]][not.use.mat] <- 0
    }
    else{
      if( j %in% c("EN", "PN", "OralIntake", "Propofol", "protocol")){
        data.wide[[j]][not.use.mat] <- "No"
      }
      else{
        if( j == "proteinCat" ) {
          data.wide[[j]][not.use.mat] <- "lower"
        }
        else stop("New unnknown variable detected, adjust NA replacement.")
      }
    }
  }

  names.list   <- names(df.list[[1]])
  classes.list <- sapply(df.list[[1]], class)
  mats         <- classes.list %in% c("matrix")
  vecs.list    <- list()

  for( j in names.list[!mats] ) {
    vecs.list[[j]] <- unlist(lapply(df.list, "[[", i = j))
  }

  data.pem <- do.call(cbind.data.frame, vecs.list)

  for(j in names.list[mats]) {
    data.pem[[j]] <- do.call(rbind, lapply(df.list, "[[", i = j))
  }

  n.obs <- sapply(lapply(df.list, "[[", i = "offset"), length)
  data.pem <- cbind(CombinedID = rep(data.wide[[idvar]], times = n.obs),
    data.pem)
  data.wide <- data.wide[rep(1:nrow(data.wide), times = n.obs), ]

  if(!all(data.pem[[idvar]] == data.wide[[idvar]]))
    stop(paste0(idvar, "s not equal"))

  data.pem <- cbind(data.pem, data.wide)
  # NOTE: rounding b/c exp/log operations make times like 30 uneven, 30 -> 29.999999..
  # Gives problems later when calculating C-index
  # Precision beyond 5th decimal not neeeded anyway
  data.pem$survtime <- round(data.pem$survtime, 5)
  data.pem[!duplicated(colnames(data.pem))]

}



#' transforms data in long format into data in ped format (Piece-wise
#' exponential data)
#'
#' @inherit preprocess_for_ped
#' @importFrom plyr join
#' @keywords internal
#' @export
transform_to_ped <- function(
  data.long,
  timevarying,
  patient,
  ICU,
  brks,
  idvar      = "CombinedID",
  timevar    = "Study_Day",
  min.day    = 0,
  survvar    = "Survdays",
  calSurvVar = "calendarSurvdays",
  max.follow = 30) {


  assert_data_frame(data.long)
  assert_character(timevarying, min.chars=1, any.missing=FALSE, min.len=1,
    unique=TRUE)
  assert_data_frame(patient)
  assert_data_frame(ICU)
  assert_numeric(brks, lower=0, finite=TRUE, any.missing=FALSE, min.len=2,
    unique=TRUE)
  assert_string(timevar)
  assert_subset(timevar, colnames(data.long))
  assert_string(idvar)
  assert_subset(idvar, colnames(data.long))
  assert_string(calSurvVar)
  assert_subset(calSurvVar, colnames(data.long))
  assert_numeric(min.day, lower=0, finite=TRUE, any.missing=FALSE, min.len=1)


  ## transform timevarying data to wide format and merge with patient and ICU
  data.wide <- preprocess_for_ped(data.long, patient, ICU, timevarying, timevar,
    idvar, calSurvVar=calSurvVar)

  ## transform timevarying variables to pem format and adjust rest of data
  wd <- wide_to_ped(data.wide, timevarying, brks, idvar, survvar=survvar,
    max.follow=max.follow)

  int.info <- int_info(brks = brks, min.int = min.day)

  # need plyr here, as dplyr::left_join can't handle matrix columns
  join(wd, int.info, by = "intmid", type="left")

}
