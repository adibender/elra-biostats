## Data import and basic preprocessing

The initial data import and basic preprocessing steps are not included 
in this respository because they contain sensitive infomation on the patient 
level.


Most important covariates are Date covariates giving information on 
	
  - `AdmissionDate`: Time of Admission to ICU (Time 0 of the survival process)
  - `icuDischargeDate`: Time of ICU discharge (used as hospital discharge when 
  hospital discharge is missing)
  - `PatientDeathTime`: Time of Death (Used to calculate survival time for non-censored as 
  `PatientDeathTime` - `AdmissionDate`)
  - `HospDischargeDate`: Time of Hospital discharge (used as censoring covariate
  in sensitivity analysis)

The data sets provided in the **`data`** folder are the preprocessed data sets, 
where date times are excluded and differences of the date times included instead.

## Hospital discharge vs. death

In our main analysis we consider survival time as time from *Admission to ICU* 
until *patient death time*, for patients that do not die until the end of follow 
up (60days max, 30 days for acute survival), we assume that they survived until 
the end of the follow up (i.e. 30 days). 

In the sensitivity analysis we consider patients without an event as beeing 
censored at time of hospital discharge. 

As the time of censoring (end of follow up vs. at hospital discharge) changes 
the data sets in piece-wise exponential format significantly (one is not simply 
a *row*-subset of the other two sets of data sets will be created and saved, 
see `../dataCurrent/` (for data where right censoring only occurs at the end 
of the follow up) and `../dataCurrentHosp/` where hospital discharge was 
considered as censoring time)


## Hospital admission/discharge
- 591 patients with missing icu discharge date (thus missing begin hospital)
- 5150 patients with missing hospital discharge date 
- 4559 of those have non missing icu-discharge date and we set hosp discharge to 
icu discharge date (2414 of those are equivalent to patient death time)
- 591 patient remain NA (same 591 without begin hosp), as these patients also 
do not have a death time, we assume that they remain in icu/hospital until the 
end of the follow up.
- variable `survhosp` is created from the difference between admission to 
icu (t=0) and hospital discharge date
- when creating the data in `../dataCurrent` we use `Survdays` as variable 
containing survival times
- when creating the data in `../dataCurrentHosp` we use `survhosp` as variable
containing survival times
