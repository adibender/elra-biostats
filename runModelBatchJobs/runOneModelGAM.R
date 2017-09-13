runModel <- function(
  combination, 
  all.combinations, 
  base.formula, 
  nutri.formulae,
  data.files,
  save.path = "") {

  library(mgcv)
  comb.z <- all.combinations[combination, ]
  data.z        <- readRDS(data.files[[comb.z["data"]]])
  nutri.formula <- nutri.formulae[[comb.z["penalty"]]]

  if( comb.z["lag"] == "L460") {
    nutri.formula <- gsub("LHartlDynf", "LSimple460f", nutri.formula)
  }
  if( comb.z["nutrition"] == "None") {
    nutri.formula <- ""
  }

  # modify formula if variables drop due to no variation (e.g. for subgroups)
  formula.string <- paste0(base.formula, nutri.formula)
  formula.z      <- as.formula(formula.string)
  if( comb.z["data"] == "medical") {
    formula.z <- update(formula.z, .~. - AdmCatID)
  }
  if( comb.z["data"] %in% 
    c("noOId4NotMedical") ) {
    formula.z <- update(formula.z, .~. - freqOId2to4)
  }

  save.string <- paste0(comb.z, collapse = "")

  m.z <- gam(formula.z, 
    data    = data.z,
    offset  = offset, family    = poisson(),
    control = gam.control(trace = TRUE),
    method  = "REML")


  ## save model
  saveRDS(m.z, file = paste0(save.path, save.string, ".Rds"))

  ## session info
  print(sessionInfo()) # will be stored in the .out files of the jobs folders of the registry

  invisible()

}
