## uses Template to create submit model files for
# 1) data where we assume that patients released from hospital survival until t=30
# 2) data where we consider release from hospital as censoring
fit.function <- c("gam")
data.sets    <- c("Current", "CurrentHosp")

registry.setups <- expand.grid(fit.function, data.sets)

colnames(registry.setups) <- c("function", "data")

base.formula <- 'event ~ s(intmid, bs = "ps") +
  ApacheIIScore + ApacheIIScore:intmid +
  s(Age, by = intmid, bs = "ps") +
  s(BMI, by = intmid, bs = "ps") +
  Year + DiagID2 + AdmCatID + Gender +
  freqInMVd2to4 + freqPFd2to4 +  freqOId2to4 + freqPNd2to4 +
  s(CombinedicuID, bs = "re", by=icuByDummy)'

for( i in seq_len(nrow(registry.setups)) ) {

  temp <- readLines("runModelsTemplate.R")

  temp[grep("^XXbaseformulaXX*", temp)] <- paste0("base.formula <- '",
    base.formula, "'")
  temp[grep("^XXdatapathXX*", temp)] <- paste0("data.folder <- '",
    paste0("data", registry.setups[i, "data"]), "'")
  suffix.data <- ifelse(registry.setups[i, "data"] == "Current", "", "Hosp")
  regname <- paste0(registry.setups[i, "function"], suffix.data)
  temp[grep("^XXregnameXX*", temp)] <- paste0("reg.name <- '", regname, "Batch'")
  temp[grep("^XXrunmodel.functionXX*", temp)] <- paste0("function.name <- '",
    ifelse(registry.setups[i, "function"] == "gam", "runOneModelGAM.R",
      "runOneModel.R"), "'")

  writeLines(temp, paste0("submitModels", toupper(registry.setups[i, "function"]),
    suffix.data, ".R"))

}
