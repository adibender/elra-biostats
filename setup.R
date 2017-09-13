# set maximal follow-up
max.follow <- 30
# for preprocessing, set survival time of patients with survival > max.follow to max.time
max.time <- 30.1
# interval break points (for piece-wise exponential model format)
brks <- c(0:floor(max.time))