## for pattern_label function
# define number of days time-dependent covariate was provided
maxdays.nutri <- 11
# define protocols for comparisons
protocols.df <- data.frame(
  low     = rep("low", maxdays.nutri),
  mid     = rep("mid", maxdays.nutri),
  full    = rep("full", maxdays.nutri),
  lowmid  = c(rep("low", 4), rep("mid", maxdays.nutri-4)),
  midfull = c(rep("mid", 4), rep("full", maxdays.nutri-4)))
# define which protocols should be compared pairwise
protocols.to.compare <- data.frame(
  pat1 = c("low"    , "low" , "lowmid" , "low"  , "mid" , "mid"),
  pat2 = c("lowmid" , "mid" , "mid"    , "full" , "midfull"    , "full"))

protocols.to.compare <- cbind(protocols.to.compare,
  comparison=paste("Comparison", LETTERS[seq_len(nrow(protocols.to.compare))]))
protocols.to.compare$label <-diff.labs <- apply(protocols.to.compare, 1,
        function(z) {
            paste(
                pattern_label(protocols.df[, z[1]], maxdays.tdc = maxdays.nutri),
                "vs. \n",
                pattern_label(protocols.df[, z[2]], maxdays.tdc = maxdays.nutri))
        })

saveRDS(protocols.to.compare, "protocolsToCompare.Rds")
saveRDS(protocols.df, "protocolsDF.Rds")
