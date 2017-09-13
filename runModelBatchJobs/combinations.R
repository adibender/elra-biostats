## create a table of combinations of penalty and lag for which models should be
# fit (we run 3 models here, first model is the application example,
# models 2 and 3 are needed b/c their coefficients will be used in simulations)
# Model 1 will also be rerun with different censoring assumptions for sensitivity
# analysis
data.names  <- c("full")
lag.schemes <- c("Expert", "L460")
penalties   <- c("t1", "t2")
nutrition   <- c("Calories")

combinations <- expand.grid(
	data      = data.names,
	lag       = lag.schemes,
	penalty   = penalties,
	nutrition = nutrition)[-4, ]

# save as matrix with row names
combinations <- as.matrix(combinations)
rownames(combinations) <- apply(combinations, 1, paste0, collapse = "")

saveRDS(combinations, file = "combinations.Rds")