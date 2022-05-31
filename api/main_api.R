library(plumber)
r <- plumb("/cnvscore/plumber.R")
r$run(host = "0.0.0.0", port=3838, swagger = TRUE)