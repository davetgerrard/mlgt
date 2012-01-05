library(roxygen2)

setwd("C:/Users/dave/HalfStarted/mlgt")

package.skeleton(name = "mlgt", code_files="mlgt.R", force=TRUE)

roxygenize("mlgt", unlink.target=T)

