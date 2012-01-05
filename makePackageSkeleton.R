library(roxygen2)

#setwd("C:/Users/dave/HalfStarted/mlgt")
setwd("C:/Users/dave/HalfStarted/mlgt/testRoxygen")

package.skeleton(name = "mlgt", code_files="mlgt.R", force=TRUE)

roxygenize("mlgt", unlink.target=T)




### to use sweave on documentation
# in R
Sweave("testSweave.R")
# then in DOS.
R CMD texify --pdf testSweave.R.tex


