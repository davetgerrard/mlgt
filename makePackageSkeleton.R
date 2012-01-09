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


setwd("C:/Users/dave/HalfStarted/mlgt/testREADME")
Sweave("../README")	# from sub-directory up to main mlgt directory
R CMD texify --pdf README.tex


