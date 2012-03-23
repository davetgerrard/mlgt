library(roxygen2)

#setwd("C:/Users/dave/HalfStarted/mlgt")
setwd("C:/Users/dave/HalfStarted/mlgt/testRoxygen")

# Set the correct dates in mlgt.R and DESCRIPTION
# NO # package.skeleton(name = "mlgt", code_files="mlgt.R", force=TRUE)
# Update DESCRIPTION FILE for VERSION
# delete contents of /man 
# copy latest .R code to /R

roxygenize("mlgt", unlink.target=T)

#edit the namespace. Remove "import(seqinr)"

# In DOS, from same directory (/testRoxygen)
#R CMD INSTALL --build mlgt	# to get windows zip
#R CMD build mlgt	# to get tar.gz
# R CMD check mlgt




## Install the new package version
###Generate the README 

setwd("C:/Users/dave/HalfStarted/mlgt/testREADME")
#Sweave("../mlgt_README")	# from sub-directory up to main mlgt directory
#R CMD texify --pdf mlgt_README.tex
#Stangle("../mlgt_README")
Sweave("../mlgt_README_0_15")	# from sub-directory up to main mlgt directory
#R CMD texify --pdf mlgt_README_0.15.tex
Stangle("../mlgt_README_0_15")



