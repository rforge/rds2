# 
# Author: johnros
###############################################################################


package.skeleton(
		name='rds2', 
		code_files = "/home/johnros/workspace/rds2/pkg/build/onlyFunctions.R", 
 		path='/home/johnros/workspace/rds2/pkg', 
		force = TRUE)


# Roxygen was used to create the intial package. Now the files should be edited manualy.
require(plyr)
require(roxygen2)
roxygenize(
 		package.dir = "/home/johnros/workspace/rds2/pkg/rds2/",
 		roxygen.dir = "/home/johnros/workspace/rds2/pkg/rds2/",
 		copy.package = FALSE,
 		unlink.target = TRUE, 
 		overwrite = TRUE)
