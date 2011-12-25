# TODO: Add comment
# 
# Author: johnros
###############################################################################


require(roxygen2)

package.skeleton(name='rds2', code_files = "/home/johnros/workspace/rds2/onlyFunctions.R", 
 		path='/home/johnros/workspace/rds2/pkg', force = TRUE)

roxygenize(
 		package.dir = "/home/johnros/workspace/rds2/pkg/rds2/",
 		roxygen.dir = "/home/johnros/workspace/rds2/pkg/rds2/",
 		copy.package = FALSE,
 		unlink.target = TRUE, 
 		overwrite = TRUE)
