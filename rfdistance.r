#################### testing RF distance #######################################

library(ape)
library(TreeDist)



## load whole autosome tree ###

nucl <- read.tree("../whole_autos_raxml.support")

## load mito consensus tree ##

mito <- read.tree("../mitotree.raxml.support")

######## apply function to compare with whole autosomal #########

sapply(list,function(z) SharedPhylogeneticInfo(nucl,z, normalize = TRUE))

######## apply function to compare with mitocondiral ########

sapply(list,function(z) SharedPhylogeneticInfo(mito,z, normalize = TRUE))


