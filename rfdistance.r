#################### testing RF distance #######################################

library(ape)
library(TreeDist)

##### creating list of tree files from the 1Mbp Loci trees ##

listcsv <- dir(pattern = "*.raxml.support")

listcsv <- read.csv("list_files.txt")

write.csv(listcsv,file= "list_files.txt")

## openning files as trees using lapply and ape  ##

t <- read.csv("list_files.txt", header = FALSE,sep = ",")

list <- lapply(t$V1,read.tree)


## load whole autosome tree ###

nucl <- read.tree("../whole_autos_raxml.support")

## load mito consensus tree ##

mito <- read.tree("../mitotree.raxml.support")

######## apply function to compare with whole autosomal #########

sapply(list,function(z) SharedPhylogeneticInfo(nucl,z, normalize = TRUE))

######## apply function to compare with mitocondiral ########

sapply(list,function(z) SharedPhylogeneticInfo(mito,z, normalize = TRUE))


