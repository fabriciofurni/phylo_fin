###### Phylogeographical structure test #####################################

NOPSAMPLES=$definevari
NATSAMPLES=$definevari
SOUSAMPLES=$definevari

### packages 

library(ape)
library(TreeDist)

##### creating list of tree files from the 1Mbp Loci trees ##

listcsv <- dir(pattern = "*.raxml.support")

listcsv <- read.csv("list_files.txt")

write.csv(listcsv,file= "list_files.txt")

## openning files as trees using lapply and ape  ##

t <- read.csv("list_files.txt", header = FALSE,sep = ",")

list <- lapply(t$V1,read.tree)

## testing monophyly by geographical localation for each tree file ##

NOP <- sapply(list,function(z) is.monophyletic(z,c("$NOPSAMPLES")))

write.csv(NOP,file= "NOP.csv")

NAT <- sapply(list,function(z) is.monophyletic(z,c("$NATSAMPLES")))

write.csv(NAT,file= "NAT.csv")

SOU <- sapply(list,function(z) is.monophyletic(z,c("$SOUSAMPLES")))

write.csv(SOU,file= "SOU.csv")

