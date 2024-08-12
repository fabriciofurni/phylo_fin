#!/bin/bash

for window in 1mb 100k 50k

do

## ASTRAL consensus trees 

## getting trees only 

cat ${window}_Consensus_BT_trees.txt | cut -f1 > ${window}_trees.txt

astral-hybrid -i ${window}_trees.txt -o ${window}_consensus_astral


## Checking trees concordance factor within the windows 

## consensus tree

iqtree2 --gcf ${window}_trees.txt -t ${window}_consensus_astral

## concatenated autosomal tree 

iqtree2 --gcf ${window}_trees.txt -t ${outputfolder}/auto/${window}_auto.min4.phy.treefile 

### mitochondrial mitochondrial tree

iqtree2 --gcf ${window}_trees.txt -t ${outputfolder}/mitochondrial/${window}_all_final_aligned.treefile 
