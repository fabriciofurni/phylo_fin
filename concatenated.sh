#!/bin/bash

module load HTSlib
module load BCFtools


## autosomal 

vcf2phylip.py -i ${outputfolder}/vcfs/${proj}_snps_auto_filtered.vcf.gz -o $outgroup --output-prefix auto

module load IQ-TREE/2.2.2.3-gompi-2022a

iqtree2 -s auto.min4.phy -m GTR+G -nt 28 -B 100 -o $outgroup


## X chromosome 

vcf2phylip.py -i ${outputfolder}/vcfs/${proj}_X_chrom_snps_filtered.vcf.gz -o $outgroup --output-prefix chrom_x

module load IQ-TREE/2.2.2.3-gompi-2022a

iqtree2 -s chromX.min4.phy -m GTR+G -nt 28 -B 1000 -o $outgroup


## Y chromosome 

vcf2phylip.py -i  ${outputfolder}/vcfs/${proj}_Y_chrom_snps_filtered.vcf.gz -o $outgroup --output-prefix chrom_y

module load IQ-TREE/2.2.2.3-gompi-2022a

iqtree2 -s chrom_y.min4.phy -m GTR+G -nt 28 -B 1000 -o $outgroup -redo

