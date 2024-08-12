#!/bin/bash

module load HTSlib
module load BCFtools

## autosomal 

iqtree2 -s chromX.min4.phy -m GTR+G -nt 28 -B 1000 -o $outgroup


## X chromosome 

bcftools view -e F_MISSING>0.2 -Oz -r CM020962.1 > chrom_x_clean.vcf.gz

../../softwares/vcf2phylip.py -i chrom_x_clean.vcf.gz -o $outgroup --output-prefix chrom_x

module load IQ-TREE/2.2.2.3-gompi-2022a

iqtree2 -s chromX.min4.phy -m GTR+G -nt 28 -B 1000 -o $outgroup



## Y chromosome 

bcftools view -S males -e F_MISSING>0.2 -Oz -r CM020963.1 > chrom_y_clean.vcf.gz

../../softwares/vcf2phylip.py -i chrom_y_clean.vcf.gz -o $outgroup --output-prefix chrom_y

module load IQ-TREE/2.2.2.3-gompi-2022a

iqtree2 -s chrom_y.min4.phy -m GTR+G -nt 28 -B 1000 -o $outgroup -redo

