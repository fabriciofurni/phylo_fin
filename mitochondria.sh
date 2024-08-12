#!/bin/bash

## getting mitochondrial genomes from BAMs 

module load angsd/0.940-GCC-11.2.0

for sample in $(cat sampleinfo | sort | uniq)
do

angsd -doCounts 1 -doFasta 2 -out ${outputfolder}/mitochondrial/${sample} -i /scratch/$USER/genome/${sample}/${sample}_dedup_reads.bam -r NC_001601.1:

## Aligning with MAFFT 

## To include published sequences include all sequences in the "${outputfolder}/mitochondrial/" directory before 'cat' 

cat ${outputfolder}/mitochondrial/*.fa.gz > ${outputfolder}/mitochondrial/all_mito_sequences.fa

module load MAFFT

mafft --auto --reorder "${outputfolder}/mitochondrial/all_mito_sequences.fa" > "${outputfolder}/mitochondrial/all_final_aligned.fa"

### Here you should manually check the alignments and trim the ends. 

## IQtree

module load IQ-TREE/2.2.2.3-gompi-2022a

iqtree2 -s all_final_aligned.fa -o ${outgroup} -m GTR+G -B 1000 -redo


