<div id="header">

</div>

## Scripts used on Furni et. al. (2023). Pervasive, genome-wide phylogeographic disruptions among fin whale (_Balaenoptera physalus_) subspecies


## Alignment and Variant Calling 

#### Raw reads alignment steps 

- Alignment of raw reads with BWA mem per library
- Conversion from SAM to BAM format with SAMTOOLS
- Merging BAM files from different libraries per sample 

#### Marking and removing duplicates 

- Mark and removing duplicates using PICARDTOOLS MarkDuplicates 
- create BAM index with PICARDTOOLS BuildBamIndex
- get alignment metrics with PICARDTOOLS CollectAlignmentSummaryMetrics
- get insert metrics with PICARDTOOLS CollectInsertSizeMetrics
- get overall metrics with PICARDTOOLS CollectWgsMetrics
- get quality metrics with PICARDTOOLS QualityScoreDistribution
- get per chromossome coverage and depth with SAMTOOLS coverage

#### Rescaling aligned reads 

- Call variants per chromossome using GATK HaplotypeCaller 
- Combine chromossome GVCF files with GATK CombineGVCF
- Genotype GVCF file with GATK GenotypeGVCF 
- Select SNPS or INDELS with GATK SelectVariants 
- Filter for high confidence variants with GATK VariantFiltration
- Apply high confidence variants with GATK BaseRecalibrator
- Recalibrating BAM files with GATK ApplyBQSR
- Get recalibration metrics with GATK AnalyzeCovariates

#### Calling Genotypes per sample by chromossome 

- Call variants using the rescaled BAM files with GATK HaplotypeCaller in ERC mode per chromossome
- Combine chromossome GVCF files using GATK CombineGVCF per sample 


## Downstream analyses 

#### Script 1

Script to create the final SNP dataset. The script does joint-genotyping samples and filters by quality, number of alleles, and others

- Combining samples GVCF files by chromossome with GATK CombineGVCF
- Genotype all samples by chromossome with GATK GenotypeGVCFs
- Annotate variants functional impact by chromossome using SNPEFF 
- Merging chromossomes VCF files with PICARDTOOL GatherVcfs
- Selecting only SNPs with BCFTOOLS filtering
- Filter SNPs only for the autosomal, and sex-chromosomes genomes


#### Script 2

Script to run whole autosomal, sex-chromosome phylogenetic estimates 

- Converting files from VCF to Phyllip using . 
- Run phylogenies using RAXML
  

#### Script 3 

Script to obtain mitochondrial sequences and run alignment with MAFFT. Code for run in RAXML after manual curation. 

- Get mitochondrial sequences from BAM files using ANGSD
- Align sequences with MAFFFT
- Run IQtree phylogenetic estimates 


## Script 4 

Script for analyses in genomic windows.

- Create input files from the whole autosomal data based on chromosome coordinates
- Run conversion from vcf to phyllip
- Run IQtree phylogenetic estimates
- Run Dsuite estimates

## Script 5

Script to test for phylogeographical structure based on monophyly and RF distances.

- Load libraries
- Input phylogeny results from the windows
- Test for monophyly using ape
- Collecting results
- Test RF distance for both mitochondrial and nuclear topologies

  



- 


 
