<div id="header">

</div>

## Scripts used on Furni et. al. (2024). Phylogenomics and pervasive genome-wide phylogenetic discordance among fin whales (_Balaenoptera physalus_)   


#### [alignment_calling.sh](https://github.com/fabriciofurni/phylo_fin/blob/main/alignment_calling.sh) <a name="alignment_calling"></a>


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

#### [filtering_snps.sh](https://github.com/fabriciofurni/phylo_fin/blob/main/filtering_snps.sh) <a name="filtering_snps"></a>

Script to create the final SNP dataset. The script does joint-genotyping samples and filters by quality, number of alleles, and others

- Combining samples GVCF files by chromossome with GATK CombineGVCF
- Genotype all samples by chromossome with GATK GenotypeGVCFs
- Annotate variants functional impact by chromossome using SNPEFF 
- Merging chromossomes VCF files with PICARDTOOL GatherVcfs
- Selecting only SNPs with BCFTOOLS filtering
- Filter SNPs only for the autosomal, and sex-chromosomes genomes


#### [concat.sh](https://github.com/fabriciofurni/phylo_fin/blob/main/concatenated.sh) <a name="concatenated"></a>

Script to run concatenated autosomal, sex-chromosome phylogenetic estimates 

- Converting files from VCF to Phyllip. 
- Run phylogenies
  

#### [mitochondria.sh](https://github.com/fabriciofurni/phylo_fin/blob/main/mitochondria.sh) <a name="mitochondria"></a>

Script to obtain mitochondrial sequences and run alignment with MAFFT. Code for run in IQTree after manual curation. 

- Get mitochondrial sequences from BAM files using ANGSD
- Align sequences with MAFFT
- Run IQtree phylogenetic estimates 


#### [windows_phylo_d.sh](https://github.com/fabriciofurni/phylo_fin/blob/main/windows_phylo_d.sh) <a name="windows"></a>

Script for analyses in genomic windows.

- Create input files from the whole autosomal data based on chromosome coordinates
- Run conversion from vcf to phyllip
- Run IQtree phylogenetic estimates
- Run Dsuite estimates

#### [astral_gcf.sh](https://github.com/fabriciofurni/phylo_fin/blob/main/astral_gcf.sh) <a name="astral_gcf"></a>

Script to run ASTRAL and GCF concordance checks per windows

- Run Astral per window
- Run GCF distant per window
  



 
