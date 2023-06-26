#!/bin/bash

## Some steps were modified from https://gencore.bio.nyu.edu/variant-calling-pipeline/
## by Fabricio Furni 

# This is a shell automated script that takes as input raw fastq files and output reads aligment, called variants, and metric files for whole-genome analyses. 
# The script uses BWA mem, SAMTOOLS, PicardTools, and the GATK best practices to produce BAM and VCF files per individual. Outputs include .BAM .GVCF  
# It creates a set of bash scripts for SLURM jobs which automatically starts one after another through job dependency. 


help_menu () {
	echo ""
	echo ""
	echo "Required arguments:"
	echo ""
	echo "-g <path to bwa indexed reference genome>"
	echo "-l <path to project folder containing sample's folders with pair-ended fastq files. e.g /data/project >"
	echo "-o <output folder without ending /. e.g. data/output >"
	echo "-c <path to chromossomes list file>"
	echo ""
	echo "-h <this helpfull message>"
	echo ""
}

arg_check () {
if [ -z "${genome}" ]
then
	echo ""
	echo "no value provided for -g"
	help_menu
	exit 1
elif [ -z "${library}" ]
then
	echo ""
	echo "no value provided for -l"
	help_menu
	exit 1
elif [ -z "${outputfolder}" ]
then
	echo ""
	echo "no value provided for -o"
	help_menu
	exit 1
elif [ -z "${chromofile}" ]
then
	echo ""
	echo "no value provided for -c"
	help_menu
	exit 1
fi
}

while getopts :g:l:o:c:h FLAG; do
	case "${FLAG}" in
		g) 	genome=${OPTARG}
			;;
		l) 	library=${OPTARG}
			;;
		o)	outputfolder=${OPTARG}
			;;
		c)	chromofile=${OPTARG}
			;;
		h) 	help_menu
			;;
		\?) echo "Invalid option: ${OPTARG}" 1>&2
			;;
		:) 	echo "Invalid option: ${OPTARG} requires an argument" 1>&2
			;;
	esac
done

shift $((OPTIND -1))

arg_check

pwd=$(pwd)

ls -p ${library} | grep / | sed 's/[/]//g' > sampleinfo

## creating output folders ##

mkdir /scratch/$USER/genome/
mkdir /scratch/$USER/genome/gvcfs/

mkdir ${outputfolder}/bams/
mkdir ${outputfolder}/fastas/
mkdir ${outputfolder}/gvcfs/
mkdir ${outputfolder}/snps/
mkdir ${outputfolder}/indels/
mkdir ${outputfolder}/metrics/
mkdir ${outputfolder}/metrics/alignment
mkdir ${outputfolder}/metrics/duplicates
mkdir ${outputfolder}/metrics/snps

###### creating automated bash scripts for each sample in the project folder #######

for sample in $(cat sampleinfo | sort | uniq)
do


mkdir ${library}/${sample}/sbatches
mkdir /scratch/$USER/genome/${sample}
mkdir /scratch/$USER/genome/${sample}/gvcfs/
mkdir /scratch/$USER/genome/${sample}/vcfs/
mkdir ${outputfolder}/${sample}
mkdir ${outputfolder}/${sample}/metrics

cd ${library}/${sample}/

## collect information from the fastq name ##

ls *1.fq.gz | sort | uniq > ${library}/${sample}/c1
ls *2.fq.gz | sort | uniq > ${library}/${sample}/c2

paste ${library}/${sample}/c1 ${library}/${sample}/c2 > ${library}/${sample}/libraries

## using readgroup library info to create alignment scripts with respective BAM HEADER ## 

cat ${library}/${sample}/libraries | while read line;

do

R1=$( echo "${line}" | basename $(cut -f1))
R2=$( echo "${line}" | basename $(cut -f2))
ID=$( echo "${line}" | basename $(cut -f1) | sed -e 's/_1.fq.gz//')
LB=$(echo "${line}" | basename $(cut -f1) | cut -d "_" -f3)
PU=$(echo "${line}" | basename $(cut -f1) | cut -d "_" -f1,2)




##################################### Script 1: Genome Alignment with BWA ################################################################################



set +H

echo "#!/bin/bash

module load BWA/0.7.17-foss-2018a

bwa mem -t 8 -M -R '@RG\tID:$ID\tLB:$LB\tPL:DNBSEQ\tPU:$PU\tSM:${sample}' ${genome} ${library}/${sample}/${R1} ${library}/${sample}/${R2} > /scratch/$USER/genome/${sample}/${ID}_aligned_reads.sam

module load SAMtools/0.1.20-foss-2018a

samtools view -bS -@ 8 -o /scratch/$USER/genome/${sample}/${ID}_${ref}_aligned_reads.bam /scratch/$USER/genome/${sample}/${ID}_aligned_reads.sam

rm /scratch/$USER/genome/${sample}/${ID}_aligned_reads.sam

module purge

### Collecting first Alignment Metrics ####

module load picard/2.18.17-Java-1.8 R/3.6.1-foss-2018a

java -jar \$EBROOTPICARD/picard.jar SortSam \
I=/scratch/$USER/genome/${sample}/${ID}_${ref}_aligned_reads.bam \
O=/scratch/$USER/genome/${sample}/${ID}_${ref}_sorted_reads.bam \
SO=coordinate \
TMP_DIR=/scratch/$USER/tempfiles/

rm /scratch/$USER/genome/${sample}/${ID}_aligned_reads.bam

java -jar \$EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics \
R=${genome} \
I=/scratch/$USER/genome/${sample}/${ID}_${ref}_sorted_reads.bam \
O=${outputfolder}/${sample}/${ID}_${ref}_alignment_metrics.txt \
TMP_DIR=/scratch/$USER/tempfiles/$HOME/.local/bin

java -jar \$EBROOTPICARD/picard.jar BuildBamIndex \
I=/scratch/$USER/genome/${sample}/${ID}_${ref}_sorted_reads.bam \
TMP_DIR=/scratch/$USER/tempfiles/

java -jar \$EBROOTPICARD/picard.jar CollectInsertSizeMetrics \
I=/scratch/$USER/genome/${sample}/${ID}_${ref}_sorted_reads.bam \
O=${outputfolder}/${sample}/${ID}_${ref}_insert_metrics.txt \
HISTOGRAM_FILE=${outputfolder}/${sample}/${ID}_${ref}_insert_size_histogram.pdf \
TMP_DIR=/scratch/$USER/tempfiles/

"> ${library}/${sample}/sbatches/${ID}_${sample}_alignment.sh

chmod 755 ${library}/${sample}/sbatches/${ID}_${sample}_alignment.sh

done

## 1st Sbatch bridge: important to create dependent jobs in a slurm cluster ##

bridge1=$(for i in $(ls ${library}/${sample}/sbatches/*_alignment.sh)
do
sbatch $i
done)

echo "submitting BAM Alignment for ${sample} job ID: $bridge1"

jobs1=$(echo $bridge1 | sed -e 's/ /:/g')



############################ 2nd SCRIPT: BAM wrap and quality filtering: merge bam files and marking duplicates ################################



echo "#!/bin/bash

##### merging BAM files per sample  #####

bamis=\$(ls /scratch/$USER/genome/${sample}/*sorted_reads.bam | sort | sed 's/.*/I= &/' | tr \"\\n\" \" \" )

module purge

module load R/3.5.0-foss-2018a-X11-20180131 picard/2.25.0-Java-11

java -Xmx96G -Xms96G -jar \$EBROOTPICARD/picard.jar MergeSamFiles \${bamis} O= /scratch/$USER/genome/${sample}/${sample}_concat.bam

### marking duplicates with PICARDTOOLS ###

java -Xmx96G -Xms96G -jar \$EBROOTPICARD/picard.jar MarkDuplicates \
I=/scratch/$USER/genome/${sample}/${sample}_concat.bam \
REMOVE_DUPLICATES=TRUE \
O=/scratch/$USER/genome/${sample}/${sample}_dedup_reads.bam \
METRICS_FILE=${outputfolder}/metrics/duplicates/${sample}_duplicates_metrics.txt \
TMP_DIR=/scratch/$USER/tempfiles/

java -Xmx96G -Xms96G -jar \$EBROOTPICARD/picard.jar BuildBamIndex \
I=/scratch/$USER/genome/${sample}/${sample}_dedup_reads.bam \
TMP_DIR=/scratch/$USER/tempfiles/

### Collecting Final Alignment metrics ###

java -Xmx96G -Xms96G -jar \$EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics \
R=${genome} \
I=/scratch/$USER/genome/${sample}/${sample}_dedup_reads.bam \
O=${outputfolder}/metrics/alignment/${sample}_concat_alignment_metrics.txt \
TMP_DIR=/scratch/$USER/tempfiles/

java -Xmx96G -Xms96G -jar \$EBROOTPICARD/picard.jar CollectInsertSizeMetrics \
I=/scratch/$USER/genome/${sample}/${sample}_dedup_reads.bam \
O=${outputfolder}/metrics/alignment/${sample}_concat_insert_metrics.txt \
HISTOGRAM_FILE=${outputfolder}/metrics/alignment/${sample}__concat_insert_size_histogram.pdf \
TMP_DIR=/scratch/$USER/tempfiles/

java -Xmx96G -Xms96G -jar \$EBROOTPICARD/picard.jar CollectWgsMetrics \
R=${genome} \
I=/scratch/$USER/genome/${sample}/${sample}_dedup_reads.bam \
O=${outputfolder}/metrics/alignment/${sample}_collect_wgs_metrics.txt \
TMP_DIR=/scratch/$USER/tempfiles/

java -Xmx96G -Xms96G -jar \$EBROOTPICARD/picard.jar QualityScoreDistribution \
I=/scratch/$USER/genome/${sample}/${sample}_dedup_reads.bam \
O=${outputfolder}/metrics/alignment/${sample}_qual_score_dist.txt \
CHART=${outputfolder}/metrics/alignment/${sample}_qual_score_dist.pdf \
TMP_DIR=/scratch/$USER/tempfiles/

module load SAMtools/1.10-GCC-9.3.0

samtools coverage /scratch/$USER/genome/${sample}/${sample}_dedup_reads.bam > ${outputfolder}/${sample}/metrics/${sample}_coverage_estimates

" > ${outputfolder}/${sample}/sbatches/${sample}_bamwrapduplicates.sh

chmod 755 ${outputfolder}/${sample}/sbatches/${sample}_bamwrapduplicates.sh

## be aware MarkDuplicates needs lots of RAM memory ##

bridge2=$(sbatch ${outputfolder}/${sample}/sbatches/${sample}_bamwrapduplicates.sh)

echo "submitting BAM quality check for ${sample} job ID:$bridge2"

job2=$(echo $bridge2 | sed -e 's/ /:/g')



#################################### 3rd Script: First Calling for recalibration by chromosssomes ##################################################################




mkdir /scratch/$USER/genome/gvcfs/${sample}/

for word in $(cat ${chromofile})

do

echo "#!/bin/bash

module load GATK/4.1.7.0-GCCcore-9.3.0-Java-11 R/3.6.1-foss-2018a

gatk --java-options \"-Xmx16G -Xms16G\" HaplotypeCaller \
-R ${genome} \
-I /scratch/$USER/genome/${sample}/${sample}_dedup_reads.bam \
-O /scratch/$USER/genome/${sample}/gvcfs/${sample}_${word}.g.vcf.gz \
 -L ${word} \
 -ERC GVCF \
--minimum-mapping-quality 30 \
--native-pair-hmm-threads 8
"> ${library}/${sample}/sbatches/${sample}_${word}_firstvariantcalling.sh

chmod 755 ${library}/${sample}/sbatches/${sample}_${word}_firstvariantcalling.sh

done

## bridge2 = creating waiting list for the posterior sbatches once all chormossomes calls are done ###

bridge3=$(for i in $(ls ${library}/${sample}/sbatches/*_firstvariantcalling.sh)
do
sbatch $i
done)

echo "Submitting per chromossome variant calling for ${sample} job IDs: $bridge3"

jobs3=$(echo $bridge3 | sed -e 's/ /:/g')




################################## 4th SCRIPT: Creating posterior recalibration job ################################################################################



echo "#!/bin/bash

module load GATK/4.1.7.0-GCCcore-9.3.0-Java-11 R/3.6.1-foss-2018a

gvcfs=\$(ls /scratch/$USER/genome/${sample}/gvcfs/*.g.vcf.gz | sort | sed 's/.*/-V= &/' | tr \"\\n\" \" \" )

gatk --java-options \"-Xmx96g -Xms96g\" CombineGVCFs \
-R ${genome} \
\${gvcfs} \
-O /scratch/$USER/genome/${sample}/gvcfs/${sample}_concat.g.vcf.gz

gatk --java-options \"-Xmx96g -Xms96g\" GenotypeGVCFs \
-R ${genome} \
-V /scratch/$USER/genome/${sample}/gvcfs/${sample}_concat.g.vcf.gz \
-O /scratch/$USER/genome/${sample}/vcfs/${sample}_variants.vcf.gz

gatk --java-options \"-Xmx96G -Xms96G\" SelectVariants \
-R ${genome} \
-V /scratch/$USER/genome/${sample}/vcfs/${sample}_variants.vcf.gz \
-select-type SNP \
-O /scratch/$USER/genome/${sample}/vcfs/${sample}_raw_snps.vcf.gz

gatk --java-options \"-Xmx96G -Xms96G\" SelectVariants \
-R ${genome} \
-V /scratch/$USER/genome/${sample}/vcfs/${sample}_variants.vcf.gz \
-select-type INDEL \
-O /scratch/$USER/genome/${sample}/vcfs/${sample}_raw_indels.vcf.gz

gatk --java-options \"-Xmx96G -Xms96G\" VariantFiltration \
-R ${genome} \
-V /scratch/$USER/genome/${sample}/vcfs/${sample}_raw_snps.vcf.gz \
--filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' \
--filter-name \"basic_snp_filter\" \
-O /scratch/$USER/genome/${sample}/vcfs/${sample}_filtered_snps.vcf.gz

gatk --java-options \"-Xmx96G -Xms96G\" VariantFiltration \
-R ${genome} \
-V /scratch/$USER/genome/${sample}/vcfs/${sample}_raw_indels.vcf.gz \
--filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' \
--filter-name \"basic_indel_filter\" \
-O /scratch/$USER/genome/${sample}/vcfs/${sample}_filtered_indels.vcf.gz

gatk --java-options \"-Xmx96G -Xms96G\" SelectVariants \
-R ${genome} \
-V /scratch/$USER/genome/${sample}/vcfs/${sample}_filtered_snps.vcf.gz \
--exclude-filtered \
-O /scratch/$USER/genome/${sample}/vcfs/${sample}_snps_recalib.vcf.gz

gatk --java-options \"-Xmx96G -Xms96G\" SelectVariants \
-R ${genome} \
-V /scratch/$USER/genome/${sample}/vcfs/${sample}_filtered_indels.vcf.gz \
--exclude-filtered \
-O /scratch/$USER/genome/${sample}/vcfs/${sample}_indels_recalib.vcf.gz

gatk --java-options \"-Xmx96G -Xms96G\" BaseRecalibrator \
-I /scratch/$USER/genome/${sample}/${sample}_dedup_reads.bam \
-R ${genome} \
--known-sites /scratch/$USER/genome/${sample}/vcfs/${sample}_snps_recalib.vcf.gz \
--known-sites /scratch/$USER/genome/${sample}/vcfs/${sample}_indels_recalib.vcf.gz \
-O /scratch/$USER/genome/${sample}/${sample}_recal_data.table

gatk --java-options \"-Xmx96G -Xms96G\" ApplyBQSR \
-R ${genome} \
-I /scratch/$USER/genome/${sample}/${sample}_dedup_reads.bam \
-bqsr /scratch/$USER/genome/${sample}/${sample}_recal_data.table \
-O /scratch/$USER/genome/${sample}/${sample}_recal_reads.bam

gatk --java-options \"-Xmx96G -Xms96G\" BaseRecalibrator \
-R ${genome} \
-I /scratch/$USER/genome/${sample}/${sample}_recal_reads.bam \
-O /scratch/$USER/genome/${sample}/${sample}_postrecal.table \
--known-sites /scratch/$USER/genome/${sample}/vcfs/${sample}_snps_recalib.vcf.gz \
--known-sites /scratch/$USER/genome/${sample}/vcfs/${sample}_indels_recalib.vcf.gz

gatk --java-options \"-Xmx96G -Xms96G\" AnalyzeCovariates \
-before /scratch/$USER/genome/${sample}/${sample}_recal_data.table \
-after /scratch/$USER/genome/${sample}/${sample}_postrecal.table \
-plots ${outputfolder}/${sample}/metrics/${sample}_recalplots.pdf

## copping the final recalibrated BAM file to the output data directory

cp /scratch/$USER/genome/${sample}/${sample}_recal_reads.bam ${outputfolder}/bams/${sample}_recal_reads.bam

cp /scratch/$USER/genome/${sample}/${sample}_recal_reads.bai ${outputfolder}/bams/${sample}_recal_reads.bai

"> ${library}/${sample}/sbatches/${sample}_posteriordep.sh

chmod 755 ${library}/${sample}/sbatches/${sample}_posteriordep.sh

## running steps after haplotype caller runs (applying bridge 2 after all chromossomes are done) ##

bridge4=$(sbatch ${library}/${sample}/sbatches/${sample}_posteriordep.sh)

echo "Recalibrating BAMS for ${sample} JobID: $bridge4"

job4=$(echo $bridge4 | sed -e 's/ /:/g')


############################ 5th Script: Applying Rescaled BAMs to call variants through the GATK in ERC mode #######################################



for word in $(cat ${chromofile})

do

echo "#!/bin/bash

module purge

module load GATK/4.1.7.0-GCCcore-9.3.0-Java-11 R/3.6.1-foss-2018a

gatk --java-options \"-Xmx16G -Xms16G\" HaplotypeCaller \
-R ${genome} \
-I /scratch/$USER/genome/${sample}/${sample}_recal_reads.bam \
-O /scratch/$USER/genome/${sample}/gvcfs/${sample}_${word}_rec_variants.g.vcf.gz \
-L ${word} \
-ERC GVCF \
--output-mode EMIT_ALL_ACTIVE_SITES \
--minimum-mapping-quality 30 \
--native-pair-hmm-threads 8

" > ${library}/${sample}/sbatches/${sample}_${word}_rec_variantcalling.sh

chmod 755 ${library}/${sample}/sbatches/${sample}_${word}_rec_variantcalling.sh

done

bridge5=$(for i in $(ls ${library}/${sample}/sbatches/*_rec_variantcalling.sh)
do
sbatch $i
done)

echo "Submitting final calling for ${sample} per chromossome $bridge5"

jobs5=$(echo $bridge5 | sed -e 's/ /:/g')


##################################### 6th SCRIPT: combining chromossomes to create a single GVCF file  ##############################################################################


echo "#!/bin/bash

module load GATK/4.1.7.0-GCCcore-9.3.0-Java-11 R/3.6.1-foss-2018a

gvcfs=\$(ls /scratch/$USER/genome/${sample}/gvcfs/*_rec_variants.g.vcf.gz | sort | sed 's/.*/-V= &/' | tr \"\\n\" \" \" )

gatk --java-options \"-Xmx96g -Xms96g\" CombineGVCFs \
-R ${genome} \
\${gvcfs} \
-O ${outputfolder}/gvcfs/${sample}_concat.g.vcf.gz

"> ${library}/${sample}/sbatches/${sample}_variantcalling.sh

chmod 755 ${library}/${sample}/sbatches/${sample}_variantcalling.sh

sbatch ${library}/${sample}/sbatches/${sample}_variantcalling.sh

echo "Submitting final calling steps for ${sample}"

rm ${library}/${sample}/c1
rm ${library}/${sample}/c2
rm ${library}/${sample}/libraries

done
