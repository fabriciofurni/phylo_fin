#!/bin/bash

help_menu () {
	echo ""
	echo ""
	echo "Required arguments:"
 	echo ""
	echo "-g <path to bwa indexed reference genome>"
	echo "-c <list with chormossome file>"
 	echo "-l <path to population folder with gvcfs files to be applied in downstream analyses>"
 	echo "-p <project code/ prefix to be assigned in the output files>"
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
elif [ -z "${chromofile}" ]
then
	echo ""
	echo "no value provided for -c"
	help_menu
	exit 1
elif [ -z "${outputfolder}" ]
then
	echo ""
	echo "no value provided for -l"
	help_menu
	exit 1
elif [ -z "${proj}" ]
then
  echo ""
  echo "no value provided for -p"
  help_menu
  exit 1
}

while getopts :g:c:l:p:d:h FLAG; do
	case "${FLAG}" in
		g) 	genome=${OPTARG}
			;;
		c) 	chromofile=${OPTARG}
			;;
		l) 	outputfolder=${OPTARG}
			;;
    p)  proj=${OPTARG}
      ;;
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

## creating outputfolder for pop analyses ##

mkdir ${outputfolder}/sbatches/
mkdir ${outputfolder}/gvcfs/chrm/
mkdir ${outputfolder}/vcfs/
mkdir ${outputfolder}/metrics/

################################# Combine sample gvcfs and genotype per chromossome ########################################################

## splitting by chromossome increases efficiency ##

for word in $(cat ${chromofile})

do

echo "#!/bin/bash

gvcfs=\$(ls ${outputfolder}/gvcfs/*.g.vcf.gz | sort | sed 's/.*/-V= &/' | tr \"\\n\" \" \" ) #working

module load GATK/4.1.7.0-GCCcore-9.3.0-Java-11

gatk --java-options \"-Xmx16g -Xms16g\" CombineGVCFs \
-R ${genome} \
\${gvcfs} \
-L ${word} \
-O ${outputfolder}/gvcfs/chrm/${proj}_${word}_concat.g.vcf.gz

gatk --java-options \"-Xmx16g -Xms16g\" GenotypeGVCFs \
-R ${genome} \
-V ${outputfolder}/gvcfs/chrm/${proj}_${word}_concat.g.vcf.gz \
-L ${word} \
--include-non-variant-sites \
-O ${outputfolder}/vcfs/${proj}_${word}.vcf.gz

"> ${outputfolder}/sbatches/${proj}_${word}_comb.sh

chmod 755 ${outputfolder}/sbatches/${proj}_${word}_comb.sh

done



#### Creating new script for filtering data #####


echo "#!/bin/bash

vcfs=\$(ls ${outputfolder}/vcfs/*.vcf.gz | sort | sed 's/.*/-I &/' | tr \"\\n\" \" \" )

module load R/3.5.0-foss-2018a-X11-20180131 picard/2.25.0-Java-11

/scratch/$USER/softwares/gatk-4.1.8/gatk GatherVcfs \
\${vcfs} \
-O ${outputfolder}/vcfs/${proj}_concat.vcf.gz \
-RI true

module load HTSlib

tabix ${outputfolder}/vcfs/${proj}_concat.vcf.gz -f

module purge

### filtering dataset 

module load BCFtools
module load HTSlib

bcftools stats -s - ${outputfolder}/vcfs/${proj}_concat.vcf.gz > ${outputfolder}/metrics/${proj}_raw_vcf_all.txt

## filtering SNPs with quality thresholds ##

bcftools view -M2 -e F_MISSING>0.2 -Oz ${outputfolder}/vcfs/${proj}_concat.vcf.gz  | bcftools view --exclude-types indels -e FORMAT/DP<3 > ${outputfolder}/vcfs/${proj}_snps_filtered.vcf.gz

tabix ${outputfolder}/vcfs/${proj}_snps_filtered.vcf.gz

## filtering autosomal SNPs ######

bcftools view --threads 16 -Oz -t ^CM020962.1,^CM020963.1,^CM018075.1${outputfolder}/vcfs/${proj}_snps_filtered.vcf.gz > ${outputfolder}/vcfs/${proj}_snps_auto_filtered.vcf.gz

tabix ${outputfolder}/vcfs/${proj}_snps_auto_filtered.vcf.gz

## filtering X chromosome SNPS #######

bcftools view ${outputfolder}/vcfs/${proj}_snps_filtered.vcf.gz -r CM020962.1 > ${outputfolder}/vcfs/${proj}_X_chrom_snps_filtered.vcf.gz

tabix ${outputfolder}/vcfs/${proj}_X_chrom_snps_filtered.vcf.gz

### filtering male Y-chrom SNPs #####

bcftools view -S ${malelist} ${outputfolder}/vcfs/${proj}_snps_filtered.vcf.gz -r CM020963.1 > ${outputfolder}/vcfs/${proj}_Y_chrom_snps_filtered.vcf.gz

tabix ${outputfolder}/vcfs/${proj}_Y_chrom_snps_filtered.vcf.gz

" > ${outputfolder}/sbatches/${proj}_filtering.sh

chmod 755 ${outputfolder}/sbatches/${proj}_filtering.sh

done



