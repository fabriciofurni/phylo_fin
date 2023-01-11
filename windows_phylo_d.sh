### Running phylogenies and Patterson's D locus by locus ####

## Runs are performed by chromossome in a loop ###


for chromo in $(cat $chromofile | cut -f1)

do


cat << EOF > ${chromo}_wind.sh
#!/bin/bash

cat $genomic_windows | grep ${chromo} | while read line

do

chr=\$(echo \$line | cut -d " " -f1)
sta=\$(echo \$line | cut -d " " -f2)
end=\$(echo \$line | cut -d " " -f3)

mkdir \${chr}

mkdir \${chr}/\${chr}_\${end}


## filtering locus from the whole variant data 

module load BCFtools

bcftools view -Oz phylogeny_auto_filt_snps.vcf.gz -r \${chr}:\${sta}-\${end}  > \${chr}/\${chr}_\${end}/\${chr}_\${end}.vcf.gz


## Converting from vcf to phylip with vcf2phylip


/data/p290631/softwares/vcf2phylip.py -i \${chr}/\${chr}_\${end}/\${chr}_\${end}.vcf.gz --output-prefix \${chr}_\${end} --output-folder \${chr}/\${chr}_\${end}/

cat \${chr}/\${chr}_\${end}/\${chr}_\${end}.min4.phy | sed 's/\*/-/g' > \${chr}/\${chr}_\${end}/\${chr}_\${end}.final.phy

rm \${chr}/\${chr}_\${end}/\${chr}_\${end}.min4.phy


### Running phylogenetic estimates with RaxML


module load RAxML-NG/1.0.2-gompi-2020b

raxml-ng-mpi -all --msa \${chr}/\${chr}_\${end}/\${chr}_\${end}.final.phy --model GTR+G --bs-tree 100 --outgroup AN020050 -prefix \${chr}/\${chr}_\${end}/\${chr}_\${end} --threads 24 --workers 4 --extra thread-nopin --force perf_threads


### Calculating D-statistic values using Dsuite Dtrios for the locus 

/data/p290631/softwares/Dsuite/Build/Dsuite Dtrios -k 100 ${chr}/\${chr}_\${end}/\${chr}_\${end}.vcf ${chr}/\${chr}_\${end}/sample_pop.txt

done

EOF
