#!/bin/bash

echo "analyzing SannichInlet_data"
refPATH=/projects/micb405/resources/project_2/2018/SaanichInlet_100m/MetaBAT2_SaanichInlet_100m/MedQPlus_MAGs/ 
for fasta in $refPATH/*.fa;
do
prefix=$(basename $fasta | sed 's/.fa//g')
domain=$( grep -w $prefix /projects/micb405/resources/project_2/2018/SaanichInlet_100m/MetaBAT2_SaanichInlet_100m/gtdbtk_output/gtdbtk.*.classification_pplacer.tsv | awk '{ print $2 }' | awk -F';' '{ print $1 }' | sed 's/d__//g' )
echo $prefix, $domain
prokka \
--outdir ~/project2/Prokka_output2/ \
--prefix $prefix \
--kingdom $domain \
--cpus 4 \
--force \
$fasta

done

