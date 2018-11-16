
Danni: 

#!/bin/bash
  
echo "analyzing SannichInlet_data"
refPATH=/projects/micb405/resources/project_2/2018/SaanichInlet_100m/MetaBAT2_SaanichInlet_100m/MedQPlus_MAGs/
for fasta in $refPATH/*.fa;
do
prefix=$(basename $fa | sed 's/.fa//g')
domain=$( grep -w $prefix path/to/gtdbtk.tsv | awk '{ print $2 }' | awk -F';' '{ print $1 }' | sed 's/d__//g' )
prokka \
--outdir Prokka_output/ \
--prefix $prefix \
--kingdom $domain \
--cpus 4 \
$fasta

done

~                                                                               
~      
