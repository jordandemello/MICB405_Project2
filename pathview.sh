#!/bin/bash

for f in /home/dzhu_mb18/project2/Prokka_output2/*ffn
do
prokka_id=$( head -1 $f | awk -F_ '{ print $1 }' | sed 's/^>//g' )
mag_id=$( echo $f | sed 's/.ffn//g')
echo $prokka_id,$mag_id
done >Prokka_MAG_map.csv
