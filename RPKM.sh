#!/bin/bash

echo "RPKM - SI042_100m"
/projects/micb405/resources/project_2/2018/rpkm \
-c /home/dzhu_mb18/project2/Prokka_output2/SI100m_concatenated_ffn.fasta \
-a /home/dzhu_mb18/project2/BWA_output/SI042_100m.qtrim.artifact.rRNA.clean.sam \
-o /home/dzhu_mb18/project2/RPKM_output/SI042_SaanichInlet_ORFs_RPKM.csv

echo "RPKM - SI048_100m"
/projects/micb405/resources/project_2/2018/rpkm \
-c /home/dzhu_mb18/project2/Prokka_output2/SI100m_concatenated_ffn.fasta \
-a /home/dzhu_mb18/project2/BWA_output/SI048_100m.qtrim.artifact.rRNA.clean.sam \
-o /home/dzhu_mb18/project2/RPKM_output/SI048_SaanichInlet_ORFs_RPKM.csv

echo "RPKM - SI072_100m"
/projects/micb405/resources/project_2/2018/rpkm \
-c /home/dzhu_mb18/project2/Prokka_output2/SI100m_concatenated_ffn.fasta \
-a /home/dzhu_mb18/project2/BWA_output/SI072_100m.qtrim.artifact.rRNA.clean.sam \
-o /home/dzhu_mb18/project2/RPKM_output/SI072_SaanichInlet_ORFs_RPKM.csv

echo "RPKM - SI074_100m"
/projects/micb405/resources/project_2/2018/rpkm \
-c /home/dzhu_mb18/project2/Prokka_output2/SI100m_concatenated_ffn.fasta \
-a /home/dzhu_mb18/project2/BWA_output/SI074_100m.qtrim.artifact.rRNA.clean.sam \
-o /home/dzhu_mb18/project2/RPKM_output/SI074_SaanichInlet_ORFs_RPKM.csv

echo "RPKM - SI075_100m"
/projects/micb405/resources/project_2/2018/rpkm \
-c /home/dzhu_mb18/project2/Prokka_output2/SI100m_concatenated_ffn.fasta \
-a /home/dzhu_mb18/project2/BWA_output/SI075_100m.qtrim.artifact.rRNA.clean.sam \
-o /home/dzhu_mb18/project2/RPKM_output/SI075_SaanichInlet_ORFs_RPKM.csv

