#!/bin/bash

echo "generating sam files"

bwa mem -t 4 ~/project2/SaanichInlet_100m_index /projects/micb405/resources/project_2/2018/Metatranscriptomes/SI042_100m.qtrim.artifact.rRNA.clean.fastq > ~/project2/BWA_output/SI042_100m.qtrim.artifact.rRNA.clean.sam 

bwa mem -t 4 ~/project2/SaanichInlet_100m_index /projects/micb405/resources/project_2/2018/Metatranscriptomes/SI048_100m.qtrim.artifact.rRNA.clean.fastq.gz > ~/project2/BWA_output/SI048_100m.qtrim.artifact.rRNA.clean.sam

bwa mem -t 4 ~/project2/SaanichInlet_100m_index /projects/micb405/resources/project_2/2018/Metatranscriptomes/SI072_100m.qtrim.artifact.rRNA.clean.fastq.gz > ~/project2/BWA_output/SI072_100m.qtrim.artifact.rRNA.clean.sam

bwa mem -t 4 ~/project2/SaanichInlet_100m_index /projects/micb405/resources/project_2/2018/Metatranscriptomes/SI074_100m.qtrim.artifact.rRNA.clean.fastq.gz > ~/project2/BWA_output/SI074_100m.qtrim.artifact.rRNA.clean.sam

bwa mem -t 4 ~/project2/SaanichInlet_100m_index /projects/micb405/resources/project_2/2018/Metatranscriptomes/SI075_100m.qtrim.artifact.rRNA.clean.fastq.gz > ~/project2/BWA_output/SI075_100m.qtrim.artifact.rRNA.clean.sam
