for fq in /projects/micb405/resources/project_2/2018/Metatranscriptomes/*100m*; 
do 
cruise=$( basename $fq | awk -F_ '{ print $1 }' ); 
echo $cruise; 
bwa mem -t 8 -p SI100m_concatenated.ffn $fq >tmp.sam; 
/projects/micb405/resources/project_2/2018/rpkm -c SI100m_concatenated.ffn -a tmp.sam -o $cruise\_100m_MetaT_rpkm.csv;
done
