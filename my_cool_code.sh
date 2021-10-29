#!/usr/bin/bash

-------------------------------------------------------------------------------------------------------------------------
###Step1-----------------------------------------------------------------------------------------------------------------

#initialise a git repository
git init

#copy the folder fastq to my working directory
cp -r /localdisk/data/BPSM/AY21/fastq .

#use fastqc to check the quality of raw data
fastqc fastq/*.fq.gz


-------------------------------------------------------------------------------------------------------------------------
###Step2-----------------------------------------------------------------------------------------------------------------






-------------------------------------------------------------------------------------------------------------------------
###Step3-----------------------------------------------------------------------------------------------------------------

# copy the folder Tcongo_genome to my working directory
cp -r /localdisk/data/BPSM/AY21/Tcongo_genome .
# build a database using bowtie2
bowtie2-build Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz T.cong

# align the read pairs to the Trypanosoma congolense genome
for i in {501..555}; 
do filename1=$(ls fastq/*${i}_1.fq.gz) ;
filename2=$( ls fastq/*${i}_2.fq.gz) ;
output_sam=${filename1/_1.fq.gz/.sam} ;
bowtie2 -x T.cong -1 $filename1 -2 $filename2 -S $output_sam; 
echo file $output_sam has been produced; 
done

#convert the “sam” to "bam" format
for i in {501..555}; 
do
input_sam=$(ls fastq/*${i}.sam) ;
output_bam=$(echo ${input_sam/sam/bam}) ;
samtools view -b -o $output_bam $input_sam ;
echo file $input_sam has been converted
done

#get indexed “.bam”
for i in {501..555}
do input=$(ls fastq/*${i}.bam) ;
output=$(echo ${input/100k/s100k}) ;
samtools sort $input > $output ;
samtools index $output ;
echo file $output has been indexed.
done


-------------------------------------------------------------------------------------------------------------------------
###Step4-----------------------------------------------------------------------------------------------------------------
#Question for Biologist: over time relative to WT controls

#copy the bedfile to my working directory
cp /localdisk/data/BPSM/AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed .

#find the filename and time value you need as a reference
awk '{FS="\t"; if($2 == "WT" && $5== "Uninduced"){print "fastq/s100k."$1".bam", $4;}}' fastq/100k.fqfiles | sort -k2,2n

#get the filelist for “.bam”
filelist=$(awk '{FS="\t"; if($2 == "WT" && $5== "Uninduced"){print "fastq/s100k."$1".bam";}}' fastq/100k.fqfiles | sort -k2,2n)

#generate counts data
bedtools multicov -bams $filelist -bed TriTrypDB-46_TcongolenseIL3000_2019.bed > WTovertime.txt

-------------------------------------------------------------------------------------------------------------------------
###Step5-----------------------------------------------------------------------------------------------------------------
##statistical mean (average) of the counts per gene

#cut the gene name and discription
cut -f 4,5 WTovertime.txt > WTovertime_gn.txt
#cut the value
cut -f 6-14 WTovertime.txt > WTovertime_v.txt

#give a header, “gene_name, gene_discription, mean_T0, mean_T24, mean_T48”
sed $'1igene_name\tgene_discription\tmean_T0\tmean_T24\tmean_T48' WTovertime_mean_data.txt > WTovertime_mean_table.txt



-------------------------------------------------------------------------------------------------------------------------
###Step6-----------------------------------------------------------------------------------------------------------------
##Fold change

#calculate the value of foldchange "T0 to T24", "T0 to T48", put them into the file "WTovertime_foldchange_data.txt"
awk '{
if(( $1+$2+$3 )!=0)
{print ( $4+$5+$6 ) / ( $1+$2+$3 )}
else
{print "NA"}
}' WTovertime_v.txt > WTovertime_T0T24.txt

awk '{
if(( $1+$2+$3 )!=0)
{print ( $7+$8+$9 ) / ( $1+$2+$3 )}
else
{print "NA"}
}' WTovertime_v.txt > WTovertime_T0T48.txt

paste WTovertime_gn.txt WTovertime_T0T24.txt WTovertime_T0T48.txt > WTovertime_foldchange_data.txt

#sort the WTovertime_foldchange_data.txt according to fold change value in reverse order
sort -t$'\t' -k3,3nr -k4,4nr WTovertime_foldchange_data.txt > WTovertime_foldchange_data_rv.txt


#give a header, “gene_name, gene_discription, foldchange_T24/T0, foldchange_T48/T0”
sed $'1igene_name\tgene_discription\tfoldchange_T24/T0\tfoldchange_T48/T0' WTovertime_foldchange_data_rv.txt > WTovertime_foldchange_table_rv.txt

