#!/bin/bash
fasta=$1
chr=$2
pos=$3

#remove PCR duplicates
seqkit fq2fa /Jurkat_4C_${fasta}_R1.fq -o ${fasta}.fa
seqkit rmdup -s -w 151 -j 10 ${fasta}.fa -o ${fasta}_nodup.fa

#select only reads with the correct primer sequence
sed -n '{/^ACTGGTGAGTACGCCAAAAATTC/{g;p}};h' \
${fasta}_nodup.fa \ 
>${fasta}_P1_R2.1.id

sed -n '{/^ACTGGTGAGTACGCCAAAAATTT/{g;p}};h' \
${fasta}_nodup.fa \
>${fasta}_P1_R2.2.id

cat ${fasta}_P1_R2.1.id ${fasta}_P1_R2.2.id >${fasta}_P1_R2.id
sed -i "s/>//g" ${fasta}_P1_R2.id
less ${fasta}_P1_R2.id|awk '{print$1}' > ${fasta}_P1_R1.id
rm ${fasta}_P1_R2*

seqkit grep --pattern-file ${fasta}_P1_R1.id \
/Jurkat_4C_${fasta}_R2.fq.gz \
-o ${fasta}_P1_R2.fq.gz 

seqkit grep --pattern-file ${fasta}_P1_R1.id \
/Jurkat_4C_${fasta}_R1.fq.gz \
-o ${fasta}_P1_R1.fq.gz

# map to hg19
bowtie2 -p 4 -x reference/index/bowtie/hg19 \
-1 ${fasta}_P1_R1.fq.gz -2 ${fasta}_P1_R2.fq.gz >${fasta}.sam
less ${fasta}.sam |grep "AS:" | grep -v "XS:" > ${fasta}.unique.sam
less ${fasta}.sam |grep @ >header
cat header ${fasta}.unique.sam >${fasta}_uunique.sam
samtools sort -O bam -@ 10 -o - ${fasta}_uunique.sam > ${fasta}_unique.sort.bam
samtools index ${fasta}_unique.sort.bam 
rm *sam

# quantify fragment coverage across the genome and in each sample
bedtools multicov -bams \
${fasta}_unique.sort.bam \
-bed frags/hg19_ApoI_NlaIII.${chr}.Frags.bed \
> ${fasta}.${chr}.unique.bedgraph

# running mean sliding window approach (21 fragends window)
Rscript running_mean.R ${fasta} ${chr}

# get the normalized bedgraph files 
less ${fasta}.${chr}.wSize21.unique.txt |awk '{print $1 "\t" $2 "\t" $3 "\t" $7}'>${fasta}.${chr}.wSize21.unique.bedgraph

#call peaks
Rscript call_peaks.R ${fasta} ${chr} ${pos}


