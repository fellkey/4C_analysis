# segment genome into restriction enzyme fragment windows 
python ~/4c/HiC-Pro/bin/utils/digest_genome.py \
        -r ^AAATTT ^AAATTC ^GAATTT ^GAATTC \
        -o RE1.frag \
        /home/data/refdir/server/reference/genome/hg19/hg19.fa

python ~/4c/HiC-Pro/bin/utils/digest_genome.py \
        -r ^CATG \
        -o RE2.frag \
        /home/data/refdir/server/reference/genome/hg19/hg19.fa

# extract the 133bp reads after RE1
cat RE1.frag | awk 'BEGIN{OFS="\t";}$2!=0{print $1,$2+1,$2+6,"RE1";}' > RE1.motif
cat RE2.frag | awk 'BEGIN{OFS="\t";}$2!=0{print $1,$2+1,$2+4,"RE2";}' > RE2.motif
cat RE1.motif RE2.motif | sort -k1,1 -k2,2n > RE1_RE2.motif
awk 'BEGIN{OFS="\t";}{if(NR!=1) print LM,$0;LM=$0;}' RE1_RE2.motif | awk 'BEGIN{OFS="\t";}$1==$5{if($4==$8) print $1,$2,$7,$4,$8,$1":"$2-1"-"$2+132,$5":"$7-133"-"$7,"Blind";else print $1,$2,$7,$4,$8,$1":"$2-1"-"$2+132,$5":"$7-133"-"$7,"nonBlind";}' > Frags.txt
awk 'BEGIN{OFS="\t";}{print $1,$3-133,$3;print $1,$2-1,$2+132;}' RE1.motif > RE1_len133.bed
bedtools getfasta -fi /home/data/refdir/server/reference/genome/hg19/hg19.fa -bed RE1_len133.bed > RE1_len133.fa

# select only unique mapping reads
refGenomeFile=/home/data/refdir/server/reference/index/bowtie/hg19
fastaFile=RE1_len133.fa
cpu=10
bamFile=RE1_len133.bam
uniquesFile=RE1_len133_uniquesID.txt
bowtie2 -p ${cpu} -x ${refGenomeFile} -f ${fastaFile} | samtools view -q 1 -hbSu - > ${bamFile}
samtools view ${bamFile} | grep -v "XS:" | cut -f1 > ${uniquesFile}

# get the final frags files
awk 'BEGIN{OFS="\t";}NR==FNR{a[$1]=$1;next}{if($6 in a || $7 in a) print $1,$2,$3,$4,$5,"uniquesFrags";else print $1,$2,$3,$4,$5,"NuniquesFrags"}' ${uniquesFile} Frags.txt > hg19_ApoI_NlaIII.Frags.bed
awk 'BEGIN{OFS="\t";}{if($4=="RE1"&&$5=="RE2") print $0,"+",$2+6;else if($4=="RE2"&&$5=="RE1") print $0,"-",$3-6;else if($4=="RE1"&&$5=="RE1") print $0,"+",$2+6"\n"$0,"-",$3-6;else print $0,"+",$2;}' hg19_ApoI_NlaIII.Frags.bed > hg19_ApoI_NlaIII.all.Frags.bed

for ((i=1; i<=22;i++)) ;do (awk -v CHR=chr$i 'BEGIN{OFS="\t";}$1==CHR{if($4=="RE1"&&$5=="RE2") print $0,"+",$2+6;else if($4=="RE2"&&$5=="RE1") print $0,"-",$3-6;else if($4=="RE1"&&$5=="RE1") print $0,"+",$2+6"\n"$0,"-",$3-6;else print $0,"+",$2;}' hg19_ApoI_NlaIII.Frags.bed > hg19_ApoI_NlaIII.chr$i.Frags.bed);done
awk -v CHR=chrY 'BEGIN{OFS="\t";}$1==CHR{if($4=="RE1"&&$5=="RE2") print $0,"+",$2+6;else if($4=="RE2"&&$5=="RE1") print $0,"-",$3-6;else if($4=="RE1"&&$5=="RE1") print $0,"+",$2+6"\n"$0,"-",$3-6;else print $0,"+",$2;}' hg19_ApoI_NlaIII.Frags.bed > hg19_ApoI_NlaIII.chrY.Frags.bed
awk -v CHR=chrX 'BEGIN{OFS="\t";}$1==CHR{if($4=="RE1"&&$5=="RE2") print $0,"+",$2+6;else if($4=="RE2"&&$5=="RE1") print $0,"-",$3-6;else if($4=="RE1"&&$5=="RE1") print $0,"+",$2+6"\n"$0,"-",$3-6;else print $0,"+",$2;}' hg19_ApoI_NlaIII.Frags.bed > hg19_ApoI_NlaIII.chrX.Frags.bed
