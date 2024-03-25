####################IS_intensity#########################
less sample_hg19_45.txt |awk '{print $2,$3-5000,$3+5000}' > sample_hg19_45_IS_10kb.txt
sed -i "s/ /\t/g" sample_hg19_45_IS_10kb.txt

#reads in peak
bedtools multicov -bams /home/data/t060307/4c/public_data/rmdup_hg19/SRR1057274_H3K27ac_trimmed.sorted.bam \
-bed ~/4c/newnew/all_sample_hg19_IS_10kb.txt \
> ~/4c/newnew/IS_analysis/IS_H3K27ac.txt 

bedtools multicov -bams /home/data/t060307/4c/public_data/rmdup_hg19/SRR4031400_H3K27me3_trimmed.sorted.bam \
-bed ~/4c/newnew/all_sample_hg19_IS_10kb.txt \
> ~/4c/newnew/IS_analysis/IS_H3K27me3.txt 

bedtools multicov -bams /home/data/t060307/4c/public_data/rmdup_hg19/SRR7782877_H3K4me1.sorted.bam \
-bed ~/4c/newnew/all_sample_hg19_IS_10kb.txt \
> ~/4c/newnew/IS_analysis/IS_H3K4me1.txt 

bedtools multicov -bams /home/data/t060307/4c/public_data/rmdup_hg19/SRR969473_RNAPolII_trimmed.sorted.bam \
-bed ~/4c/newnew/all_sample_hg19_IS_10kb.txt \
> ~/4c/newnew/IS_analysis/IS_RNAPolII.txt 

bedtools multicov -bams /home/data/t060307/4c/public_data/bulkATAC/align/SRR5063984_ATAC_NC1.last.bam \
-bed ~/4c/newnew/all_sample_hg19_IS_10kb.txt \
> ~/4c/newnew/IS_analysis/IS_ATAC_NC1.txt 

##################4C_peak_intensity####################
sample=$1
chr=$2
#find overlaps
bedtools intersect -a ${sample}_${chr}_peakC_4mb_unique.txt \
        -b SRR1057274_H3K27ac_peaks.narrowPeak > ${sample}_${chr}_peakC_H3K27ac.txt

bedtools intersect -a ${sample}_${chr}_peakC_4mb_unique.txt \
        -b SRR4031400_H3K27me3_peaks.narrowPeak > ${sample}_${chr}_peakC_H3K27me3.txt

bedtools intersect -a ${sample}_${chr}_peakC_4mb_unique.txt \
        -b SRR7782877_H3K4me1_peaks.narrowPeak > ${sample}_${chr}_peakC_H3K4me1.txt

bedtools intersect -a ${sample}_${chr}_peakC_4mb_unique.txt \
        -b SRR969473_RNAPolII_peaks.narrowPeak > ${sample}_${chr}_peakC_RNAPolII.txt

bedtools intersect -a ${sample}_${chr}_peakC_4mb_unique.txt \
        -b SRR5063984_ATAC_NC1_peaks.narrowPeak > ${sample}_${chr}_peakC_ATAC_NC1.txt

#reads in peak
bedtools multicov -bams SRR1057274_H3K27ac_trimmed.sorted.bam \
-bed ${sample}_${chr}_peakC_H3K27ac.txt \
> ${sample}_${chr}_peakC_H3K27ac_intensity.txt

bedtools multicov -bams SRR4031400_H3K27me3_trimmed.sorted.bam \
-bed ${sample}_${chr}_peakC_H3K27me3.txt \
> ${sample}_${chr}_peakC_H3K27me3_intensity.txt

bedtools multicov -bams SRR7782877_H3K4me1.sorted.bam \
-bed ${sample}_${chr}_peakC_H3K4me1.txt \
> ${sample}_${chr}_peakC_H3K4me1_intensity.txt

bedtools multicov -bams SRR969473_RNAPolII_trimmed.sorted.bam \
-bed ${sample}_${chr}_peakC_RNAPolII.txt \
> ${sample}_${chr}_peakC_RNAPolII_intensity.txt

bedtools multicov -bams SRR5063984_ATAC_NC1.last.bam \
-bed ${sample}_${chr}_peakC_ATAC_NC1.txt \
> ${sample}_${chr}_peakC_ATAC_NC1_intensity.txt
