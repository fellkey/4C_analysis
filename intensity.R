rm(list=ls())
options(stringsAsFactors = F)

library(GenomicRanges)

color=c("#D9352B","#4286B7","#5CAB58","#FF874A","#B9BB5A","#8679BE")

###############Function_elements
CHIP=c('SRR1057274_H3K27ac_peaks_intensity.narrowPeak',
       'SRR4031400_H3K27me3_peaks_intensity.narrowPeak',
       'SRR7782877_H3K4me1_peaks_intensity.narrowPeak',
       'SRR969473_RNAPolII_peaks_intensity.narrowPeak')

CHIPname=c('H3K27ac',
           'H3K27me3',
           'H3K4me1',
           'RNAPolII')

ATAC=c('SRR5063984_ATAC_NC1_peaks_intensity.narrowPeak')


ATACname=c('ATAC_NC1')

data=c(CHIP,ATAC,UCSC)
name=c(CHIPname,ATACname)
allreads=c(
  81198561,
  119751590,
  39476061,
  17793371,
  36553616
)

#####################################IS_analysis############################
IS=data.frame()
for (i in 1:5) {
  IS1=read.table(paste0("~/4c/newnew/IS_analysis/IS_",name[i],"_41.txt"))
  IS1$func=name[i]
  IS1$allreads=allreads[i]
  IS1$intensity=IS1$V4/(IS1$V3-IS1$V2)/IS1$allreads*1e9
  IS1=cbind(sample,IS1)
  IS=rbind(IS,IS1)
}

IS=IS[,c(3,4,5,6,7,1,2,8,9)]
colnames(IS)=c("chr","start","end","reads","func","sample","status","allreads","intensity")

#####################################4C_analysis############################
df=read.table("sample_hg19_17.txt",sep='\t')
sample=df$sample
chr=df$chr

setwd("~/4c/newnew/4C_analysis/")
peak=data.frame()
for(a in 1:17){
  for (i in 1:5) {
    peak1=read.table(paste0("~/4c/newnew/4C_analysis/",sample[a],"_",chr[a],"_4C_",name[i],"_intensity.txt"))
    peak1$sample=sample[a]
    peak1$status=status[a]
    peak1$func=name[i]
    peak1$allreads=allreads[i]
    peak1$intensity=peak1$V4/(peak1$V3-peak1$V2)/peak1$allreads*1e9
    peak1$intensity_mean=mean(peak1$intensity)
    peak=rbind(peak,peak1)
  }
}
