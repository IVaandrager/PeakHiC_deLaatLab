#!/bin/bash

cd /delaat/group/geert/WAPL_HAP1/chipseq/FASTQ
bwa aln -t 6 /data2/repository/genome/hg19/nohaplotype/hg19.fa SRR5266523.fastq.gz | bwa samse -n 1 /data2/repository/genome/hg19/nohaplotype/hg19.fa - SRR5266523.fastq.gz | samtools view -q 1 -ut /data2/repository/genome/hg19/nohaplotype/hg19.fa - | samtools sort -T /tmp/SRR5266523 -o /delaat/group/geert/WAPL_HAP1/chipseq/MAPPING/Wapl3_3_IgG.bam -

cd /delaat/group/geert/WAPL_HAP1/chipseq/MAPPING/
java -jar /delaat/group/iwan/programmes/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar MarkDuplicates --INPUT Wapl3_3_IgG.bam --OUTPUT Wapl3_3_IgG_marked.bam --METRICS_FILE Wapl3_3_IgG_picard.txt

cd /delaat/group/geert/WAPL_HAP1/chipseq/FASTQ
bwa aln -t 6 /data2/repository/genome/hg19/nohaplotype/hg19.fa SRR5266527.fastq.gz | bwa samse -n 1 /data2/repository/genome/hg19/nohaplotype/hg19.fa - SRR5266527.fastq.gz | samtools view -q 1 -ut /data2/repository/genome/hg19/nohaplotype/hg19.fa - | samtools sort -T /tmp/SRR5266527 -o /delaat/group/geert/WAPL_HAP1/chipseq/MAPPING/Wapl3_3_CTCF.bam -

cd /delaat/group/geert/WAPL_HAP1/chipseq/MAPPING/
java -jar /delaat/group/iwan/programmes/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar MarkDuplicates --INPUT Wapl3_3_CTCF.bam --OUTPUT Wapl3_3_CTCF_marked.bam --METRICS_FILE Wapl3_3_CTCF_picard.txt --ASSUME_SORT_ORDER=coordinate

cd /delaat/group/geert/WAPL_HAP1/chipseq/FASTQ
bwa aln -t 6 /data2/repository/genome/hg19/nohaplotype/hg19.fa SRR5266531.fastq.gz | bwa samse -n 1 /data2/repository/genome/hg19/nohaplotype/hg19.fa - SRR5266531.fastq.gz | samtools view -q 1 -ut /data2/repository/genome/hg19/nohaplotype/hg19.fa - | samtools sort -T /tmp/SRR5266531 -o /delaat/group/geert/WAPL_HAP1/chipseq/MAPPING/Wapl3_3_SMC1.bam -

cd /delaat/group/geert/WAPL_HAP1/chipseq/MAPPING/
java -jar /delaat/group/iwan/programmes/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar MarkDuplicates --INPUT Wapl3_3_SMC1.bam --OUTPUT Wapl3_3_SMC1_marked.bam --METRICS_FILE Wapl3_3_SMC1_picard.txt --ASSUME_SORT_ORDER=coordinate


####### MACS2 PEAK CALLING
cd /delaat/group/geert/WAPL_HAP1/chipseq/MAPPING/
macs2 callpeak -t Wapl3_3_CTCF_marked.bam -c Hap1_IgG_marked.bam -f BAM -g hs -n Wapl3_3_CTCF --outdir /delaat/group/geert/WAPL_HAP1/chipseq/MACS2/ -B --SPMR --nomodel -q 0.01

cd /delaat/group/geert/WAPL_HAP1/chipseq/MACS2/
macs2 bdgcmp -t Wapl3_3_CTCF_treat_pileup.bdg -c Wapl3_3_CTCF_control_lambda.bdg -o Wapl3_3_CTCF_fc.bdg -m FE

cd /delaat/group/geert/WAPL_HAP1/chipseq/MAPPING/
macs2 callpeak -t Wapl3_3_SMC1_marked.bam -c Wapl3_3_IgG_marked.bam -f BAM -g hs -n Wapl3_3_SMC1 --outdir /delaat/group/geert/WAPL_HAP1/chipseq/MACS2/ -B --SPMR --nomodel -q 0.01

cd /delaat/group/geert/WAPL_HAP1/chipseq/MACS2/
macs2 bdgcmp -t Wapl3_3_SMC1_treat_pileup.bdg -c Wapl3_3_SMC1_control_lambda.bdg -o Wapl3_3_SMC1_fc.bdg -m FE


export PATH=$PATH:/delaat/group/geert/WAPL_HAP1/chipseq/MACS2/
####### convert bedgraph to bigwig
./bdg2bw Wapl3_3_CTCF_fc.bdg hg19.len

./bdg2bw Wapl3_3_SMC1_fc.bdg hg19.len



#########################
#java -jar /delaat/group/iwan/programmes/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar MarkDuplicates --INPUT Hap1_SMC1.bam --OUTPUT Hap1_SMC1_marked.bam --METRICS_FILE Hap1_SMC1_picard.txt --ASSUME_SORT_ORDER=coordinate

####### MACS2 PEAK CALLING

#macs2 callpeak -t Hap1_SMC1_marked.bam -c Hap1_IgG_marked.bam -f BAM -g hs -n Hap1_SMC1 --outdir /delaat/group/geert/WAPL_HAP1/chipseq/MACS2/ -B --SPMR --nomodel -q 0.01

#export PATH=$PATH:/delaat/group/geert/WAPL_HAP1/chipseq/MACS2/

#macs2 bdgcmp -t Hap1_SMC1_treat_pileup.bdg -c Hap1_SMC1_control_lambda.bdg -o Hap1_SMC1_fc.bdg -m FE

####### convert bedgraph to bigwig
#export PATH=$PATH:/home/geert/localdev/prog/bedtools2/bin/

#./bdg2bw Hap1_SMC1_fc.bdg hg19.len

