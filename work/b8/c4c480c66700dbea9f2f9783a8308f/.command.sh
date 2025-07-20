#!/bin/bash -ue
bwa index Cauris.fa
bwa mem Cauris.fa data/processed/Cauris_trimmed.fastq > Cauris.sam
samtools view -bS Cauris.sam > Cauris.bam
samtools sort Cauris.bam -o Cauris.sorted.bam
samtools index Cauris.sorted.bam
freebayes -f Cauris.fa Cauris.sorted.bam > data/processed/Cauris.vcf
