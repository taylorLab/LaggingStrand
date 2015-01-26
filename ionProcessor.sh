#!/bin/bash
# Depends on (but should work with other version numbers):
#  BEDTools 2.16.2
#  samtools 0.1.18
#  bowtie 2.0.0
#####
# Use as:
# ./ionProcessor.sh <runName>
#
# Where <runName> is the prefix name of the fastq
#####
export NM=$1
cd ${NM}
# Assumes bowtie index of sacCer3 reference genome is in the parent directory.
bowtie2 -p 16 -x ../sacCer3 -U *.fastq -S ${NM}_sacCer3_bowtie2.e2e.sam
samtools view -b -q 30 -S ${NM}_sacCer3_bowtie2.e2e.sam > ${NM}_sacCer3_bowtie2.q30.bam
samtools sort ${NM}_sacCer3_bowtie2.q30.bam ${NM}_sacCer3_bowtie2.q30.sorted
samtools index ${NM}_sacCer3_bowtie2.q30.sorted.bam
bamToBed -cigar -i ${NM}_sacCer3_bowtie2.q30.bam > ${NM}_sacCer3_bowtie2.q30.bed
perl -w ../riboReadEndBedToGidChrM.pl ${NM}_sacCer3_bowtie2.q30.bed > ${NM}_sacCer3_bowtie2.q30.gid.withChrM


