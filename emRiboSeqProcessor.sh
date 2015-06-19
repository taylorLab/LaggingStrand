#!/bin/bash
# Queries to martin.taylor@igmm.ed.ac.uk
#####
# Depends on (but should work with other version numbers):
#  BEDTools 2.16.2
#  samtools 0.1.18
#  bowtie 2.0.0
#####
# Use as:
# ./emRiboSeqProcessor.sh <runName>
#
# Where <runName> is the prefix name of the fastq file obtained from
# the sequencing run.
#
# runName.tmp.bed is the final output file that should be used for
# browser visulisation and downstream analysis.
#####

export NM=$1
mkdir ${NM}
cd ${NM}
ln -s ../${NM}.fastq
# Assumes bowtie index of sacCer3 reference genome is in the starting directory.

# Use -p N to parallelise the alignment if more processor cores are available.
# Set by default to N=2 below.
bowtie2 -p 2 -x ../sacCer3 -U ${NM}.fastq -S ${NM}.sam
samtools view -b -S ${NM}.sam > ${NM}.master.bam
samtools view -b -q 30 -S ${NM}.sam > ${NM}.bam
bamToBed -i ${NM}.bam > ${NM}.bed
# This Perl code does the 5' end coordinate tansform. The sort component is required for ordering
# assumed in subsequent steps.
perl -alne 'if($F[5]=~s/\+/-/){$F[2]=$F[1];$F[1]--;}else{$F[5]=~s/-/\+/;$F[1]=$F[2];$F[2]++}print join "\t", @F' < ${NM}.bed | sort -k1,1 -k2,2n -k 6 > ${NM}.ribo.bed
# The group by counts the number of 5' ends per nucleotide and does so separately for each strand.
bedtools groupby -i ${NM}.ribo.bed -g 1,3,6 -c 1 -full -o count > ${NM}.counts.bed
# This function counts the number of genome mapped reads
export totalfpe=`wc -l ${NM}.ribo.bed | awk '{print $1}'`
cat ${NM}.counts.bed | perl -alne '$F[6]=($F[6]/$ENV{totalfpe})*1e6; print join "\t", @F;' > ${NM}.tmp.bed


