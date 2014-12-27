#!/bin/bash
# Assumes that this is run from the directory /array09a/martin/ionTorrent/ionProton
module add apps/gcc/BEDTools/2.16.2 apps/gcc/samtools/0.1.18 apps/gcc/bowtie/2.0.0
export NM=$1
cd ${NM}
bowtie2 -p 16 -x ../../sacCer3 -U *.fastq -S ${NM}_sacCer3_bowtie2.e2e.sam
samtools view -b -q 30 -S ${NM}_sacCer3_bowtie2.e2e.sam > ${NM}_sacCer3_bowtie2.q30.bam
samtools sort ${NM}_sacCer3_bowtie2.q30.bam ${NM}_sacCer3_bowtie2.q30.sorted
samtools index ${NM}_sacCer3_bowtie2.q30.sorted.bam
bamToBed -cigar -i ${NM}_sacCer3_bowtie2.q30.bam > ${NM}_sacCer3_bowtie2.q30.bed
perl -w ../../riboReadEndBedToGidChrM.pl ${NM}_sacCer3_bowtie2.q30.bed > ${NM}_sacCer3_bowtie2.q30.gid.withChrM
cat ${NM}_sacCer3_bowtie2.q30.gid.withChrM | perl -walne 'next if($F[0]<12071327);$F[0]-=12071326;print join"\t",@F;' > ${NM}_sacCer3_bowtie2.q30.gid.justChrM
cp  ${NM}_sacCer3_bowtie2.q30.gid.withChrM ../../RUNS/
