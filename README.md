LaggingStrand
=============

Exploring the relationships between DNA replication and mutation. The FIG*
subdirectories contain the code used to generate figures in the manuscript 
Reijns et al, Nature 2015 (doi:10.1038/nature14183). Associated 
dependencies and code for data preparation are located in the main 
directory.

Queries
-------
The code here was written to address specific questions, it has not been
polished or optimised for reuse. It is provided as is, as a record of
analyses undertaken. However, some aspects of the work do have potential
for reuse and we are keen to encourage and support that reuse. Please email
martin.taylor@igmm.ed.ac.uk if you have specific queries.

Data
----
Data used in this analysis and it's sources are summarised in the file
dataSources.txt. EmRiboSeq data generated in this work is available from
GEO under accession number GSE64521. The script ionProcessor.sh documents
primary processing of the emRiboSeq fastq files and depends on the script
riboReadEndBedToGidChrM.pl.

Yeast polymorphism, between-species relative substitution rate, genomic
masking and transcript annotation are available as a pre-processed table
(doi:).

GID
---
Yeast genome analysis was performed in R using large data-frames with
one row per nucleotide of the yeast genome. GID (gid) used throughout
scripts denotes an integer index of the sacCer3 reference genome,
effectively concatenating each of the chromosome reference sequences.
GID is a 1-based coordinate system with no padding between chromosomes.

Memory requirements
-------------------
The large data-frame approach to analysis used here requires substantial
RAM for analysis. Loading all relevant data sources including the
emRiboSeq data requires approximately 250GB RAM.


