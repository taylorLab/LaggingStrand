
File: gm2.sacCer3summaryAnnotationTable.Rdat
============================================
# Intended for loading into R for analysis.
# Each row corresponds to a single nucleotide of the sacCer3 reference genome.
# Coordinates are 1-based.
# Columns mean..
chr	  Chromosome number (chrM = 17).
chromPos  Position in chromosome.
gid	  Position in genome (coordinate in concatenated chromosomes).
forNuc	  Forward strand nucleotide [ATCG].
gcFor150  Mean G+C content for a window of this nucleotide and the next 149.
gcRev150  Mean G+C content for a window of this nucleotide and the preceeding 149.
polCount  Yeast polymorphism observed at this site [01].
polCall   Read coverage to call yeast polymorphisms at this site [01].
zoo	  Yeast between-species relative substitution rate (baseml under HKY) 1=genome wide average.
repTime	  Replication time measure (higher value, earlier replicating).
transcriptFor	Transcription annotated on the forward strand [01].
transcriptRev	Transcription annotated on the reverse strand [01].
exclude	  Regions excluded from analysis for non-uniqueness of mapping.

