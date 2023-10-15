#!/bin/sh

# 15.06.2016 19:41:47 EDT
# Harm van Bakel <hvbakel@gmail.com>

# Convert file to transcript-level gff, with exon and CDS included
gtf2gff.pl -f exon,CDS gencode.v23.annotation.gtf > gencode.v23.annotation.gff

# Convert gff to bed file


# Get bedfile with 5' UTRs


# Get bedfile with 3' UTRs


# Convert to UTR lengths


# Merge everything together and include gene IDs

