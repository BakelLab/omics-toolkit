# omics-toolkit
Collection of scripts and utilities to analyze various types of 'omics data, organized by categories (Assembly tools, NGS tools, genome browser track tools, general text file processing utilities, and vcf tools).



## Installation

#### Github

To install directly from github, first clone this repository.

```bash
git clone git@github.com:BakelLab/omics-toolkit.git
```

Next, add the subdirectories within the `omics-toolkit` repository to your `PATH` environment variable. You can easily collect the ':' delimited string of all subdirectories by executing the following command from within the repository folder:

```bash
realpath */ | xargs | perl -pe 's/ /:/g'
```



## Assembly tools

| Script | Description |
| ---- | ---- |
| amos-layCheck.pl | Checks an amos layout file (such as those produced by the cabog pacBioToCA pipeline) for a variety of errors. |
| amos-layDropContigs.pl | Removes LAY entries for supplied contig identifiers so that they are not used to build amos banks |
| parse-hapcut2-blocks.pl | Parser to extract hapcut2 haplotype blocks. |
| soap2agp.pl | Convert soap files to AGP/FSA files suitable for submission to genbank. |

## NGS tools


| Script | Description |
| ---- | ---- |
| bam-insert-sizes.pl | Print insert sizes for paired reads in a bam file. Note that the input bam  file must be sorted by read pairs, not by coordinates. |
| bioanalyzer-pdf-parse-peak-tables.pl | Utility script to extract tab-delimited tables with peak sizes from bioanalyzer PDF files. |
| blat-coverage.pl | Gives information on the proportion of each query sequence that is covered by blat hits. |
| bowtie-insert-sizes.pl | Print insert sizes for paired reads in a bowtie output file. Note that the input bam  file must be sorted by read pairs, not by coordinates. |
| check-nucleotide-barcode-differences.pl | Count the minimum amount of base changes between all sequence pairs in an inputted barcode list. Used to check whether barcodes are sufficiently diverse. |
| extract-contigs_filter-blast.pl | Extract contigs based on a set of blast alignment filters. |
| fasta2tabbed.pl | Convert a fasta file to a tab-delimited format with the header in the first column and the sequence in the second column. |
| fasta-filter-by-id.pl | Filter a multi-fasta file by a list of identifiers. The identifiers are matched against the first word of the fasta header (i.e. all text before the first space, tab or line ending). |
| fasta-filter-by-size.pl | Filter a multi-fasta file by the length of the sequences it contains, to retain sequences above and/or below a certain length threshold. |
| fasta-get-info.pl | When provided with a multifasta file, this script will return the length, gc content, and N-content of each sequence, and, optionally, the length of the largest ORF on the forward and reverse strand |
| fasta-mask-pilon-failed.pl | Mask nucleotides that could not be corrected by 'pilon' by changing them to 'N' |
| fasta-orient-to-landmark.pl | Reorient a fasta file to a provided landmark sequence (e.g. an ORI sequence). |
| fasta-reflow.pl | Reflow sequences in a fasta file so that each sequence line is set to specified max length. |
| fasta-rename-ids.pl | Renames the identifiers in the fasta headers according to a reference file that contains a mapping between the current and a new identifier. |
| fasta-shredder.pl | Shred sequences in a (multi) fasta file into smaller segments with a specified amount of overlap between segments. |
| fasta-shuffle-order.pl | Randomly reorder sequences in a multifasta file. |
| fasta-sliding-window-gc.pl | Calculate the GC content of sequences within a multi-fasta file across a sliding window applied to each sequence. |
| fasta-splitter-bysize.pl | Split a file containing multiple fasta sequences into a set of smaller multifasta files, each containing a subset of sequences that add up to the approximate specified length. |
| fasta-splitter.pl | Split a file containing multiple fasta sequences into a set of smaller files, each containing no more than the specified number of sequences |
| fastq2tabbed.pl | Convert a fastq file to a tab-delimited format with columns containing read name, sequence and quality strings. |
| fastq_all2std.pl | Convert between a variety of different fastq quality encodings and sequencing file formats. |
| fastq-check-pairing.pl | Validate whether paired-end fastq files are properly paired and have no mismatching forward/reverse reads. |
| fastq-demultiplex-single.pl | Very basic script to demultiplex a single-end fastq file that has barcodes at the start of each read. |
| fastq-filter-by-id-paired.pl | Filter specific reads from a set of paired-end fastq files based on the common part of the read identifier. |
| fastq-filter-by-id.pl | Filter specific reads in a fastq file based on read ID matching |
| fastq-filter-by-seq.pl | Filter specific reads in a fastq file based on nucleotide (sub)sequence matching |
| fastq-filter-by-size.pl | Filter specific reads in a fastq file based on read lengths. |
| fastq-get-maxlength.pl | Find the maximum length of reads in a fastq file. |
| fastq-length-trimmer.pl | Trim reads in a fastq input file to a specified length by removing bases at the 5' and/or 3' ends. |
| fastq-multi-to-standard.pl | Converts a multi-line fastq file to the standard 4-lines per entry, standard form |
| fastq-quality-trimmer-paired.pl | Process paired-end fastq files to remove low-quality read ends. Assumes an illumina-type chemistry where the phred quality of bases progressively declines in a 5' to 3' direction. |
| fastq-quality-trimmer.pl | Process single-end fastq files to remove low-quality read ends. Assumes an illumina-type chemistry where the phred quality of bases progressively declines in a 5' to 3' direction. |
| fastq-standardize-headers.pl | Converts any fastq header to a standard format where each read ends in /1 or /2 to denote the orientation. |
| fastq-stats.pl | Gets read and base count stats from a fastq file. |
| megan2classification.pl | Convert megan output files to a tab-delimited file with taxonomic classifications. |
| pacbio_get-max-subreads.pl | Filter longest subreads from a PacBio filtered_subreads file. |
| parse-blastn-vector-screening.pl | Parse blastn output of a set of sequences against UniVec. |
| phymmblScoreReadsShm.pl | Wraps the phymmbl scoreReads script to run on /dev/shm to speed things up and to avoid temp files from mixing |
| qseq2fastq.pl | Converts illumina qseq files to fastq format. |
| riboseq-get-psite-offset-counts.pl | Estimate p-site offset counts in RiboSeq data |
| riboseq-get-psite-offset-tracks.pl | Create tracks of p-sites from RiboSeq data using offsets from the above script. |
| riboseq-get-transcriptome-psite-offset-tracks.pl | Create transcriptome tracks of p-sites from RiboSeq data. |
| riboseq-normalize-psite-offset-tracks.pl | Normalize psite offset tracks. |
| sff-extract-454-mates.pl | Extract 454 mate pairs from ssf formatted files. |
| calc-N.R | Calculate N statistics from a tab-delimited file with sequence lengths (e.g. the output of fasta-get-info.pl) |
| groseqHMM-predict-transcripts.R | Predict transcripts from GroSeq data |
| groseqHMM-QC-transcripts.R | Perform QC analysis of GroSeq transcripts. |
| pacbio_estimate-coverage-cost.R | Calculates the cost of pacbio sequencing for a genome of a particular size. |
| pairwise-hypergeometric.R | Perform pairwise hypergeometric tests for a set of overlaps. |
| plot-coverageBed.R | Plot bedfile coverage. |
| plot-expression-distributions.R | Plot gene expression distribution for limma outputs. |
| plot-gviz-genome-region.R | Make Gviz plots of genomic regions based on provided track data. |
| plot-heatmap.R | Make an expression heatmap plot from a tab-delimited file with heatmap data. |
| plot-hsmetrics.R | Plot HS metrics |
| plot-module-FET-matrix.R | Create a heatmap of FET enrichments based on overlaps of two sets of lists. |
| plot-riboseq-histograms.R | Plot psite histograms for riboseq data. |
| plot-travelling-ratio.R | Make travelling ratio plots for GroSeq data. |
| readstack-duplicate-stats.R | Assess the level of PCR duplicates by identifying stacks of reads with the exact start/end coordinates. |
| bam2fastq.sh | Wrapper script to convert a bam file back to fastq format. Reads are shuffled randomly to avoid biases in STAR and minimap2 mapping. |
| fasta-splitter-random.sh | Split a multi-fasta file, where each subset file contains a defined number of sequences that are randomly picked from the multi-fasta file. |
| flncbam2fastq.sh | Convert an flnc bam file to fastq format. |
| interproscan-distributed.sh | Wrapper script that takes a multi-fasta file as input and starts multiple interproscan jobs on the minerva cluster that each contain a defined subset of sequences. |
| pacbio_isoform-capture-enrichment.sh | Calculate the fold-enrichment in pacbio isoform capture experiments. |
| pfamscan-distributed.sh | Wrapper script that takes a multi-fasta file as input and starts multiple interproscan jobs on the minerva cluster that each contain a defined subset of sequences. |
| riboseq-profile-tracks.sh | Wrapper script to produce riboseq profile tracks. |
| split-bam-by-cell-barcodes.sh | Split a bam file by cell barcodes. |

## Genome browser track tools

| Script | Description |
| ---- | ---- |
| bar2wig.pl | Convert Affymetrix bar-formatted files to wig format. |
| bed2fasta.pl | Extract bed exon, transcript, gene sequences from a reference fasta file. |
| bed2gbrowse.pl | Takes a bed file and converts it into a file that can be loaded as a custom annotation track in gbrowse. |
| bed2gff_no-terminal-exons.pl | Takes a bed file and converts it into a gff file that only contains internal exons. |
| bed2gff.pl | Converts a bed file to a gff file |
| bed2intronexongff.pl | Takes a bed file and converts it into a gff file that contains both intron and exon features. |
| bed2intronflank.pl | Takes a bed file and extracts introns with defined flanking sequences. |
| bedgraph-extract-singleton-regions.pl | Extract singleton regions (i.e. areas with coverage=1) from a bedgraph file. |
| calc-gc-content.pl | Calculate gc content in a sliding window across a provided fasta sequence. |
| check-primers-against-library.pl | Check a set of primers against a library to see if they uniquely amplify only a single library sequence. |
| collapseBedGraph.pl | Merges consecutive regions in a bedgraph file together if they are within a specified distance of eachother. |
| collapseBed.pl | Takes one or more bed files with multiple overlapping transcripts and prepares a  non-redundant version such that each exon annotation only occurs in one transcript. |
| create-igb-genome.pl | Create a genome index that can be used for the IGB genome browser |
| fasta2bnib.pl | Converts fasta-formatted chromosome sequence files into binary sequence files for loading into IGB. |
| fixedWig2Sgr.pl | Converts a fixed-wig formatted file to affymetrix sgr format. |
| genpept2fasta.pl | Convert genpept files to fasta format. |
| get-antisense-measurements-matrix.pl | Get an expression measurement matrix of antisense strand features. |
| get-feature-data.pl | Intersect genomic features with genome-wide measurements to plot data values relative to the feature coordinates. |
| get-feature-overlap-pvals.pl | Calculate pvalues for the overlap between two feature sets |
| get_gene_orientations.pl | Calculate relative orientation between neighboring genes in a gff input file (i.e. diverging, converging, tandem, etc.) |
| get-relative-uas-gene-distances.pl | Get relative distances between genes and upstream activating sequences. |
| get-sense-measurements-matrix.pl | Get an expression measurement matrix of sense strand features. |
| gff2bed.pl | Convert a gff file to bed format. |
| gff2gbrowse.pl | Convert a gff file to gbrowse format. |
| gff2igb-php-links.pl | Convert a gff file to a php file of gene features, where clicking on a gene feature will automatically focus on that feature in a running IGB genome browser instance. |
| gff2length.pl | Calculate the processed/spliced length of each feature in a gff file. Also returns outermost genomic coordinates of each feature. |
| gff3-cds-to-full.pl | Convert a basic gff3 file containing only CDS annotations (i.e. as produced by RAST) to full specification with gene and exon descriptions and parent links. |
| gff3tobeddetail.pl | Convert a gff3 file to a beddetail formatted file. |
| gff-add-flanking-region.pl | Add flanking regions to transcripts within a gff file. |
| gff-count-attributes.pl | Count how many times a particular attribute value or combination of attribute values occurs in a gff file. |
| gff-get-terminal-exons.pl | Extract only the terminal exons from a gff file. |
| glimmer2bed.pl | Convert glimmer gene predictions to bed format. |
| gtf2gff.pl | Convert a gtf file to a gff file. |
| gtf-add-gene-entries.pl | Add 'gene' entry lines to a gtf file that lacks them |
| gtf-addreplace-attributes.pl | Add or replace attributes (key/value pairs) within a gtf file |
| gtf-count-attributes.pl | Count occurrences of one or more gtf attribute (key/value pairs) combinations. |
| gtf-filter-attributes.pl | Filter a gtf file on certain attributes (key/value pairs), e.g transcript_id. |
| gtf-get-gene-flanking-regions.pl | Extract gene-flanking regions from a gtf file. |
| gtf-get-gene-regions.pl | Extract gene regions (i.e. start/end coordinates of genes) from a gtf file. |
| intersectBedSgr.pl | Get Affymetrix sgr or bar file data points that overlap bed features. |
| iprscan2go.pl | Script to convert an interproscan file to a gaf formatted GO reference file to use for GO term enrichment analysis. |
| iprscan2oneline.pl | Compress interproscan results into a single line of interpro and GO annotations to facilitate unique file joins. |
| iprscan2topgoinput.pl | Convert interproscan output into a format that can be read by the bioconductor TopGO package to provide a custom GO reference file. |
| map-array-probes_arraydata2sgr.pl | Map sgr expression data onto a reference genome. |
| map-array-probes_blat2bpmap.pl | Convert blat files to bpmap format. |
| map-array-probes_blat2gff.pl | Convert microarray probe mapping blat output files to gff format. |
| map-array-probes_bpmap2tab.pl | Convert bpmap files to tab-delimited files. |
| map-array-probes_run-blat.pl | Map the location of microarray probes onto a reference genome using blat. |
| merge-feature-measurements-matrices.pl | Merge a set of feature-measurement matrices together to create a new matrix. |
| merge-gff-features-by-name.pl | Merge gff features by feature name. |
| peakmatcher.pl | Find overlap between peaks in the two files supplied |
| reads2sgr.pl | Convert sequence read coordinates in basic bed format to sgr format. |
| reads2wig.pl | Convert sequence read coordinates in basic bed format to wig format. |
| sgr2bar.pl | Convert sgr formatted files to bar format. |
| sgrcounts2fixedwincounts.pl | Convert an sgr file with count data at mapped positions in the genome to a file with counts in non-overlapping windows of a specified width. |
| shuffle-bed-features.pl | Shuffle feature coordinates within a specified set of regions (e.g. whole chromosomes or intergenic regions). |
| sjtab2bed.pl | Convert a STAR SJ.out.tab file to a bed-formatted junction file |
| sync-cluster-diagrams.pl | Sync cluster diagrams by rows and columns. |
| synteny-map_blocks.pl | Prepares a gff file that can be loaded into IGB to show the synteny blocks for another genome. |
| synteny-map_probes.pl | Map tiling array probe measurements from a syntenic target genome onto a reference chromosome. |
| tas_interval.pl | Tilecore interval analysis |
| tas_normalize.pl | Normalize affy cel files |
| variablestepwig2bedgraph.pl | Convert variable step wig files to bedgraph format. |
| wig2bar.pl | Convert wig files to Affymetrix bar format. |
| zeropad-binmatrix.pl | Zero pad a bin-matrix formatted file |
| gtf-remove-stopcodon-from-CDS.R | Convert a gtf file to remove the stop codon from CDS features. Needed for e.g. IsoSwitchAnalyzeR |
| plink-mds-mclust.R | Plink cluster analysis |
| plot-matrix-averagograms.R | Make averagograms across genomic features based on matrices of expression data. |
| plot-matrix-image.R | Make a matrix image plot across genomic features based on matrices of expression data. |
| plot-moving-averagograms.R | Plot a moving averagogram across genomic features based on input bam or expression files. |
| batch-gtf-remove-stopcodon-from-CDS.sh | Batch convert a gtf file to remove the stop codon from CDS features. Needed for e.g. IsoSwitchAnalyzeR |

## Utility

| Script | Description |
| ---- | ---- |
| add-header | Adds a header to any text file. The default separator for the header fields is a tab. Use -s to specify other separators. |
| author-contribution-generator | Convert tab-delimited file of author contributions to formatted text. First row must have a header that contains the contribution descriptions. |
| count-specific-characters | Count the number of occurrences of specific character(s) in a text file. |
| cut-by-ids | Extract selected columns from a delimited text file based on column header names. |
| diff-by-ids | Filters a file to return only the lines that match a set of identifiers specified in an ID file. |
| diff-cols-by-ids | Filters a file to return only the columns that match a set of identifiers specified in an ID file. |
| dim | Returns the line and column dimensions of a delimited text file |
| duplicated-ids | Returns a list of identifiers that occur more than once in a specified column of a delimited text file. |
| FET-compare-lists.R | Assess FET significance of overlaps between two lists picked from a set of background genes. |
| filter-matrix-cols | Filters the columns of a delimited count matrix text file based on a set of thresholds. Only those columns matching the filters will be retained. |
| filter-matrix-rows | Filters the rows of a delimited count matrix text file based on a set of thresholds. Only those columns matching the filters will be retained. |
| fjoin | FJoin finds overlaps (or, in general, proximity-based) pairs of features, given two feature sets. |
| groupBy_freqtab | Addon for the bedtools 'groupBy' utility that converts the freqdesc or freqasc output to a tabular format that is easier to interpret. |
| ids-common2all | Returns a list of identifiers that are shared between all input files |
| ids-common2both | Returns a list of identifiers that are shared between a column or row in two delimited input files |
| ids-uniq2left | Returns a list of identifiers found only in the first input file. |
| images2pdf | Convert a set of images to a pdf file that has one page per image. |
| intersect-by-ids | Filters a file to return only the lines that match a set of identifiers specified in an ID file. |
| intersect-cols-by-ids | Filters a file to return only the columns that match a set of identifiers specified in an ID file. |
| join-by-ids | Combines two files (line by line) based on a shared identifier. |
| join-by-position | Merge selected columns from the join file with the source file based on an overlap in position. |
| lesscol | Format a tab-delimited text file as a padded space-delimited text file. |
| matrix-to-boxplot | Converts the columns of a data matrix to a two-colum format (value, label) that can be used to make R boxplots. |
| multi-join-by-ids | Joins one column in a set of files together in a large matrix based on identifiers shared between all files |
| numbered-header | Get listing of header fields with column numbers. Utility script to easily identify the column numbers to process with awk etc. |
| reduce-pdf-filesize | Reduce the file size of a pdf document by downsampling images. |
| scale-column | Adjust the values of one column in a tab-delimited count matrix file by multiplying with a scaling factor. |
| scale-matrix | Adjust the values of a tab-delimited count data matrix by multiplying with a scaling factor. |
| sorted-intersect-by-ids | Filters a file to return only the lines that match a set of identifiers specified in an ID file. The lines in the filter file are re-sorted based on the order of the identifiers in the ID file. |
| split-on-column | Splits one or more files according to the identifiers in a specified column. One file is created for each unique identifier. |
| transpose | Transpose a delimited matrix text file, swapping rows and columns. |
| uniq-ids | Returns a list of unique identifiers in a row or column of a delimited text file. |

## Vcf tools

| Script | Description |
| ---- | ---- |
| vcf-count-variants.pl | Count variants in a vcf file |
| vcf-get-variants.pl | Extract variants from a vcf file. |