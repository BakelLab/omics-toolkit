#!/usr/bin/env Rscript

# 29.11.2018 19:55:16 EST
# Harm van Bakel <hvbakel@gmail.com>

#############
# ARGUMENTS #
#############

#  Col 1: the long flag name. A multi-character string.
#  Col 2: short flag alias of Column 1. A single-character string.
#  Col 3: 0=no argument, 1=required argument, 2=optional argument.
#  Col 4: Possible values: logical, integer, double, complex, character.
#  Col 5: A brief description of the purpose of the option.
library(getopt)
args = matrix(c('input'     , 'i', 1, "character", "Input file with track file information (2 cols, track-name and track-file, no header)",
                'list'      , 'l', 1, "character", "Provide either geneIDs (1 col), geneIDs and geneNames (2 cols), or coordinates (4 cols bed format) to plot",
                'reference' , 'r', 1, "character", "Name of the reference track (optional)",
                'ucscgenome', 'g', 1, "character", "Set genome for optional UCSC tracks to add to plots. Default: hg38",
                'ucsctracks', 'u', 1, "character", "Optional comma-separated list of UCSC tracks to add. Options are SegDups,KnownGenes,RefGenes,EnsGenes,CpgIslands,Conservation,GCcontent",
                'zoom'      , 'z', 1, "numeric",   "Flanking region to include as a fraction of the size of the selected gene locus (value between 0 and 100). Default: 0.3",
                'datatype'  , 't', 1, "character", "Gviz data plot type. Default: histogram",
                'output'    , 'o', 2, "character", "Name of output file (optional)",
                'help'      , 'h', 0, "logical",   "Brief help message"
               ), ncol=5,byrow=T);
opt  = getopt(args);

# Specify default argument values
if ( is.null(opt$output)     ) { opt$output     = "track-plots.pdf" }
if ( is.null(opt$zoom)       ) { opt$zoom       = 0.3               }
if ( is.null(opt$datatype)   ) { opt$datatype   = "h"               }
if ( is.null(opt$ucscgenome) ) { opt$ucscgenome = "hg38"            }

# Help message
if ( !is.null(opt$help) || is.null(opt$input) || is.null(opt$list) ) {
   self=sub("^.*/", "", commandArgs()[4], perl=T)
   cat("\n", getopt(args, command=self, usage=T), "\n")
   q(status=1);
}

#############
# LIBRARIES #
#############

suppressMessages(library(Gviz, quietly=T));

track.colors = c("#FFD58A","#ADD8E6","#BFDB72","#D08EC0","#929EC8","#DBC496")

#############
# FUNCTIONS #
#############

ucscSegDups <- function( genome=character(), chr=character(), from=integer(), to=integer() ) {
   sd = UcscTrack(genome = genome, chromosome = chr, 
                  track = "genomicSuperDups", from = from, to = to,
                  trackType = "AnnotationTrack", 
                  start = "chromStart", end = "chromEnd", 
                  id = "name", shape = "box", fill = "#006400", 
                  stacking = "dense", name = "SegDups")
   return(sd)
}

ucscKnownGenes <- function( genome=character(), chr=character(), from=integer(), to=integer() ) {
   kg = UcscTrack(genome = genome, chromosome = chr, 
                  track = "knownGene", from = from, to = to,
                  trackType = "GeneRegionTrack", 
                  rstarts = "exonStarts", rends = "exonEnds", 
                  gene = "name", symbol = "name", 
                  transcript = "name", strand = "strand", 
                  fill = "#8282d2", name = "UCSC Genes")
   return(kg)
}

ucscRefGenes <- function( genome=character(), chr=character(), from=integer(), to=integer() ) {
   rg = UcscTrack(genome = genome, chromosome = chr,
                  track = "xenoRefGene", from = from, to = to,
                  trackType = "GeneRegionTrack", 
                  rstarts = "exonStarts", rends = "exonEnds", 
                  gene = "name",  symbol = "name2", 
                  transcript = "name", strand = "strand",
                  fill = "#8282d2", stacking = "dense", 
                  name = "Other RefSeq")
   return(rg)
}

ucscEnsGenes <- function( genome=character(), chr=character(), from=integer(), to=integer() ) {
   eg = UcscTrack(genome = genome, chromosome = chr, 
                  track = "ensGene", from = from, to = to,
                  trackType = "GeneRegionTrack", 
                  rstarts = "exonStarts", rends = "exonEnds",
                  gene = "name", symbol = "name2", 
                  transcript = "name", strand = "strand", 
                  fill = "#960000", name = "Ensembl Genes")
   return(eg)
}

ucscCpgIslands <- function( genome=character(), chr=character(), from=integer(), to=integer() ) {
   ci = UcscTrack(genome = genome, chromosome = chr, 
                  track = "cpgIslandExt", from = from, to = to,
                  trackType = "AnnotationTrack", 
                  start = "chromStart", end = "chromEnd", 
                  id = "name", shape = "box", fill = "#006400", 
                  name = "CpG Islands")
   return(ci)
}

ucscConservation <- function( genome=character(), chr=character(), from=integer(), to=integer() ) {
   co = UcscTrack(genome = genome, chromosome = chr, 
                  track = "Conservation", 
                  table = "phyloP30wayPlacental",
                  from = from, to = to, trackType = "DataTrack", 
                  start = "start", end = "end", data = "score",
                  type = "hist", window = "auto", 
                  col.histogram = "darkblue", 
                  fill.histogram = "darkblue",
                  ylim = c(-3.7, 4), name = "Conservation")
   return(co)
}

ucscGCcontent <- function( genome=character(), chr=character(), from=integer(), to=integer() ) {
   gc = UcscTrack(genome = genome, chromosome = chr, 
                  track = "GC Percent", table = "gc5Base",
                  from = from, to = to, trackType = "DataTrack", 
                  start = "start", end = "end", data = "score",
                  type = "hist", window = -1, windowSize = 1500, 
                  fill.histogram = "black", col.histogram = "black",
                  ylim = c(30, 70), name = "GC Percent")
   return(gc)
}


########
# MAIN #
########

# opt = opt=list(input="track_inputs.txt", list='reclocus-cluster-regions_subset.txt')

# Load track input file and gene list
d.track = read.delim(opt$input, header=F, stringsAsFactors=F, check.names=F)
d.list  = read.delim(opt$list, header=F, stringsAsFactors=F, check.names=F)

# Format UCSC
if ( !is.null(opt$ucsctracks) ){
   opt$ucsctracks = strsplit( opt$ucsctracks ,",")[[1]]
   stopifnot( all( opt$ucsctracks %in% c("SegDups","KnownGenes","RefGenes","EnsGenes","CpgIslands","Conservation","GCcontent") ) )
}

# Check input
if (ncol(d.track) != 2){
   stop("Track file should have two columns (name, track-file)\n");
}
if ( any(duplicated(d.track[,1])) ){
   stop("Duplicate track names found\n");
}
if ( any(duplicated(d.track[,2])) ){
   stop("Duplicate track file names found\n");
}

# Read track list object
col.count = 1;
options(ucscChromosomeNames=FALSE)
track.list = list();
for (trackid in 1:nrow(d.track)){
   track.name = d.track[trackid,1];
   track.file = d.track[trackid,2];
   
   if (file.exists(track.file)){
      if (col.count < length(track.colors)){
         fill = track.colors[col.count]
      } else {
         fill = "darkgrey"
      }
      
      # Load track type based on file extension
      if ( length( grep("(\\.gtf|\\.gff3)$", track.file, ignore.case=T, perl=T) ) > 0 ){
         track.list[[track.name]] = GeneRegionTrack(track.file, name=track.name, collapse=F, min.width=0.5, lwd=0.5, col=fill, fill=fill);
         col.count = col.count + 1
      } else if ( length( grep("(\\.bed|\\.gff|\\.gff1|\\.gff2)$", track.file, ignore.case=T, perl=T) ) > 0 ) {
         track.list[[track.name]] = AnnotationTrack(track.file, name=track.name, collapse=F, min.width=0.5, lwd=0.5, col=fill, fill=fill);
         col.count = col.count + 1
      } else if ( length( grep("(\\.bedgraph|\\.wig|\\.bigwig|\\.bw)$", track.file, ignore.case=T, perl=T) ) > 0 ) {
         track.list[[track.name]] = DataTrack(track.file, name=track.name, type=opt$datatype, baseline=0, col.histogram=NULL, col=NULL);
      } else if ( length( grep("(\\.fa|\\.fasta|\\.2bit)$", track.file, ignore.case=T, perl=T) ) > 0 ) {
         track.list[[track.name]] = SequenceTrack(track.file, name=track.name);
      } else if ( length( grep("(\\.bam)$", track.file, ignore.case=T, perl=T) ) > 0 ) {
         track.list[[track.name]] = AlignmentsTrack(track.file, name=track.name);
      } else {
         warning( paste("Skipping unsupported track file:", track.file, "\n", sep=" ") );
      }
   } else {
      warning( paste("Track file", track.file, "does not exist; skipping\n", sep=" ") );
   }
}

# Select the reference track
if ( is.null(opt$reference) ){
   reference.track = names(track.list)[1]
} else {
   if ( opt$reference %in% names(track.list) ) {
      reference.track = opt$reference
   } else {
      stop( paste("Reference track", opt$reference, "does not exist in track list\n", sep=" ") );
   }
}

# Cycle through provided gene IDs and make plot
track.list.bkp = track.list;  # Keep the current state of the track list since we will be modifying it for plotting
pdf(opt$output, width=14, height=20)
for (i in 1:nrow(d.list)){

   # Subset the reference track to just the region to be plotted
   ref     = track.list[[reference.track]];
   if ( ncol(d.list) == 4 ){          # Coordinates provided (chr, start, end, name)
      main = d.list[i,4]
      ref  = ref[seqnames(ref)==d.list[i,1] & end(ref)>=d.list[i,2] & start(ref)<=d.list[i,3] & feature(ref) %in% c("exon")]
   } else if ( ncol(d.list) == 5 ) {  # Coordinates provided with higlight genes (chr, start, end, name, highlight genes)
      main = d.list[i,4]
      ref  = ref[seqnames(ref)==d.list[i,1] & end(ref)>=d.list[i,2] & start(ref)<=d.list[i,3] & feature(ref) %in% c("exon")]
      # Highlight features in tracks
      hlgt = strsplit( d.list[i,5] ,"\\|")[[1]]
      for(track in 1:length(track.list)){
         if( class(track.list[[track]])[1] %in% c("GeneRegionTrack", "AnnotationTrack")  ) {
            sub.sel = group(track.list[[track]]) %in% hlgt;
            feature(track.list[[track]])[sub.sel] = "highlight"
         }
      }
   } else if ( ncol(d.list) == 1 ) {  # Gene IDs provided
      gene = d.list[i,1];
      main = gene;
      ref  = ref[gene(ref) == gene & feature(ref) %in% c("exon")]
      sub.sel = gene(ref) == gene
      feature(track.list[[reference.track]])[sub.sel] = "highlight"
   } else if ( ncol(d.list) == 2 ) {  # Gene IDs and display name provided
      gene = d.list[i,1];
      main = paste(d.list[i,1], " (", d.list[i,2], ")", sep="");
      ref  = ref[gene(ref) == gene & feature(ref) %in% c("exon")]
      sub.sel = gene(ref) == gene
      feature(track.list[[reference.track]])[sub.sel] = "highlight"
   } else {
      cat("Skipping\n");
   }

   # Plot progress
   cat(paste("Preparing plot for: ", main, "\n", sep=""));   
   
   # Proceed if we have any hits
   if ( length( seqlevels(ref)) == 1 ) {
      ref.chr   = seqlevels(ref)[1];
      chromosome(ref) = ref.chr;
      ref.start = min( c(start(ref), end(ref)) )
      ref.end   = max( c(start(ref), end(ref)) )
      ref.size  = abs(ref.end - ref.start)
      zoom.start = floor(ref.start-(opt$zoom*ref.size))
      zoom.end   = ceiling(ref.end+(opt$zoom*ref.size))
      if(zoom.start < 0) zoom.start = 0
      cat(opt$ucscgenome, ref.chr, zoom.start, zoom.end, "\n")
      
      # Make a plot if we have any results, otherwise insert empty page listing results not found
      tryCatch(
      {
         
         # Gather data for optional UCSC tracks
         ucsc.tracks = list()
         if ( !is.null(opt$ucsctracks) ){
            for( ut in 1:length(opt$ucsctracks) ){
               if (opt$ucsctracks[ut]=="SegDups"){
                  ucsc.tracks[[ut]] = ucscSegDups( genome=opt$ucscgenome, chr=ref.chr, from=zoom.start, to=zoom.end )
               } else if (opt$ucsctracks[ut]=="KnownGenes"){
                  ucsc.tracks[[ut]] = ucscKnownGenes( genome=opt$ucscgenome, chr=ref.chr, from=zoom.start, to=zoom.end )
               } else if (opt$ucsctracks[ut]=="RefGenes"){
                  ucsc.tracks[[ut]] = ucscRefGenes( genome=opt$ucscgenome, chr=ref.chr, from=zoom.start, to=zoom.end )
               } else if (opt$ucsctracks[ut]=="EnsGenes"){
                  ucsc.tracks[[ut]] = ucscEnsGenes( genome=opt$ucscgenome, chr=ref.chr, from=zoom.start, to=zoom.end )
               } else if (opt$ucsctracks[ut]=="CpgIslands"){
                  ucsc.tracks[[ut]] = ucscCpgIslands( genome=opt$ucscgenome, chr=ref.chr, from=zoom.start, to=zoom.end )
               } else if (opt$ucsctracks[ut]=="Conservation"){
                  ucsc.tracks[[ut]] = ucscConservation( genome=opt$ucscgenome, chr=ref.chr, from=zoom.start, to=zoom.end )
               } else if (opt$ucsctracks[ut]=="GCcontent"){
                  ucsc.tracks[[ut]] = ucscGCcontent( genome=opt$ucscgenome, chr=ref.chr, from=zoom.start, to=zoom.end )
               } else {
                  # Uknown track type
               }
            }
            track.list = append(ucsc.tracks, track.list)
         }
         track.list = append( list(axis=GenomeAxisTrack()), track.list )
         ht = HighlightTrack(trackList = track.list, start = ref.start, width=ref.size, chromosome = ref.chr, col="red", fill=NA)
         plotTracks(ht, chromosome=ref.chr, from=zoom.start, to=zoom.end, col=NULL, main=main, collapse=F, min.width=0.5, highlight="red", groupAnnotation="group", transcriptAnnotation = "symbol");         
      },
      error=function(e){ # Just print exception message
         message(paste("Error generating plot for: ", main, "\n", sep=""));
         plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n', xaxt='n', yaxt='n', xlab='', ylab='')
         text(1,4, paste("Plotting error for: ", main, "\n", sep=" ") , pos=4)
         print(e);
      }
      );
   
   } else {
      # We didn't get any hits so we make an empty plot stating that
      plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n', xaxt='n', yaxt='n', xlab='', ylab='')
      text(1,4, paste("No data for: ", main, "\n", sep=" ") , pos=4)
   }
   
   # Restore the state of the track list
   track.list = track.list.bkp;
}
graphics.off()
