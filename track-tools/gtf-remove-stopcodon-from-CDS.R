#!/usr/bin/env Rscript

# 18.10.2020 18:05:07 EDT
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
args = matrix(c('input'  , 'i', 1, "character", "Input file name",
                'output' , 'o', 2, "character", "Name of output file (optional)",
                'help'   , 'h', 0, "logical",   "Brief help message"
               ), ncol=5,byrow=T);
opt  = getopt(args);

# Specify default argument values
if ( is.null(opt$output) ) { opt$output = "" }

# Help message
if ( !is.null(opt$help) || is.null(opt$input)) {
   self=sub("^.*/", "", commandArgs()[4], perl=T)
   cat("\n", getopt(args, command=self, usage=T), "\n")
   q(status=1);
}


########
# MAIN #
########

gff0      <- rtracklayer::import(opt$input)
parentcol <- gff0$transcript_id

# for exon-exon junctions in the middle of a stop codon
removerows <- integer(0)

# fill in for top strand
CDS_top <- which(as.logical(gff0$type == "CDS" & BiocGenerics::strand(gff0) == "+"))
CDS_top_split <- split(CDS_top, parentcol[CDS_top])
modifyrows <- integer(length(CDS_top_split))
tosubtract <- integer(length(CDS_top_split))
for(i in seq_along(CDS_top_split)){
   cds <- CDS_top_split[[i]]
   last <- cds[which.max(BiocGenerics::start(gff0[cds]))]
   if(BiocGenerics::width(gff0[last]) > 3){
      # whole stop codon is in this exon (vast majority of cases)
      modifyrows[i] <- last
      tosubtract[i] <- 3L
   } else {
      # stop codon goes into the previous exon
      tosubtract[i] <- 3L - BiocGenerics::width(gff0[last]) # remaining nucleotides to remove
      removerows <- c(removerows, last)
      cds <- cds[cds != last]
      replacement = cds[which.max(start(gff0[cds]))]
      if ( length(replacement)>0 ){
         if ( !is.na(replacement) ){
            modifyrows[i] <- replacement
         }
      }
   }
}
BiocGenerics::end(gff0[modifyrows]) <- BiocGenerics::end(gff0[modifyrows]) - tosubtract

# fill in for bottom strand
CDS_bot <- which(as.logical(gff0$type == "CDS" & BiocGenerics::strand(gff0) == "-"))
CDS_bot_split <- split(CDS_bot, parentcol[CDS_bot])
modifyrows <- integer(length(CDS_bot_split))
tosubtract <- integer(length(CDS_bot_split))
for(i in seq_along(CDS_bot_split)){
   cds <- CDS_bot_split[[i]]
   last <- cds[which.min(BiocGenerics::start(gff0[cds]))]
   if(BiocGenerics::width(gff0[last]) > 3){
      # whole stop codon is in this exon (vast majority of cases)
      modifyrows[i] <- last
      tosubtract[i] <- 3L
   } else {
      # stop codon goes into the previous exon
      tosubtract[i] <- 3L - BiocGenerics::width(gff0[last]) # remaining nucleotides to remove
      removerows <- c(removerows, last)
      cds <- cds[cds != last]
      replacement = cds[which.min(BiocGenerics::start(gff0[cds]))]
      if ( length(replacement)>0 ){
         if ( !is.na(replacement) ){
            modifyrows[i] <- replacement
         }
      }
   }
}
BiocGenerics::start(gff0[modifyrows]) <- BiocGenerics::start(gff0[modifyrows]) + tosubtract

if(length(removerows) > 0) gff0 <- gff0[-removerows] # CDS that only contain stop codon

rtracklayer::export(gff0, opt$output, format = "gtf")
