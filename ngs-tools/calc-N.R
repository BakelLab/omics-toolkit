#!/usr/bin/env Rscript

# 30.08.2010 22:22:38 EDT
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
                'col'    , 'c', 2, "integer",   "Column with fragment lengths. Default: 2",
                'raw'    , 'r', 0, "logical",   "Brief help message",
                'help'   , 'h', 0, "logical",   "Brief help message",
                'thresh' , 't', 1, "numeric",   "Threshold"
               ), ncol=5,byrow=T);
opt  = getopt(args);

# Specify default argument values
if ( is.null(opt$output) ) { opt$output = "" }
if ( is.null(opt$thresh) ) { opt$thresh = 150}
if ( is.null(opt$col)    ) { opt$col    = 2  }

# Help message
if ( !is.null(opt$help) || is.null(opt$input)) {
   self=sub("^.*/", "", commandArgs()[4], perl=T)
   cat("\n", getopt(args, command=self, usage=T), "\n")
   q(status=1);
}

########
# MAIN #
########

data = read.delim(file=opt$input);
subset.sort = as.numeric(sort(data[data[,opt$col]>=opt$thresh,opt$col]));
subset.csum = cumsum(subset.sort);
coverage = subset.csum[length(subset.csum)];
n90 = subset.sort[which.min(abs(subset.csum-(coverage * 0.1)))];
n50 = subset.sort[which.min(abs(subset.csum-(coverage * 0.5)))];
n10 = subset.sort[which.min(abs(subset.csum-(coverage * 0.9)))];
n05 = subset.sort[which.min(abs(subset.csum-(coverage * 0.95)))];
if ( is.null(opt$raw) ) {
   cat(sprintf("%-10s %s nt", "threshold", opt$thresh), "\n");
   cat(sprintf("%-10s %s ", "contigs", length(subset.csum)), "\n");
   cat(sprintf("%-10s %.2f Mb", "coverage", coverage/1000000), "\n");
   cat(sprintf("%-10s %.2f Kb", "n90", n90/1000), "\n");
   cat(sprintf("%-10s %.2f Kb", "n50", n50/1000), "\n");
   cat(sprintf("%-10s %.2f Kb", "n10", n10/1000), "\n");
   cat(sprintf("%-10s %.2f Kb", "n05", n05/1000), "\n");
   cat(sprintf("%-10s %.2f Kb", "min", min(subset.sort)/1000), "\n");
   cat(sprintf("%-10s %.2f Kb", "max", max(subset.sort)/1000), "\n");
   cat(sprintf("%-10s %.2f Kb", "mean", mean(subset.sort)/1000), "\n");
   cat(sprintf("%-10s %.2f Kb", "median", median(subset.sort)/1000), "\n");
} else {
   cat(sprintf("%s\t%s", "threshold", opt$thresh), "\n");
   cat(sprintf("%s\t%s", "contigs", length(subset.csum)), "\n");
   cat(sprintf("%s\t%s", "coverage", coverage), "\n");
   cat(sprintf("%s\t%s", "n90", n90), "\n");
   cat(sprintf("%s\t%s", "n50", n50), "\n");
   cat(sprintf("%s\t%s", "n10", n10), "\n");
   cat(sprintf("%s\t%s", "n05", n05), "\n");
   cat(sprintf("%s\t%s", "min", min(subset.sort)), "\n");
   cat(sprintf("%s\t%s", "max", max(subset.sort)), "\n");
   cat(sprintf("%s\t%.2f", "mean", mean(subset.sort)), "\n");
   cat(sprintf("%s\t%s", "median", median(subset.sort)), "\n");
}
