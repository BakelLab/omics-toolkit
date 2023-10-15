#!/usr/bin/env Rscript

# 12.10.2016 08:29:58 EDT
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
if ( is.null(opt$output) ) { opt$output = "plot.pdf" }

# Help message
if ( !is.null(opt$help) || is.null(opt$input)) {
   self=sub("^.*/", "", commandArgs()[4], perl=T)
   cat("\n", getopt(args, command=self, usage=T), "\n")
   q(status=1);
}

########
# MAIN #
########

library(ggplot2)
library(reshape2)

d       = read.delim(opt$input, header=T);
if ("PCT_SELECTED_BASES" %in% colnames(d) ){
   bar.lab = gsub("\\/.*","",d[,1]);
   bar.lab = gsub("_S\\d+_L\\d+","",bar.lab);
   
   # Selected bases; on + near bait
   bar.dat = rbind(d$PCT_SELECTED_BASES, 1-d$PCT_SELECTED_BASES);
   colnames(bar.dat) = bar.lab;
   rownames(bar.dat) = c("On target","Off target");
   bar.dat = melt(bar.dat);
   colnames(bar.dat) = c("Capture","Rank","value");
   bar.all.plt = ggplot(bar.dat, aes(x = Rank, y = value, fill = Capture)) + 
               geom_bar(stat = "identity") + theme_bw() + scale_fill_grey(start=0, end=0.9) +
               theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
               xlab("Sample") +
               ylab("Fraction") +
               ggtitle("Selected bases (on + near bait)");
   
   # Selected bases; on bait only
   on.bait = d$ON_BAIT_BASES / d$PF_BASES_ALIGNED;
   bar.dat = rbind(on.bait, 1-on.bait);
   rownames(bar.dat) = c("On target","Off target");
   bar.dat = melt(bar.dat);
   colnames(bar.dat) = c("Capture","Rank","value");
   bar.onb.plt = ggplot(bar.dat, aes(x = Rank, y = value, fill = Capture)) + 
               geom_bar(stat = "identity") + theme_bw() + scale_fill_grey(start=0, end=0.9) +
               theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
               xlab("Sample") +
               ylab("Fraction") +
               ggtitle("Selected bases (on bait only)");
   
   pdf(file=opt$out);
   print(bar.all.plt);
   print(bar.onb.plt);
   graphics.off();

} else{
   stop("Error: missing PCT_SELECTED_BASES column\n");
}
