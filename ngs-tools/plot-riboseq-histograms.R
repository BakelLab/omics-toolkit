#!/usr/bin/env Rscript

# 04.12.2019 10:28:20 EST
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
                'output' , 'o', 2, "character", "Prefix of output file (optional)",
                'help'   , 'h', 0, "logical",   "Brief help message"
               ), ncol=5,byrow=T);
opt  = getopt(args);

# Specify default argument values
if ( is.null(opt$output) ) { opt$output = "histogram" }

# Help message
if ( !is.null(opt$help) || is.null(opt$input)) {
   self=sub("^.*/", "", commandArgs()[4], perl=T)
   cat("\n", getopt(args, command=self, usage=T), "\n")
   q(status=1);
}

#############
# LIBRARIES #
#############

suppressMessages(library(ggplot2));
suppressMessages(library(dplyr));

########
# MAIN #
########

# opt = list(input="ribo_ltm_vir_riboseq-histogram.txt", output="ribo_ltm_vir")

# Read histogram data
d = read.delim(opt$input, header=T)

# Make a faceted plot on linear scale
pdf(file=paste(opt$output, "_riboseq_facet_counts.pdf", sep=""))
ggplot(d, aes(x=Position, y=Count, fill=Segment)) +
   geom_bar(stat="identity") + 
   facet_wrap( ~ Segment, ncol=2, scales="free") + 
   theme(legend.position="none") + 
   xlim(-20, 80) + 
   geom_vline(xintercept=0, linetype="dashed") + 
   geom_vline(xintercept=-16, linetype="dotted") + 
   geom_vline(xintercept=-10, linetype="dotted") +
   ggtitle(paste(opt$output, "- Read counts"))

# Make faceted plot, log scale with all axes the same
d$Count = log10(d$Count)
ggplot(d, aes(x=Position, y=Count, fill=Segment)) + 
   geom_bar(stat="identity") + 
   facet_wrap( ~ Segment, ncol=2) + 
   theme(legend.position="none") + 
   xlim(-20, 80) +
   geom_vline(xintercept=0, linetype="dashed") + 
   geom_vline(xintercept=-16, linetype="dotted") + 
   geom_vline(xintercept=-10, linetype="dotted") +
   ggtitle(paste(opt$output, "- Log10(Read counts)"))

graphics.off()



