#!/usr/bin/env Rscript

# 26.03.2013 18:57:58 EDT
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
args = matrix(c('input'  , 'i', 1, "character", "Tab-delimited input with mapped read locations (chr, start, end, strand)",
                'output' , 'o', 2, "character", "Name of output file (optional)",
                'V'      , 'v', 2, "integer",   "V parameter value",
                'tB'     , 'b', 2, "integer",   "tB parameter value",
                'window' , 'w', 2, "integer",   "Window size (default: 50)",
                'help'   , 'h', 0, "logical",   "Brief help message"
               ), ncol=5,byrow=T);
opt  = getopt(args);

# Specify default argument values
if ( is.null(opt$output) ) { opt$output = sub("\\..*$", ".hmm", opt$input, perl=T) }
if ( is.null(opt$window) ) { opt$window = 50 }

# Help message
if ( !is.null(opt$help) || is.null(opt$input)) {
   self=sub("^.*/", "", commandArgs()[4], perl=T)
   cat("\n", getopt(args, command=self, usage=T), "\n")
   q(status=1);
}

#############
# LIBRARIES #
#############
require(GROseqCalls)


#############
# FUNCTIONS #
#############

# iterateHMM
#
# Run an iteration of the HMM model with a particular set of
# parameters and write the output to a bed file
iterateHMM <- function(E2alldata, LtProbB, UTS, prefix) {
   BothStrands = DetectTranscriptsEM(E2alldata, LtProbB=LtProbB, UTS=UTS, thresh=1, debug=TRUE)

   ## Write out hg18 transcripts.
   write.table(BothStrands[[4]], file=pipe(paste("gzip -c > ", prefix, ".B", LtProbB, ".V", UTS, ".bed.gz", sep="")), quote=F, row.names=F, col.names=F, sep="\t")
}


########
# MAIN #
########

# Read input bed file and order by probe position
print(paste("Doing all adata in one var; tB=",opt$tB,"V=",opt$V))
data = read.table(opt$input, header=F, colClasses=c("character","integer","integer","character"), sep="\t", check.names=F, stringsAsFactors=F)
print("Read data")
data = data[ order(data[,1], data[,2], data[,3]) ,];
print("Ordered data")

# Run the transition from a transcript --> no transcript over a grid
iterateHMM(data, LtProbB=opt$tB, UTS=opt$V,prefix=opt$output)

