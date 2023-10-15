#!/usr/bin/env Rscript

# 29.03.2013 15:10:16 EDT
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

#############
# FUNCTIONS #
#############

RUNTEST <- function(O100) {
   tID = paste(O100[[1]], O100[[2]], O100[[3]], O100[[6]], sep="-")
   gID = paste(O100[[7]], O100[[8]], O100[[9]], O100[[12]], sep="-")

   ## First, simply determine the number of times that a transcript is repeated.
   gREP = NROW(gID) -NROW(unique(gID))
   tREP = NROW(tID) -NROW(unique(tID))

   ## As transcripts are repeated many times, report the 
   ##  number of genes with at least ONE transcript.
   ngID = unlist(lapply(as.character(unique(gID)), function(x) NROW(which(gID == x)) ))
   ntID = unlist(lapply(as.character(unique(tID)), function(x) NROW(which(tID == x)) ))
   return(c(gREP, sum(ngID>1), NROW(unique(gID)), tREP, sum(ntID>1), NROW(unique(tID))))
}

########
# MAIN #
########

data = read.table(opt$input, sep="")
cat(c(opt$input, RUNTEST(data)), file=opt$output, sep="\t")
cat("\n", file=opt$output, append=T)
