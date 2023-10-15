#!/usr/bin/env Rscript

# 25.06.2013 18:14:43 EDT
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
                'length' , 'l', 2, "integer",   "Read length",
                'help'   , 'h', 0, "logical",   "Brief help message"
               ), ncol=5,byrow=T);
opt  = getopt(args);

# Specify default argument values
if ( is.null(opt$output) ) { opt$output = opt$input }

# Help message
if ( !is.null(opt$help) || is.null(opt$input) || is.null(opt$length) ) {
   self=sub("^.*/", "", commandArgs()[4], perl=T)
   cat("\n", getopt(args, command=self, usage=T), "\n")
   q(status=1);
}


########
# MAIN #
########

data = read.table(opt$input, colClasses=c("character","integer","integer","numeric"), sep="\t", comment.char = "#", stringsAsFactors=F);
data = data[data[,3]-data[,2]<=opt$length,4];
pct.single = sprintf("%.2f", 100*(sum(data==1)/length(data)));
pct.multi  = sprintf("%.2f", 100*(sum(data>1)/length(data)));
medcov     = median(data);
opt$output = sub("\\.png|\\.PNG", "", opt$output, perl=T);
png(file=paste(opt$output, "png", sep="."), width=800, height=800)
hist(data[data<=50], breaks=seq(1,50), xlim=c(0,50), xlab="coverage", main="Readstack coverage histogram", freq=F)
legend("topright", legend=c(paste("% singleton\t", pct.single), paste("% multi read\t", pct.multi), paste("Median cov\t", medcov)));
dev=dev.off()
cat(paste("% singleton\t", pct.single, "\n", sep=""), paste("% multi read\t", pct.multi, "\n", sep=""), paste("Median cov\t", medcov, "\n", sep=""), sep="");
