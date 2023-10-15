#!/usr/bin/env Rscript

# 19.07.2012 14:37:45 EDT
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
args = matrix(c('input'     , 'i', 1, "character", "Input file name",
                'output'    , 'o', 2, "character", "Name of output file (optional)",
                'help'      , 'h', 0, "logical",   "Brief help message",
                'start'     , 's', 1, "numeric",   "Start of fragment length range (bp). Default: 200",
                'end'       , 'e', 1, "numeric",   "End of fragment length range (bp). Default: 10000",
                'stepsize'  , 'p', 1, "numeric",   "Step size (bp). Default: 100",
                'smrtcells' , 'c', 1, "numeric",   "Number of SMRT cells used for dataset. Default: 5",
                'cellcost'  , 'd', 1, "numeric",   "Cost per SMRT cell. Default: 300",
                'maxcov'    , 'm', 1, "numeric",   "Maximum coverage to plot. Default: 10",
                'genomesize', 'g', 1, "numeric",   "Genome size in megabase. Default: 820",
                'ymax'      , 'y', 1, "numeric",   "Maximum number of SMRT cells to plot. Default 100"
               ), ncol=5,byrow=T);
opt  = getopt(args);

# Specify default argument values
if ( is.null(opt$output)     ) { opt$output     = "smrt-coverage.pdf" }
if ( is.null(opt$start)      ) { opt$start      = 200 }
if ( is.null(opt$end)        ) { opt$end        = 10000 }
if ( is.null(opt$stepsize)   ) { opt$stepsize   = 100 }
if ( is.null(opt$smrtcells)  ) { opt$smrtcells  = 5 }
if ( is.null(opt$cellcost)   ) { opt$cellcost   = 300 }
if ( is.null(opt$maxcov)     ) { opt$maxcov     = 10 }
if ( is.null(opt$genomesize) ) { opt$genomesize = 820 }
if ( is.null(opt$ymax)       ) { opt$ymax       = 100 }

# Help message
if ( !is.null(opt$help) || is.null(opt$input)) {
   self=sub("^.*/", "", commandArgs()[4], perl=T)
   cat("\n", getopt(args, command=self, usage=T), "\n")
   q(status=1);
}

########
# MAIN #
########

# Read input
data = read.delim(opt$input)

# Get cumulative coverage plot for input data
x = c();
y = c();
for(i in seq(opt$start, opt$end, opt$stepsize)){
   x = append(x, i);
   y = append(y, sum(as.numeric(data[data[,2]>=i,2])));
}
x = x/1000;
y = (y/1000000);
opt$output = sub("\\.pdf$|\\.PDF$", "", opt$output, perl=T);
pdf(file=paste(opt$output, ".pdf", sep=""))
plot(x,y/opt$genomesize, type="l", main="Cumulative coverage", xlab="Read length threshold (kb)", ylab="Coverage (x-fold)")
legend("topright", legend=paste(opt$smrtcells, "SMRT cells"))

# Get estimates for the number of SMRT cells needed at various coverage levels
coveragelevels = seq(0.5,opt$maxcov, 0.5);
coveragecolors = rainbow(length(coveragelevels));
plot(1,1, type="n", xlim=c(0,max(x)), ylim=c(0,opt$ymax), main="SMRT coverage", xlab="Read length threshold (kb)", ylab="SMRT cell count")
for(i in 1:length(coveragelevels)){
   smrtcount = ((coveragelevels[i]*opt$genomesize)/y)*opt$smrtcells;
   lines(x, smrtcount, col=coveragecolors[i]);
}
legend("bottomright", legend=rev(coveragelevels), col=rev(coveragecolors), lty=1, cex=0.5)
abline(h=46)
dev=dev.off()
