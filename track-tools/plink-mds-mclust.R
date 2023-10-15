#!/home/dpinto/bin/Rscript

# 18.01.2012 13:08:27 EST
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
args = matrix(c('input'      , 'i', 1, "character", "Input file name (.mds plink output file)",
                'output'     , 'o', 2, "character", "Name of output file (optional)",
                'maxclust'   , 'm', 2, "integer",   "Maximum number of clusters to consider. Default: 70",
                'components' , 'c', 2, "integer",   "No. of principal components to consider. Default: from plink file",
                'help'       , 'h', 0, "logical",   "Brief help message"
               ), ncol=5,byrow=T);
opt  = getopt(args);

# Specify default argument values
if ( is.null(opt$output)   ) { opt$output   = sub("\\..*$", "", opt$input, perl=T) }
if ( is.null(opt$maxclust) ) { opt$maxclust = 70 }

# Help message
if ( !is.null(opt$help) || is.null(opt$input)) {
   self=sub("^.*/", "", commandArgs()[4], perl=T)
   cat("\n", getopt(args, command=self, usage=T), "\n")
   q(status=1);
}

#############
# LIBRARIES #
#############
library(mclust)


########
# MAIN #
########

# Read data
MDS = read.table(opt$input,header=T);

# Use mclust to cluster data based on the number of principal components in the .mds file from plink
if ( is.null(opt$components) ) { opt$components = ncol(MDS)-3; }
MCLUST   = Mclust(MDS[, 4:(opt$components+3)], G=1:opt$maxclust)
write.table(cbind(as.character(MDS[,2]),MCLUST$classification), file=paste(opt$output, ".mclust.classify.txt", sep=""), quote=F, row.names=F, col.names=F, sep="\t")

# Plot 1st and 2nd principal components
pdf(file=paste(opt$output, ".mclust.pc12.pdf", sep=""), width=7, height=7)
plot(MDS[,4],MDS[,5],xlab="Coordinate on First Dimension",ylab="Coordinate on Second Dimension",main="MDS analysis", col=MCLUST$classification)
legend("bottomright", legend=unique(MCLUST$classification), col=unique(MCLUST$classification), pch=1);
dev=dev.off()

# Plot all other components against each other
pdf(file=paste(opt$output, ".mclust.pcALL.pdf", sep=""), width=7, height=7)
plot(MDS[,4:(opt$components+3)], col=MCLUST$classification, cex=0.8)
dev=dev.off()
