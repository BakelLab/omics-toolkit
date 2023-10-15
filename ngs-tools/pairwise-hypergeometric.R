#!/usr/bin/env Rscript

# 26.02.2018 08:38:47 EST
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
args = matrix(c('input'      , 'i', 1, "character", "Input file name",
                'output'     , 'o', 2, "character", "Name of output file (optional)",
                'universe'   , 'u', 2, "numeric",   "Size of the gene universe",
                'correction' , 'c', 2, "character", "Multiple testing correction approach. Default: Bonferroni",
                'help'       , 'h', 0, "logical",   "Brief help message"
               ), ncol=5,byrow=T);
opt  = getopt(args);

# Specify default argument values
if ( is.null(opt$correction) ) { opt$correction = "bonferroni" }


# Help message
if ( !is.null(opt$help) || is.null(opt$input) || is.null(opt$universe) ) {
   self=sub("^.*/", "", commandArgs()[4], perl=T)
   cat("\n", getopt(args, command=self, usage=T), "\n")
   q(status=1);
}

#############
# LIBRARIES #
#############


#############
# FUNCTIONS #
#############


########
# MAIN #
########

# Set output name
if ( is.null(opt$output) ) { 
   opt$output     = sub('\\.txt$', '', basename(opt$input))
}

# Read input matrix
d   = read.delim(opt$input, row.names=1)

# Make output matrix of the same dimension
out = matrix( nrow=nrow(d), ncol=ncol(d) )
rownames(out) = rownames(d)
colnames(out) = colnames(d)

# Calculate hypergeometric p-values
for (row in (1:nrow(d)) ){
   for (col in (1:ncol(d)) ){
      if( row != col){
         row.list.size = d[row,row];
         col.list.size = d[col,col];
         overlap       = d[row,col];
         if (row > col){
            out[row,col] = phyper(overlap-1, row.list.size, opt$universe-row.list.size, col.list.size, lower.tail = FALSE, log.p = FALSE);
         } else {
            out[row,col] = phyper(overlap-1, col.list.size, opt$universe-col.list.size, row.list.size, lower.tail = FALSE, log.p = FALSE);
         }
      }
   }
}

# Do correction if requested
if ( !is.null(opt$correction) ){
   if ( tolower(opt$correction) == "bonferroni"){
      out = sum(lower.tri(d)) * out
      out[out > 1] = 1
   }
}

# Write output
write.table(out,         file=paste(opt$output, "_hypergeo_pvals.txt", sep=""),       row.names=T, col.names=T, sep="\t", quote=F )
write.table(-log10(out), file=paste(opt$output, "_hypergeo_nlog10pvals.txt", sep=""), row.names=T, col.names=T, sep="\t", quote=F )

