#!/usr/bin/env Rscript

# 16.05.2013 12:56:20 EDT
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
args = matrix(c('input'  , 'i', 1, "character", "Tab-delimited input file (name, position, coverage)",
                'output' , 'o', 2, "character", "Plot output file prefix",
                'lowcov' , 'l', 2, "numeric",   "Highlight areas with coverage below this value",
                'rows'   , 'r', 2, "numeric",   "Number of plot rows per page. Default: 2",
                'cols'   , 'c', 2, "numeric",   "Number of plot columns per page. Default: 1",
                'bars'   , 'b', 2, "logical",   "Plot bars instead of a polygon (useful for sparse data).",
                'ylim'   , 'y', 2, "character", "Y axis range specified as 'max,min'. Default: automatic",
                'help'   , 'h', 0, "logical",   "Brief help message"
               ), ncol=5,byrow=T);
opt  = getopt(args);

# Specify default argument values
if ( is.null(opt$lowcov) ) { opt$lowcov = 0  }
if ( is.null(opt$rows)   ) { opt$rows   = 2  }
if ( is.null(opt$cols)   ) { opt$cols   = 1  }

# Help message
if ( !is.null(opt$help) || is.null(opt$input) || is.null(opt$output)) {
   self=sub("^.*/", "", commandArgs()[4], perl=T)
   cat("\n", getopt(args, command=self, usage=T), "\n")
   q(status=1);
}

########
# MAIN #
########

data = read.table(opt$input, colClasses=c("character","integer","numeric"), sep="\t", comment.char="#", header=F)
#median.cov = tapply(data[,3],data[,1], median);
#seglist = names(median.cov)[order(median.cov, decreasing=T)]
seglist = sort(unique(as.character(data[,1])));
pdf(file=paste(opt$output, "pdf", sep="."));
par(mfrow=c(opt$rows,opt$cols));
for(seg in seglist){
   set  = data[,1]==seg;
   x    = data[set,2];
   y    = data[set,3];
   if (!is.null(opt$ylim)){
      ylim = as.numeric(unlist(strsplit(opt$ylim, split=",")));
      if (length(ylim) != 2) stop("Error: ylim must have two values 'min,max'");
      ylim = sort(ylim);
   }else{
      ylim = c(0, max(y));
   }
   plot(x,y, type="n",xlab="Position (bp)", ylab="Coverage", ylim=ylim, main=seg);
   if(opt$lowcov){
      lowcov = x[y<opt$lowcov]
      lines(lowcov,rep(ylim[2],length(lowcov)),type="h", col="red");
   }
   if (is.null(opt$bars)){
      polygon(c(min(x), x, max(x)), c(0, y, 0),  col = "blue");
   }
   else{
      lines(x,y,type="h");
   }
}
dev=dev.off()
