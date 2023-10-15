#!/usr/bin/env Rscript

# 30.07.2013 16:20:53 EDT
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
args = matrix(c('input'  , 'i', 1, "character", "Tab-delimited input file listing Name and Filename of datasets to plot (gene_id, transcript_id, promoter size, promoter count, gene body size, gene body count)",
                'output' , 'o', 2, "character", "Name of output file (optional)",
                'thresh' , 't', 2, "integer",   "Promoter proximal read count threshold for active genes",
                'set'    , 's', 2, "character", "Subset of gene IDs to plot traveling ratio for",
                'log10'  , 'l', 0, "logical",   "Export log10 traveling ratio data",
                'xmin'   , 'a', 2, "numeric",   "X axis minimum. Default: -1",
                'xmax'   , 'b', 2, "numeric",   "X axis maximum. Default: 4",
                'invert' , 'v', 0, "logical",   "Invert ratio calculation",
                'help'   , 'h', 0, "logical",   "Brief help message"
               ), ncol=5,byrow=T);
opt  = getopt(args);

# Specify default argument values
if ( is.null(opt$output) ) { opt$output = gsub(".txt$", "", opt$input) }
if ( is.null(opt$xmin)   ) { opt$xmin   = -1                           }
if ( is.null(opt$xmax)   ) { opt$xmax   = 4                            }
if ( is.null(opt$log10)  ) { opt$log10  = FALSE                        }

# Help message
if ( !is.null(opt$help) || is.null(opt$input) || is.null(opt$thresh)) {
   self=sub("^.*/", "", commandArgs()[4], perl=T)
   cat("\n", getopt(args, command=self, usage=T), "\n")
   q(status=1);
}

########
# MAIN #
########

# Load list of files to plot
file.list = read.delim(opt$input, header=T, comment.char="#", stringsAsFactors=F)

# Check input
if(!all(c("Name","Filename") %in% colnames(file.list))){
   stop("Error: The input file must have columns with header names 'Name' and 'Filename'\n");
}
if ( any(duplicated(file.list$Name)) ){
   stop("Error: Duplicate values found in the 'Name' column\n");
}
if ( nrow(file.list) > 12 ){
   stop("Error: Cannot plot more than 12 datasets\n");
}

# Build list object for plotting
plot.data = list();
for( name in file.list$Name ){
   # Read data
   counts    = read.delim(file.list$Filename[file.list$Name == name], header=F, comment.char="#");
   if (ncol(counts)==8){
      cat("subtracting background\n")
      b = (counts[,8]/counts[,7])*1000
      densities = as.matrix(data.frame(p=((counts[,4]/counts[,3])*1000)-b, g=((counts[,6]/counts[,5])*1000)-b ));
      cat(sum(densities<0),"\n")
   } else {
      densities = as.matrix(data.frame(p=(counts[,4]/counts[,3])*1000, g=(counts[,6]/counts[,5])*1000));
   }
   
   # Filter active genes
   if ( is.null(opt$invert) ){
      set.gdens   = densities[,2]>0 & densities[,1]>=0;
   } else {
      set.gdens   = densities[,1]>0 & densities[,2]>=0;
   }
   set.maxdens = apply(densities, 1, max)>median(densities);
   if ( is.null(opt$invert) ){
      set.pdens   = counts[,4]>=opt$thresh;
   } else {
      set.pdens   = counts[,6]>=opt$thresh;
   }
   set.all     = set.gdens & set.maxdens & set.pdens;
   
   # Optionally select a subset of genes to plot TR for
   if(!is.null(opt$set)){
      select     = read.delim(opt$set, header=F);
      set.select = counts[,1] %in% select[,1]
      set.all    = set.all & set.select;
   }
   
   # calculate traveling ratios for active genes
   if ( is.null(opt$invert) ){
      x    = sort(log10(densities[set.all,1]/densities[set.all,2]));
   } else{
      x    = sort(log10(densities[set.all,2]/densities[set.all,1]));
   }
   y    = (1:length(x))/length(x)
   plot.data[[name]] = list(x=x, y=y);
   
   # Get the closest TR to y=0.4
   ymin    = abs(y-0.4)
   TRcomp  = x[ which( ymin==min(ymin) )[1] ]
   cat(name, "\t", TRcomp, "\n")
   
   # Write a table of log10 TRs
   if (opt$log10){
      if ( is.null(opt$invert) ){
         out = data.frame(geneids=counts[set.all,1], log10TR=log10(densities[set.all,1]/densities[set.all,2]));
      }
      else{
         out = data.frame(geneids=counts[set.all,1], log10TR=log10(densities[set.all,2]/densities[set.all,1]));
      }
      write.table(out, file=paste("log10TR_", name, ".txt", sep=""), row.names=F, col.names=T, quote=F, sep="\t");
   }
}

# Draw plot
colors = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928");
colors = colors[1:nrow(file.list)]
names(colors) = names(plot.data)

pdf(paste(opt$output, "pdf", sep="."));
plot(1,1, type="n", xlim=c(opt$xmin,opt$xmax), ylim=c(0,1), xlab="Traveling ratio", ylab="Fraction of active genes", main=opt$output);
for( name in names(plot.data) ){
   lines( plot.data[[name]]$x,  plot.data[[name]]$y, col=colors[[name]])
}
legend("topleft", legend=names(plot.data), col=colors, lty=1);
dev=dev.off();
