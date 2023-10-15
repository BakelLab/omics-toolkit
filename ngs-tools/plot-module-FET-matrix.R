#!/usr/bin/env Rscript

# 21.02.2018 12:33:46 EST
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
args = matrix(c('input'      , 'i', 1, "character", "gProfiler output file name",
                'output'     , 'o', 2, "character", "Name of output file (optional)",
                'pval'       , 'p', 2, "numeric",   "P-value threshold. Default: 0.05",
                'topN'       , 't', 2, "numeric",   "Restrict to the top N categories per cluster. Default: 3",
                'maxlogp'    , 'm', 2, "numeric",   "Set upper threshold to -logP score. Default: none",
                'summary'    , 's', 2, "character", "Optional module summary file to get module path strings.",
                'modulelist' , 'l', 2, "character", "Optional module list to print matrix for",
                'setlist'    , 'e', 2, "character", "Optional set list to print matrix for",
                'color'      , 'c', 2, "character", "Plot color. Default: darkgreen",
                'width'      , 'w', 2, "numeric",   "Plot width in inches. Default: 8",
                'height'     , 'u', 2, "numeric",   "Plot height in inches. Default: 11",
                'help'       , 'h', 0, "logical"  , "Brief help message"
               ), ncol=5,byrow=T);
opt  = getopt(args);

# Specify default argument values
if ( is.null(opt$output) ) { opt$output = "FET-matrix" }
if ( is.null(opt$domain) ) { opt$pval   = 0.05         }
if ( is.null(opt$topN) )   { opt$topN   = 3            }
if ( is.null(opt$color) )  { opt$color  = "darkgreen"  }
if ( is.null(opt$width) )  { opt$width  = 8            }
if ( is.null(opt$height) ) { opt$height = 11           }

# Help message
if ( !is.null(opt$help) || is.null(opt$input)) {
   self=sub("^.*/", "", commandArgs()[4], perl=T)
   cat("\n", getopt(args, command=self, usage=T), "\n")
   q(status=1);
}


#############
# LIBRARIES #
#############

library(ggplot2)
library(reshape2)
library(dplyr,quietly=T);
library(NMF,quietly=T);
library(heatmaply,quietly=T);
library(htmlwidgets,quietly=T);
library(data.tree,quietly=T);

########
# MAIN #
########

# opt = list(input="DEG_FET-Table.txt", pval=0.05, output="test", topN=3, summary="module_summary.txt", modulelist="module-list-final.txt", setlist="DEG-lists-of-interest.txt")

# Read file and check for required header info
d = read.delim(opt$input, header=T, stringsAsFactors=F)
if ( !all(c("FET_pvalue", "corrected.FET.pvalue", "set1_Name", "set2_Name", "intersects") %in% colnames(d)) )
   stop("Error: required column names not found\n");

# Subset to pval thresh
d = d[ d$corrected.FET.pvalue < opt$pval ,]

# Optionally, read list of modules to include in the matrix
module.list = NULL;
if ( !is.null(opt$modulelist) ){
   module.list = as.character(read.delim(opt$modulelist, header=F)[,1])
}

# Optionally, read module summary file and remap cluster IDs to full pathstring hierarchy
if ( !is.null(opt$summary) ){
   module.summary = read.delim(opt$summary, header=T)
   module.pc      = module.summary[,c("module.id", "module.parent")]
   module.tree    = FromDataFrameNetwork(module.pc)
   module.pstring = ToDataFrameTree(module.tree, "pathString")[,2]
   pstring.map    = data.frame(module=gsub(".*\\/","",module.pstring), pathstring=gsub("comp1_"," ",module.pstring))

   
   pstring.lut    = match(as.character(d$set1_Name), pstring.map$module)
   if ( any( is.na(pstring.lut) ) ){
      stop("Error: modules could not be remapped because not all module names are represented in the module_summary file\n");
   } else {
      d$set1_Name = as.character(pstring.map$pathstring[pstring.lut])
   }
   
   if ( !is.null(module.list) ){
      pstring.lut    = match(as.character(module.list), pstring.map$module)
      if ( any( is.na(pstring.lut) ) ){
         stop("Error: modulelist could not be remapped because not all module names are represented in the module_summary file\n");
      } else {
         module.list = as.character(pstring.map$pathstring[pstring.lut])
      }
   }
   
}

# Get all module names
if ( is.null(module.list) ){
   out.modules = as.character( unique(d$set1_Name) );
} else {
   out.modules = module.list;
   d = d[ d$set1_Name %in% out.modules, ]
}

# Select only the top N terms per cluster
if ( !is.null(opt$topN) ){
   d.subset = NULL;
   for ( i in out.modules ){
      sub = d[ d$set1_Name == i, ];
      sub = sub[ order(sub$FET_pvalue, decreasing=F), ]
      if ( nrow(sub) > opt$topN ){
         sub = sub[ 1:opt$topN, ]
      }
      if (is.null(d.subset)){
         d.subset = sub;
      } else {
         d.subset = rbind(d.subset, sub)
      }
   }
} else {
   d.subset = d;
}


# Get all unique categories
if ( is.null(opt$setlist) ){
   out.terms  = as.character( unique(d.subset$set2_Name) )
} else {
   out.terms  = as.character(read.delim(opt$setlist, header=F)[,1])
   d = d[ d$set2_Name %in% out.terms, ]
}

# Set up output matrix
out.matrix       = matrix( nrow=length(out.modules), ncol=length(out.terms), data=0)
rownames(out.matrix) = out.modules
colnames(out.matrix) = out.terms

# Cycle through terms to fill the matrix
for ( i in 1:ncol(out.matrix) ){
   term         = colnames(out.matrix)[i]
   term.modules = as.character( d$set1_Name[ d$set2_Name == term ] )
   term.pvals   = -log10( d$corrected.FET.pvalue[ d$set2_Name == term ] )
   
   if ( any(row.names(out.matrix) %in% term.modules) ){
      out.matrix[ row.names(out.matrix) %in% term.modules, i ] = term.pvals;
   }
}

# Write output matrix
write.table(out.matrix, file=paste(opt$output, "_matrix.txt", sep=""), row.names=T, col.names=T, quote=F, sep="\t")

# Write output table
write.table(d, file=paste(opt$output, "_filtered.txt", sep=""), row.names=T, col.names=T, quote=F, sep="\t")
write.table(d.subset, file=paste(opt$output, "_filtered_top", opt$topN ,".txt", sep=""), row.names=F, col.names=T, quote=F, sep="\t")

# Restrict max score for plotting if requested
if ( !is.null(opt$maxlogp) ){
   out.matrix[out.matrix > opt$maxlogp] = opt$maxlogp
}

# Image plot the output matrix
out.melt = melt(out.matrix[nrow(out.matrix):1,])
out.melt = out.melt[order(as.character(out.melt$Var1)),]
out.melt$Var1 = as.character(out.melt$Var1)
out.melt$Var2 = factor(out.melt$Var2, levels=out.terms)
out.plot = ggplot(data = out.melt, aes(x=Var2, y=Var1, fill=value)) +
               geom_tile(color="grey") +
               scale_fill_gradient2(low = "white", high = opt$color, mid="white") +
               xlab( paste("Set", opt$domain, "category")) +
               ylab("Module") + 
               theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0, size=5), axis.text.y=element_text(size=5))

pdf( paste(opt$output, "pdf", sep="."), width=opt$width, height=opt$height)
print(out.plot)
dev.off()

# Make a pdf heatmap
nmf.options(grid.patch=TRUE) # Patch to avoid generating a blank page before heatmap
pdf( paste(opt$output, "_heatmap.pdf", sep=""), width=opt$width, height=opt$height )
aheatmap(t(out.matrix), Colv=T, Rowv=T, distfun=function(x) dist(x, method="euclidian") );
nmf.options(grid.patch=FALSE);
dev.off()

# Make an interactive heatmap
heatmap = heatmaply(t(out.matrix), Rowv=T, Colv=T, distfun=function(x) dist(x, method="euclidian"), margins = c(120, 250) );
saveWidget(heatmap, paste(opt$output, "_heatmap.html", sep=""), selfcontained=F);
unlink("Rplot001.jpeg");

