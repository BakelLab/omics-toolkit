#!/usr/bin/env Rscript

# 16.02.2017 09:54:53 EST
# Harm van Bakel <hvbakel@gmail.com>
# Ana S. Gonzalez-Reiche

#############
# ARGUMENTS #
#############

#  Col 1: the long flag name. A multi-character string.
#  Col 2: short flag alias of Column 1. A single-character string.
#  Col 3: 0=no argument, 1=required argument, 2=optional argument.
#  Col 4: Possible values: logical, integer, double, complex, character.
#  Col 5: A brief description of the purpose of the option.
library(getopt)
args = matrix(c('input'    , 'i', 1, "character", "Tab-delimited file with heatmap data",
                'genelist' , 'g', 2, "character", "Name of file containing list of genes to plot (optional)",
                'colanno'  , 'a', 2, "character", "Name of file with column annotations. Must contain header. (optional)",
                'output'   , 'o', 2, "character", "Name of output file (optional)",
                'cluster'  , 'c', 2, "character", "Indicate row clustering: euclidian, pearson, none. Default: pearson",
                'scale'    , 's', 2, "character", "Indicate scaling: none, row, column. Default: row",
                'zmax'     , 'u', 2, "numeric",   "Maximum color scale value. Default: autoscale to 99th percentile",
                'zmin'     , 'l', 2, "numeric",   "Minimum color scale value. Default: autoscale to 99th percentile",
                'psi'      , 'p', 0, "logical",   "When dealing with PSI data, select maximum intron per cluster to plot heatmap",
                'help'     , 'h', 0, "logical",   "Brief help message"
               ), ncol=5,byrow=T);
opt  = getopt(args);

# Specify default argument values
if ( is.null(opt$output) ) { opt$output  = "plot"    }
if ( is.null(opt$cluster)) { opt$cluster = "pearson" }
if ( is.null(opt$scale)  ) { opt$scale   = "row"     }

# Help message
if ( !is.null(opt$help) || is.null(opt$input)) {
   self=sub("^.*/", "", commandArgs()[4], perl=T)
   cat("\n", getopt(args, command=self, usage=T), "\n")
   q(status=1);
}

#############
# LIBRARIES #
#############

suppressMessages(library(dplyr,quietly=T));
suppressMessages(library(NMF,quietly=T));
suppressMessages(library(heatmaply,quietly=T));
suppressMessages(library(htmlwidgets,quietly=T));

########
# MAIN #
########

tryCatch( {

# Don't treat strings as factors by default
options(stringsAsFactors = FALSE)

# Check clustering options
opt$cluster = tolower(opt$cluster);
if(!all(opt$cluster %in% c("none","euclidian","pearson") )){
   stop("Unknown cluster option selected. Choose from none, euclidian, or pearson\n");
}

# Check clustering options
opt$scale = tolower(opt$scale);
if(!all(opt$scale %in% c("none","row","column") )){
   stop("Unknown scale option selected. Choose from none, row, or column\n");
}

# Read matrix data
if ( is.null(opt$psi) ){
   d = read.delim(opt$input, check.names=F);
} else {
   d            = read.delim(opt$input, check.names=F, row.names=NULL);
   leafclusters = d$row.names
   d            = d[,-1]
}

# Read column annotations if provided
if (!is.null(opt$colanno)) {
   annCol = read.delim(opt$colanno, check.names=F);
   
   # Check if column annotations match and then remove the sample name column
   if(!all(colnames(d) == annCol[,1])){
      stop("First column in the column annotation matrix does not match the data file header\n");
   } else {
      annCol = annCol[,2:ncol(annCol)];
   }
   
} else {
   annCol = NULL;
}

# When dealing with leafcutter PSI data, select the intron with the largest delta-PSI between sample groups
if ( !is.null(opt$psi) ){
   groups = unique(annCol[,1]);
   if ( length(groups) == 2){
      colset.A   = annCol[,1]==groups[1];
      colset.B   = annCol[,1]==groups[2];
      mean.psi.A = apply( d[, colset.A], 1, mean, na.rm=T);
      mean.psi.B = apply( d[, colset.B], 1, mean, na.rm=T);
      delta.psi  = abs(mean.psi.A - mean.psi.B);
      
      dat.psi = data.frame(Clusters=leafclusters, delta.psi=delta.psi, row.id=seq(1,length(leafclusters))); 
      
      max.psi = sort(unlist(lapply( split(dat.psi, dat.psi$Clusters), function(y) { rowid = which ( y$delta.psi == max(y$delta.psi) )[1]; return( y$row.id[rowid] ) } )));
      
      # Check if nothing went wrong
      if ( !all(names(max.psi) == leafclusters[max.psi]) ){
         stop("Error during selection of maximum delta-psi intron among leafcutter clusters\n");
      }
      #d = d[max.psi,]
      #rownames(d) = names(max.psi)
      
   } else{
      stop("Only two sample groups are supported for leafcutter PSI data\n");
   }
}

# Read gene list if provided and use it to subset the data matrix
if (!is.null(opt$genelist)) {
   genes.list = read.delim(opt$genelist, header=F);
   genes.list = as.character(genes.list[,1]);
   d = d[ rownames(d) %in% genes.list, ];
   if ( nrow(d) < 1){
      stop("No data remains after gene filtering\n");
   }
}

# Get a scaled version of the matrix
if ( opt$scale == "column" ){
   d.scale = scale(d);
} else if ( opt$scale == "row" ){
   d.scale = t(scale(t(d), center=T, scale=F));
} else{
   d.scale = d;
}
# remove rows with no values (NA) after scaling to avoid errors when plotting
if (any(is.na(d.scale))){
  na.row = apply(d.scale, 1, function(x){all(is.na(x))})
  d.scale = d.scale[!na.row,]
}

# Set up color pallette
if (!is.null(opt$zmin)){
   zmin = opt$zmin;
} else {
   zmin = -max(abs(d.scale), na.rm=T);
}
if (!is.null(opt$zmax)){
   zmax = opt$zmax;
} else {
   zmax = max(abs(d.scale), na.rm=T);
}
breaks     = seq(zmin, zmax, length.out=50);
my_palette = colorRampPalette(c("blue", "white", "red"))(n = length(breaks) );


# Make heatmap plot
nmf.options(grid.patch=TRUE) # Patch to avoid generating a blank page before heatmap
tryCatch(
   {
      png(file=paste(opt$output, "_heatmap.png", sep = ""), width=1000, height=1000);
      if (opt$cluster == "pearson") {
         aheatmap(d.scale, cexRow=0, cexCol=0, Colv=T, Rowv=T, annCol=annCol, legend=F, annLegend=F, color=my_palette, breaks=breaks, distfun=function(x) as.dist((1-cor(t(x),method="pearson", use="complete.obs"))/2) );
      } else if (opt$cluster == "euclidian") {
         aheatmap(d.scale, Colv=T, Rowv=T, annCol=annCol, color=my_palette, breaks=breaks, distfun=function(x) dist(x, method="euclidian") );
      } else {
         aheatmap(d.scale, Colv=T, Rowv=NA, annCol=annCol, color=my_palette, breaks=breaks );
      }
   },
   error=function(e){ # Just print exception message
      message("Error generating heatmap: ");
      print(e);
   },
   finally=dev.off()
);
nmf.options(grid.patch=FALSE);

# Make an interactive heatmap plot
outfile   = paste(opt$output,"_heatmap.html", sep="");

if (opt$cluster == "pearson") {
    heatmap = heatmaply(d.scale, Rowv=T, Colv=T, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red", midpoint = 0 ), col_side_colors=annCol, distfun=function(x) as.dist((1-cor(t(x),method="pearson"))/2), margins = c(120, 120), limits=c(zmin,zmax) );
} else if (opt$cluster == "euclidian") {
    heatmap = heatmaply(d.scale, Rowv=T, Colv=T, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red", midpoint = 0 ), col_side_colors=annCol, distfun=function(x) dist(x, method="euclidian"), margins = c(120, 120), limits=c(zmin,zmax) );
} else {
    heatmap = heatmaply(d.scale, Rowv=F, Colv=F, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red", midpoint = 0 ), col_side_colors=annCol, margins = c(120, 120), limits=c(zmin,zmax) );
}
saveWidget(heatmap, outfile, selfcontained=F);
unlink("Rplot001.jpeg");

# Make violin plots per gene
if( nrow(d) > 0 ){

   filename = paste(opt$output, "_barplot-SE.pdf", sep="");
   pdf(file=filename);
   gene.groups = annCol[,1];
   for (gene in 1:nrow(d)){
      gene.data     = d[gene,];
      gene.plotdata = data.frame(PSI=as.numeric(gene.data), Group=gene.groups);
      
      # Fix control position
      unique.groups  = unique(gene.plotdata$Group)
      ctl.pos        = match("CTL", unique.groups)
      if ( !is.na(ctl.pos) ){
         unique.groups = c(unique.groups[-ctl.pos], unique.groups[ctl.pos])
      }
      gene.plotdata$Group = factor(gene.plotdata$Group, levels=unique.groups)
      
      gene.summary  = gene.plotdata %>%
         group_by(Group) %>%
         summarize(PSI_mean = mean(PSI),
                   PSI_se = sqrt(var(PSI)/length(PSI)));
      gg = ggplot(gene.summary) + aes(x = Group, y = PSI_mean, fill = Group) + geom_bar(stat="identity") +
         geom_errorbar(aes(x = Group, ymin = PSI_mean-PSI_se, ymax = PSI_mean+PSI_se), color = "black", width = 0.2, data = gene.summary)+
         ggtitle(rownames(d)[gene]) +
         theme(legend.position = "none")
      if (!is.null(opt$psi)){
         gg = gg + ylim(0,1)
      }
      print(gg);
   }
   dev.off();

   filename = paste(opt$output, "_barplot-SD.pdf", sep="");
   pdf(file=filename);
   gene.groups = annCol[,1];
   for (gene in 1:nrow(d)){
      gene.data     = d[gene,];
      gene.plotdata = data.frame(PSI=as.numeric(gene.data), Group=gene.groups);
      
      # Fix control position
      unique.groups  = unique(gene.plotdata$Group)
      ctl.pos        = match("CTL", unique.groups)
      if ( !is.na(ctl.pos) ){
         unique.groups = c(unique.groups[-ctl.pos], unique.groups[ctl.pos])
      }
      gene.plotdata$Group = factor(gene.plotdata$Group, levels=unique.groups)
      
      gene.summary  = gene.plotdata %>%
         group_by(Group) %>%
         summarize(PSI_mean = mean(PSI),
                   PSI_se = sd(PSI));
      gg = ggplot(gene.summary) + aes(x = Group, y = PSI_mean, fill = Group) + geom_bar(stat="identity") +
         geom_errorbar(aes(x = Group, ymin = PSI_mean-PSI_se, ymax = PSI_mean+PSI_se), color = "black", width = 0.2, data = gene.summary)+
         ggtitle(rownames(d)[gene]) +
         theme(legend.position = "none")
      if (!is.null(opt$psi)){
         gg = gg + ylim(0,1)
      }
      print(gg);
   }
   dev.off();

   filename = paste(opt$output, "_violinplot.pdf", sep="");
   pdf(file=filename);
   gene.groups = annCol[,1];
   for (gene in 1:nrow(d)){
      gene.data     = d[gene,];
      gene.plotdata = data.frame(PSI=as.numeric(gene.data), Group=gene.groups);
      
      # Fix control position
      unique.groups  = unique(gene.plotdata$Group)
      ctl.pos        = match("CTL", unique.groups)
      if ( !is.na(ctl.pos) ){
         unique.groups = c(unique.groups[-ctl.pos], unique.groups[ctl.pos])
      }
      gene.plotdata$Group = factor(gene.plotdata$Group, levels=unique.groups)
      
      gene.summary  = gene.plotdata %>%
         group_by(Group) %>%
         summarize(PSI_mean = mean(PSI),
                   PSI_se = sqrt(var(PSI)/length(PSI)));
      gg = ggplot() + 
         geom_violin(aes(x = Group, y = PSI, fill = Group), data = gene.plotdata, , alpha = 0.1) + 
         geom_point(aes(x = Group, y = PSI, color = Group), position = position_jitter(width = 0.1, height = 0.0), alpha = 0.6, data = gene.plotdata) +
         geom_point(aes(x = Group, y = PSI_mean), color = "black", size = 2, data = gene.summary) + 
         geom_errorbar(aes(x = Group, ymin = PSI_mean-PSI_se, ymax = PSI_mean+PSI_se), color = "black", width = 0.2, data = gene.summary)+
         ggtitle(rownames(d)[gene]) +
         theme(legend.position = "none")
      if (!is.null(opt$psi)){
         gg = gg + ylim(0,1)
      }
      print(gg);
   }
   dev.off();
   
} else {
   cat( "## Warning: no highlight genes found for individual gene plots\n#\n", sep="");
}

## Save environment for debugging when encountering any errors
},
error=function(e){ # Just print exception message
   message("\n### An error occurred during plotting:\n\n");
   print(e);
   save.image(file = "error.RData");
   }
);
