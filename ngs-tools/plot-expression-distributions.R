#!/usr/bin/env Rscript

# 24.11.2018 17:37:00 EST
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
args = matrix(c('input'  , 'i', 1, "character", "Input file name with expression levels",
                'anno'   , 'a', 1, "character", "Annotations in anno format",
                'list'   , 'l', 1, "character", "List of gene IDs",
                'groups' , 'g', 1, "character", "Comma separated list of annotation groups to plot. Default: protein_coding,lincRNA,lncRNA_RefSeq,lncRNA_Cabili",
                'output' , 'o', 2, "character", "Name of output file (optional)",
                'help'   , 'h', 0, "logical",   "Brief help message"
               ), ncol=5,byrow=T);
opt  = getopt(args);

# Specify default argument values
if ( is.null(opt$output) ) { opt$output = "plot.pdf" }
if ( is.null(opt$groups) ) { opt$groups = "protein_coding,lincRNA,lncRNA_RefSeq,lncRNA_Cabili" }

# Help message
if ( !is.null(opt$help) || is.null(opt$input) || is.null(opt$anno) || is.null(opt$list)) {
   self=sub("^.*/", "", commandArgs()[4], perl=T)
   cat("\n", getopt(args, command=self, usage=T), "\n")
   q(status=1);
}

#############
# LIBRARIES #
#############

library(ggplot2)
library(reshape)
library(plyr)
library(ggrepel)

#############
# FUNCTIONS #
#############


########
# MAIN #
########

#opt = list(input="logCPM-normalized.txt", anno="gencode.v19.cabili.refseq.annotation.ERCC.anno", list="lncRNA-IDs.txt");

# Read data
d.count = read.delim(opt$input, stringsAsFactor=F, check.names=F)
d.anno  = read.delim(opt$anno, stringsAsFactor=F, check.names=F)
d.list  = read.delim(opt$list, stringsAsFactor=F, check.names=F)

# Check annotation file format
if( !all( c("GeneID","Chr","Length","Symbol","type_of_gene") %in% colnames(d.anno)) ){
   stop("Annotation file must be tab-delimited and contain at least the following column headers: GeneID Chr Length Symbol type_of_gene\n");
}

# Allow for extra column with info on list IDs
if (ncol(d.list)==2){
   colnames(d.list) = c("ID","List_type");
}

# Split annotation sets
opt$groups = unlist(strsplit(opt$groups, split=","));

# Make plots for every column in the data file
pdf(opt$output)
for (i in 1:ncol(d.count) ){
   d.plot = d.count[,i, drop=F]
   plot.name = colnames(d.count)[i]
   
   # Get values for d.list
   list.dat = d.plot[rownames(d.plot) %in% d.list[,1],, drop=F];
   list.dat = data.frame(names = rownames(list.dat), values=list.dat[,1]);
   list.dat = merge(list.dat, d.anno, by.x="names", by.y="GeneID", all.x=T)
   if (ncol(d.list)== 2){
      list.dat = merge(list.dat, d.list, by.x="names", by.y="ID", all.x=T)
      list.dat$plotid = paste(list.dat$names, list.dat$List_type, sep=" - ")
   } else {
      list.dat$plotid = paste(list.dat$names, list.dat$type_of_gene, sep=" - ")
   }
   
   # Gather plot data
   plot.dat = NULL;
   for (group in opt$groups){
      group.ids = d.anno$GeneID[d.anno$type_of_gene == group];
      group.dat = d.plot[rownames(d.plot) %in% group.ids,];
      if (is.null(plot.dat)){
         plot.dat = data.frame(group = rep(group,length(group.dat)), values=group.dat);
      } else {
         plot.dat = rbind(plot.dat, data.frame(group = rep(group,length(group.dat)), values=group.dat) );
      }
   }
   
   # Get median and upper/lower quartiles
   mu = ddply(plot.dat, "group", summarise, grp.md=median(values), grp.dn=quantile(values, probs = 0.15), grp.up=quantile(values, probs = 0.85))
   
   # Build ggplot object
   p = ggplot(plot.dat, aes(x=values, fill=group)) +
         geom_density(alpha=0.2) +
         geom_vline(data=mu, aes(xintercept=grp.dn, color=group), linetype="dotted") +
         geom_vline(data=mu, aes(xintercept=grp.md, color=group), linetype="solid") +
         geom_vline(data=mu, aes(xintercept=grp.up, color=group), linetype="dashed")
   
   # Get yrange of plot
   ylim = layer_scales(p)$y$range$range  
   
   # Add list points with labels
   p = p +
         geom_point(data=list.dat, aes(x=values, y=0)) +
         geom_label_repel(data=list.dat, aes(x=values, y=0, label = plotid), box.padding = 0.5, point.padding = 0.5, segment.color = 'grey50', size=2, ylim=ylim, direction="y", fill="white") +
         ggtitle(plot.name) +
         theme_classic()
   print(p)
}
graphics.off()

