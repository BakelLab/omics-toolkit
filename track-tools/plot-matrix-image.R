#!/usr/bin/env Rscript

# 25.04.2013 13:12:32 EDT
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
args = matrix(c('input'       , 'i', 1, "character", "Input file name (required)",
                'colorscale'  , 's', 1, "character", "Color scale mapping file (required)",
                'plotname'    , 'o', 2, "character", "Name of output file (required)",
                'plotformat'  , 'p', 2, "character", "Plot type (png or pdf). Default: png",
                'drow'        , 'r', 2, "integer",   "First data row. Default: 2",
                'dcol'        , 'c', 2, "integer",   "First data col. Default: 3",
                'theta'       , 't', 2, "integer",   "Smoothing theta (typically between 1-3). Default: no smoothing",
                'zmin'        , 'l', 2, "double",    "Z minimum. Default: auto",
                'zmax'        , 'u', 2, "double",    "Z maximum. Default: auto",
                'nlevel'      , 'n', 2, "integer",   "Number of color levels",
                'help'        , 'h', 0, "logical",   "Brief help message"
               ), ncol=5,byrow=T);
opt  = getopt(args);

# Specify default argument values
if ( is.null(opt$output)     ) { opt$output     = "" }
if ( is.null(opt$plotformat) ) { opt$plotformat = "png"        }
if ( is.null(opt$drow)       ) { opt$drow       = 1  }
if ( is.null(opt$dcol)       ) { opt$dcol       = 3  }

# Help message
if ( !is.null(opt$help) || is.null(opt$input) || is.null(opt$plotname) || is.null(opt$colorscale)) {
   self=sub("^.*/", "", commandArgs()[4], perl=T)
   cat("\n", getopt(args, command=self, usage=T), "\n")
   q(status=1);
}

#############
# LIBRARIES #
#############
library(fields)

########
# MAIN #
########

# Read input data
d   = read.table(opt$input, header=F, stringsAsFactors=F, check.names=F, skip=opt$drow-1, sep="\t");
d   = data.matrix(d[,opt$dcol:ncol(d)]);
d   = t(d)[,nrow(d):1];
col = read.table(opt$colors)

# Optionally smooth the image
if (!is.null(opt$theta)){
   d = image.smooth(d, theta = opt$theta);
   d = d$z;
}

# Set zmin
if (is.null(opt$zmin)){
   zmin = min(d);
} else{
   zmin = opt$zmin;
   d[d<zmin] = zmin;
}

# set zmax
if (is.null(opt$zmax)){
   zmax = max(d);
} else{
   zmax = opt$zmax;
   d[d>zmax] = zmax;
}

# Set up the plot output files
if (opt$plotformat == "png"){
   plot.name = sub("\\.png|\\.PNG", "", opt$plotname, perl=T);
   png(file=paste(plot.name, ".png", sep=""), width=900, height=900);
}else{
   plot.name = sub("\\.pdf|\\.PDF", "", opt$plotname, perl=T);
   pdf(file=paste(plot.name, ".pdf", sep=""), width=7, height=4);
}

# Generate the image
image.plot(d,col=rgb(col), zlim=c(zmin,zmax));
dev=dev.off()
