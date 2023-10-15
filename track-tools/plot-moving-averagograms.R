#!/usr/bin/env Rscript

# 30.11.2011 17:09:16 EST
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
args = matrix(c('input'      , 'i', 1, "character", "Comma-delimited list of input files",
                'plotname'   , 'o', 2, "character", "Name of plot output file (required)",
                'plotformat' , 'p', 2, "character", "Plot type (png or pdf). Default: png",
                'bw'         , 'b', 1, "integer",   "Bandwidth for moving window (in bp). Default: 30",
                'spacing'    , 's', 1, "integer",   "Spacing between points (in bp). Default: 10",
                'measure'    , 'm', 1, "character", "Averaging method (mean, geomean or median). Default: mean",
                'type'       , 't', 2, "character", "Drawing type (points, lines, both). default: lines",
                'interval'   , 'v', 0, "logical",   "Plot confidence interval for the mean or IQR for the median",
                'raw'        , 'r', 0, "logical",   "Add the raw data points to the plot",
                'density'    , 'd', 0, "logical",   "Add a density plot of the raw data points using smoothScatter",
                'xlim'       , 'x', 2, "character", "X axis range specified as 'max,min'. Default: automatic",
                'ylim'       , 'y', 2, "character", "Y axis range specified as 'max,min'. Default: automatic",
                'hlines'     , 'k', 2, "character", "Optional comma-separated list of coordinates for horizontal lines",
                'vlines'     , 'l', 2, "character", "Optional comma-separated list of coordinates for vertical lines",
                'invert'     , 'n', 2, "character", "Optional comma separated list of 1's and 0's indicating which lines to invert",
                'help'       , 'h', 0, "logical",   "Brief help message"
               ), ncol=5,byrow=T);
opt  = getopt(args);

# Specify default argument values
if ( is.null(opt$plotformat) ) { opt$plotformat  = "png"   }
if ( is.null(opt$bw)         ) { opt$bw          = 30      }
if ( is.null(opt$spacing)    ) { opt$spacing     = 10      }
if ( is.null(opt$measure)    ) { opt$measure     = "mean"  }
if ( is.null(opt$type)       ) { opt$type        = "lines" }

# Help message
if ( !is.null(opt$help) || is.null(opt$input) || is.null(opt$plotname) ) {
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

movingGeomean <- function(x, bw, xvals, yvals, featurecount, func){
   dn = x-bw;
   up = x+bw;
   subset = (xvals >= dn) & (xvals <= up);
   d = yvals[subset];
   setsize = (((2*bw)+1)*featurecount);
   if (length(d) < setsize){
      d = c(d, rep(0,  (setsize - length(d)) ));
   }
   return(mean(log(d+1), na.rm=T));
}

movingMean <- function(x, bw, xvals, yvals, featurecount, func){
   dn = x-bw;
   up = x+bw;
   subset = (xvals >= dn) & (xvals <= up);
   d = yvals[subset];
   setsize = (((2*bw)+1)*featurecount);
   if (length(d) < setsize){
      d = c(d, rep(0,  (setsize - length(d)) ));
   }
   return(mean(d, na.rm=T));
}

movingMedian <- function(x, bw, xvals, yvals, featurecount, func){
   dn = x-bw;
   up = x+bw;
   subset = (xvals >= dn) & (xvals <= up);
   d = yvals[subset];
   setsize = (((2*bw)+1)*featurecount);
   if (length(d) < setsize){
      d = c(d, rep(0,  (setsize - length(d)) ));
   }
   return(median(d, na.rm=T));
}

movingMeanError <- function(x, bw, xvals, yvals, featurecount) {
   dn = x-bw;
   up = x+bw;
   subset = (xvals >= dn) & (xvals <= up);
   d = yvals[subset];
   setsize = (((2*bw)+1)*featurecount);
   if (length(d) < setsize){
      d = c(d, rep(0,  (setsize - length(d)) ));
   }
   n = length(d)
   if (n>1){
      m = mean(d, na.rm=T);
      s = sd(d, na.rm=T);
      return(qt(0.975,df=n-1)*s/sqrt(n));
   }else{
      return(NA);
   }
}

movingMedianQuantiles <- function(x, bw, xvals, yvals, featurecount, probs) {
   dn = x-bw;
   up = x+bw;
   subset = (xvals >= dn) & (xvals <= up);
   d = yvals[subset];
   setsize = (((2*bw)+1)*featurecount);
   if (length(d) < setsize){
      d = c(d, rep(0,  (setsize - length(d)) ));
   }
   return(ifelse(any(subset), as.numeric(quantile(d,probs=probs, na.rm=T)),NA));
}

readRawData <- function(files=c(), header=T){
   out = list(ymax=NULL, ymin=NULL, xmax=NULL, xmin=NULL, data=NULL);
   for(i in 1:length(files)){
      out$data[[i]] = read.delim(files[i], header=header, check.names=F, stringsAsFactors=F);
      out$featurecount[[i]] = length(unique(out$data[[i]][,1]))
      out$data[[i]] = out$data[[i]][,2:3];
      
      # Find ymax and ymin for dataset
      tmp.ymax       = max(out$data[[i]][,2]);
      tmp.ymin       = min(out$data[[i]][,2]);
      if(is.null(out$ymax)){
         out$ymax = tmp.ymax;
      }
      else{
         out$ymax = ifelse(tmp.ymax > out$ymax, tmp.ymax, out$ymax);
      }
      if(is.null(out$ymin)){
         out$ymin = tmp.ymin;
      }
      else{
         out$ymin = ifelse(tmp.ymin < out$ymin, tmp.ymin, out$ymin);
      }
      
      # Find xmax and xmin for dataset
      tmp.xmax       = max(out$data[[i]][,1]);
      tmp.xmin       = min(out$data[[i]][,1]);
      if(is.null(out$xmax)){
         out$xmax = tmp.xmax;
      }
      else{
         out$xmax = ifelse(tmp.xmax > out$xmax, tmp.xmax, out$xmax);
      }
      if(is.null(out$xmin)){
         out$xmin = tmp.xmin;
      }
      else{
         out$xmin = ifelse(tmp.xmin < out$xmin, tmp.xmin, out$xmin);
      }
   }
   return(out);
}

smoothData <- function(self=list(), bw=30, spacing=10, measure="mean", interval=F, xlim=NULL){
   positions = seq(xlim[1],xlim[2],spacing);
   out       = list(plot.x=positions, plot.y=NULL, plot.yup=NULL, plot.ydn=NULL, xmin=min(positions), xmax=max(positions), ymin=NULL, ymax=NULL);
   for(i in 1:length(self$data)){
      if(measure == "mean"){
         out$plot.y[[i]] = sapply(positions, movingMean, bw=bw, xvals=self$data[[i]][,1], yvals=self$data[[i]][,2], featurecount=self$featurecount[[i]] );
         if (interval){
            plot.error = sapply(positions, movingMeanError, bw=bw, xvals=self$data[[i]][,1], yvals=self$data[[i]][,2], featurecount=self$featurecount[[i]] );
            out$plot.yup[[i]]   = out$plot.y[[i]] - plot.error;
            out$plot.ydn[[i]]   = out$plot.y[[i]] + plot.error;
         }
      }else if (measure == "geomean"){
        out$plot.y[[i]] = sapply(positions, movingGeomean, bw=bw, xvals=self$data[[i]][,1], yvals=self$data[[i]][,2], featurecount=self$featurecount[[i]] );
         if (interval){
            stop(paste("Error: confidence interval plotting not supported for measure '", measure, "'.\n"));
         } 
      }else{
         out$plot.y[[i]]    = sapply(positions, movingMedian, bw=bw, xvals=self$data[[i]][,1], yvals=self$data[[i]][,2], featurecount=self$featurecount[[i]] );
         if (interval){
            out$plot.yup[[i]] = sapply(positions, movingMedianQuantiles, bw=bw, xvals=self$data[[i]][,1], yvals=self$data[[i]][,2], featurecount=self$featurecount[[i]], probs=0.25);
            out$plot.ydn[[i]] = sapply(positions, movingMedianQuantiles, bw=bw, xvals=self$data[[i]][,1], yvals=self$data[[i]][,2], featurecount=self$featurecount[[i]], probs=0.75);
         }
      }
      
      # Exponentiate the geometric means
      if (measure == "geomean"){
         out$plot.y[[i]] = exp(out$plot.y[[i]])-1;
      }

      # Set ylim
      tmp.ymin = min(c(out$plot.y[[i]], out$plot.yup[[i]], out$plot.ydn[[i]]), na.rm=T);
      tmp.ymax = max(c(out$plot.y[[i]], out$plot.yup[[i]], out$plot.ydn[[i]]), na.rm=T);
      if(is.null(out$ymax)){
         out$ymax = tmp.ymax;
      }else{
         out$ymax = ifelse(tmp.ymax > out$ymax, tmp.ymax, out$ymax);
      }
      if(is.null(out$ymin)){
         out$ymin = tmp.ymin;
      }else{
         out$ymin = ifelse(tmp.ymin < out$ymin, tmp.ymin, out$ymin);
      }
   }
   return(out);
}

########
# MAIN #
########

# Check arguments
if(!( opt$type %in% c('lines', 'points', 'both') )) stop("Unknown plot type, must be either 'lines', 'points' or 'both'");

# Read raw datafiles
filenames   = unique(unlist(strsplit(opt$input, split=",")));
data.raw    = readRawData(filenames, header=TRUE);

# Set xlim here since we need it for the smoothing
if (!is.null(opt$xlim)){
   xlim = as.numeric(unlist(strsplit(opt$xlim, split=",")));
   if (length(xlim) != 2) stop("Error: xlim must have two values 'min,max'");
   xlim = sort(xlim);
}else{
   xlim = c(data.raw$xmin, data.raw$xmax);
}

# Smooth the data
interval    = ifelse(is.null(opt$interval), FALSE, TRUE);
data.smooth = smoothData(data.raw, bw=opt$bw, spacing=opt$spacing, measure=opt$measure, interval=interval, xlim=xlim);

# Optionally invert the sign of plotted data
if(!is.null(opt$invert)){
   invert = as.numeric(unlist(strsplit(opt$invert, split=",")));
   if(!all(invert %in% c(1,0))) stop("Error: invert argument must be comma-separated list of 1s and 0s");
   if(length(filenames)!=length(invert)) stop("Error: invert list must match the number of datasets to plot");
   for (i in 1:length(invert)){
      if(invert[i]){ 
         data.raw$data[[i]]    = -1 * data.raw$data[[i]];
         data.smooth$plot.y[[i]] = -data.smooth$plot.y[[i]];
         if(interval){
            tmp.ydn = data.smooth$plot.ydn[[i]];
            data.smooth$plot.ydn[[i]] = -data.smooth$plot.yup[[i]];
            data.smooth$plot.yup[[i]] = -tmp.ydn;
         }
         data.raw$ymin    = min( c(data.raw$data[[i]][,1], data.raw$ymin), na.rm=T);
         data.smooth$ymin = min( c(data.smooth$ymin, data.smooth$plot.y[[i]], data.smooth$plot.yup[[i]], data.smooth$plot.ydn[[i]]), na.rm=T);
      }
   }
}

# Now we can set up ylim
plot.raw     = ifelse(is.null(opt$raw), FALSE, TRUE);
plot.density = ifelse(is.null(opt$density), FALSE, TRUE);
if (!is.null(opt$ylim)){
   ylim = as.numeric(unlist(strsplit(opt$ylim, split=",")));
   if (length(ylim) != 2) stop("Error: ylim must have two values 'min,max'");
   ylim = sort(ylim);
}else{
   if(plot.raw | plot.density){
      ylim=c( min(c(data.raw$ymin,data.smooth$ymin), na.rm=T), max(c(data.raw$ymax,data.smooth$ymax), na.rm=T) );
   }else{
      ylim=c(data.smooth$ymin, data.smooth$ymax);
   }
}

# Set up the plot colors, add transparant colors for the polygons
colors.op = rainbow(length(filenames));
colors.tr = rainbow(length(filenames), alpha=0.3);
if(length(colors.op)==1){
   colors.op = "#000000FF";
   colors.tr = "#00000044";
}

# Set up the plot output files
if (opt$plotformat == "png"){
   plot.name = sub("\\.png|\\.PNG", "", opt$plotname, perl=T);
   png(file=paste(plot.name, ".png", sep=""), width=900, height=500);
   layout(matrix(c(1,2), nrow = 1), widths = c(0.8, 0.2));
   par(cex=1,mar = c(5, 4, 4, 2) + 0.1);
}else{
   plot.name = sub("\\.pdf|\\.PDF", "", opt$plotname, perl=T);
   pdf(file=paste(plot.name, ".pdf", sep=""), width=7, height=4);
   layout(matrix(c(1,2), nrow = 1), widths = c(0.8, 0.2));
   par(cex=0.7, mar = c(5, 4, 4, 2) + 0.1);
}

# Set up the empty plot area with axes
plot(1,1, xlim=xlim, ylim=ylim, type="n", main=plot.name, xlab="Distance", ylab="Score");
abline(h=0, col="grey");
abline(v=0, col="grey");

# Do we need to plot any ablines?
if (!is.null(opt$vlines)){
   vlines = unique(unlist(strsplit(opt$vlines, split=",")));
   for(vline in vlines){
      abline(v=vline, col="grey");
   }
}
if (!is.null(opt$hlines)){
   hlines = unique(unlist(strsplit(opt$hlines, split=",")));
   for(hline in hlines){
      abline(h=hline, col="grey");
   }
}


# Make density plots first to make sure rest is overplotted
if (plot.density){
   if(length(filenames) > 1){
      cat("Warning: more than one data set provided, only the first set was used for the density plot!\n");
   }
   smoothScatter(data.raw$data[[1]], add=T, nbin=1000, transformation = function(x) x^0.35, colramp = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")));
}

# Plot raw data first to make sure the lines & polygons are overplotted
if(plot.raw){
   for(i in 1:length(filenames)){
      points(data.raw$data[[i]], pch=".", col=colors.tr[i]);
   }
}

# Next plot the polygons for confidence intervals
if (interval){
   for(i in 1:length(filenames)){
      # The polygons get all screwed up if there is missing data for the confidence intervals
      # so we need to make sure to only draw polygons for ranges with data points!
      ranges     = c();
      selection  = !(is.na(data.smooth$plot.yup[[i]]) | is.na(data.smooth$plot.ydn[[i]]));
      for(index in 1:length(selection)){
         if ( (index == 1) | (index == length(selection)) ){
            if (selection[index]) ranges = append(ranges,index);
         } else{
            if (selection[index]==T & selection[index-1]==F & selection[index+1]==F ) ranges = append(ranges,c(index,index));
            if (selection[index]==T & selection[index-1]==T & selection[index+1]==F ) ranges = append(ranges,index);
            if (selection[index]==T & selection[index-1]==F & selection[index+1]==T ) ranges = append(ranges,index);
         }
      }
      if(length(ranges) %% 2 != 0) stop("Assertion 'ranges==even' failed in polygon plotting\n");
      ranges = matrix(ranges,ncol=2,byrow=T);
      
      # Plot the ranges
      for(index in 1:nrow(ranges)){
         rs = ranges[index,1];
         re = ranges[index,2];
         if ( re-rs > 0 ){
            pol.x = c(data.smooth$plot.x[rs:re], rev(data.smooth$plot.x[rs:re]));
            pol.y = c(data.smooth$plot.yup[[i]][rs:re],rev(data.smooth$plot.ydn[[i]][rs:re]));
            polygon(pol.x, pol.y, col=colors.tr[i], border=NA);
         }
      }
   }
}

# And finally, plot the smoothing lines
for(i in 1:length(filenames)){
   plot.type = switch(opt$type, lines="l", points="p", both="o");
   lines(data.smooth$plot.x, data.smooth$plot.y[[i]], col=colors.op[i], type=plot.type, pch=20);
}

# Add the legend and plotting info
if (opt$plotformat == "png"){
   par(mar = c(5, 0, 4, 2) + 0.1, cex=0.8);
}else{
   par(mar = c(5, 0, 4, 2) + 0.1, cex=0.5);
}
plot(1,1, type="n", xlim=c(0,10), ylim=c(0,10), axes=F, ann = FALSE);
legend("topleft", legend=filenames, lty=1, col=colors.op);
text(0,0.3, paste("measure: ", opt$measure, "\nbw: ", opt$bw, "\nspacing: ", opt$spacing, sep=""), pos=4);


# Close the plot
dev = dev.off();
