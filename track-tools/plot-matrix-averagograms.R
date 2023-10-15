#!/usr/bin/env Rscript

# 21.01.2010 16:48:48 EST
# Harm van Bakel <hvbakel@gmail.com>

#############
# LIBRARIES #
#############
library(getopt);
library(splus2R);

#############
# ARGUMENTS #
#############

#  Col 1: the long flag name. A multi-character string.
#  Col 2: short flag alias of Column 1. A single-character string.
#  Col 3: 0=no argument, 1=required argument, 2=optional argument.
#  Col 4: Possible values: logical, integer, double, complex, character.
#  Col 5: A brief description of the purpose of the option.
args = matrix(c('input'       , 'i', 1, "character", "Comma separated list of standard feature matrix files to plot",
                'plotname'    , 'o', 2, "character", "Name of plot output file (required)",
                'plotformat'  , 'p', 2, "character", "Plot type (png or pdf). Default: png",
                'measure'     , 'm', 2, "character", "Measure to plot (mean, geomean, median or sum). Default: mean",
                'bw'          , 'b', 2, "integer",   "Bandwith in number of bins for smoothing. Default: 1 (no smoothing)",
                'interval'    , 'v', 0, "logical",   "Plot confidence interval for the mean or IQR for the median",
                'fit'         , 'f', 2, "character", "Drawing type (points, lines or spline). default: lines",
                'xlim'        , 'x', 2, "character", "X axis range specified as 'max,min'. default: automatic",
                'ylim'        , 'y', 2, "character", "Y axis range specified as 'max,min'. default: automatic",
                'zeroy'       , 'z', 2, "logical",   "Force y-axis scale to include zero",
                'extremes'    , 'e', 2, "logical",   "Highlight local minima/maxima in plot",
                'extremevals' , 't', 2, "logical",   "Write local minima/maxima to file",
                'raw'         , 'r', 2, "logical",   "Add the raw data points to the plot",
                'datacol'     , 'd', 2, "integer",   "Number of first column with data. Default: 2",
                'hlines'      , 'k', 2, "character", "Optional comma-separated list of coordinates for horizontal lines",
                'vlines'      , 'l', 2, "character", "Optional comma-separated list of coordinates for vertical lines",
                'invert'      , 'n', 2, "character", "Optional comma separated list of 1's and 0's indicating which lines to invert",
                'ftscale'     , 's', 2, "numeric",   "Scaling factor for feature length relative to flanking regions. Default: 1",
                'help'        , 'h', 2, "logical",   "Brief help message"
               ), ncol=5,byrow=T);
opt  = getopt(args);

# Help message
if ( !is.null(opt$help) || is.null(opt$input) || is.null(opt$plotname)) {
   self=sub("^.*/", "", commandArgs()[4], perl=T);
   cat("\n", getopt(args, command=self, usage=T), "\n");
   q(status=1);
}

# Specify default argument values
if ( is.null(opt$plotformat)  ) { opt$plotformat  = "png"        }
if ( is.null(opt$fit)  )        { opt$fit         = "lines"      }
if ( is.null(opt$measure)  )    { opt$measure     = "mean"       }
if ( is.null(opt$datacol)  )    { opt$datacol     = 2            }
if ( is.null(opt$bw)       )    { opt$bw          = 1            }
if ( is.null(opt$ftscale)  )    { opt$ftscale     = 1            }


###############
# SUBROUTINES #
###############

geomean <- function (x, na.rm=T){
   return(mean(log(x+1),na.rm=na.rm));
}

meanError <- function (x, na.rm=T){
   n = sum(!is.na(x));
   m = mean(x, na.rm=na.rm);
   s = sd(x, na.rm=na.rm);
   if(n>1){
      return(qt(0.975,df=n-1)*s/sqrt(n));
   } else{
      return(NA);
   }
}

medianQuantile <- function(x, probs=0.25, na.rm=T){
   return( ifelse(any(!is.na(x)), as.numeric(quantile(x, probs=probs, na.rm=na.rm)), NA) );
}

bwAverage <- function(x, d=matrix(), measure="mean", bw=1, probs=0.25, na.rm=T){
   bin.flank = bw - 1;
   bin.s = x - bin.flank;
   bin.e = x + bin.flank;
   if (bin.s < 1)       bin.s = 1;
   if (bin.e > ncol(d)) bin.e = ncol(d);
   if ( as.character(measure) == "medianQuantile"){
      return(get(measure, mode="function")(as.numeric(d[,bin.s:bin.e]), probs=probs, na.rm=T));
   }
   else{
      return(get(measure, mode="function")(as.numeric(d[,bin.s:bin.e]), na.rm=T));
   }
}

# Read raw data matrices
readRawData <- function(filenames=c(), datacol=2, coordinates=F){
   out = list(plot.x=NULL, plot.y=NULL, xmin=NULL, xmax=NULL, ymax=NULL, ymin=NULL);
   for(i in 1:length(filenames)){
      tmp.matrix = read.delim(filenames[i], header=T, check.names=F, stringsAsFactors=F);
      tmp.matrix = tmp.matrix[,datacol:ncol(tmp.matrix)];
      
      # Check if the matrix contains feature data (with or without flanking)
      feature.match = grep("^ft\\.\\d+$",names(tmp.matrix), perl=T);
      feature.nbins = length(feature.match);
      if(feature.nbins>0){
         # No flanking data
         if (feature.nbins == ncol(tmp.matrix)){
            names(tmp.matrix) = seq((1/feature.nbins)/2,1-((1/feature.nbins)/2), 1/feature.nbins)
         } 
         else {  # Flanking data present
            flank.start = max(feature.match) + 1;
            flank.end   = max(feature.match) + feature.nbins;
            if(flank.end == ncol(tmp.matrix)){
               shift = as.numeric(names(tmp.matrix)[flank.start]) + as.numeric(names(tmp.matrix)[flank.end]) * opt$ftscale;
               names(tmp.matrix)[feature.match]         = as.numeric(names(tmp.matrix)[flank.start:flank.end]) * opt$ftscale;
               names(tmp.matrix)[flank.start:flank.end] = as.numeric(names(tmp.matrix)[flank.start:flank.end]) + shift; 
            }
            else{
               stop(paste("Invalid feature header for file", filenames[i], "\n"));
            }
         }
      }
      
      # Set up X and Y coordinate data
      out$plot.x[[i]] = as.numeric(names(tmp.matrix));
      out$plot.y[[i]] = as.matrix(tmp.matrix);
      if (any(is.na(out$plot.x[[i]]))) stop(paste("Could not parse coordinate header for", filenames[i], "\n"));
      
      # Set xmin and xmax
      tmp.max       = max(out$plot.x[[i]]);
      tmp.min       = min(out$plot.x[[i]]);
      if(is.null(out$xmax)){
         out$xmax = tmp.max;
      }
      else{
         out$xmax = ifelse(tmp.max > out$xmax, tmp.max, out$xmax);
      }
      if(is.null(out$xmin)){
         out$xmin = tmp.min;
      }
      else{
         out$xmin = ifelse(tmp.min < out$xmin, tmp.min, out$xmin);
      }
      
      # Find ymax and ymin for raw data
      tmp.max       = max(out$plot.y[[i]], na.rm=T);
      tmp.min       = min(out$plot.y[[i]], na.rm=T);
      if(is.null(out$ymax)){
         out$ymax = tmp.max;
      }
      else{
         out$ymax = ifelse(tmp.max > out$ymax, tmp.max, out$ymax);
      }
      if(is.null(out$ymin)){
         out$ymin = tmp.min;
      }
      else{
         out$ymin = ifelse(tmp.min < out$ymin, tmp.min, out$ymin);
      }
   }
   return(out);
}

# Smooth matrices
smoothData <- function(rawdata, measure="mean", interval=F, bw=1){
   out = list(plot.x=NULL, plot.y=NULL, plot.yup=NULL, plot.ydn=NULL, xmin=NULL, xmax=NULL, ymin=NULL, ymax=NULL);
   for(i in 1:length(rawdata$plot.y)){
      # Smooth data
      out$plot.x[[i]] = rawdata$plot.x[[i]];
      positions       = 1:length(out$plot.x[[i]]);
      out$plot.y[[i]] = sapply(positions, bwAverage, rawdata$plot.y[[i]], measure=measure, bw=bw, na.rm=T);
      
      # Calculate confidence intervals
      if(interval){
         if(measure == "mean"){
            plot.error = sapply(positions, bwAverage, rawdata$plot.y[[i]], measure="meanError", bw=bw, na.rm=T);
            out$plot.ydn[[i]]   = out$plot.y[[i]] - plot.error;
            out$plot.yup[[i]]   = out$plot.y[[i]] + plot.error;
         }else if (measure == "geomean"){
            plot.error = sapply(positions, bwAverage, log(rawdata$plot.y[[i]]+1), measure="meanError", bw=bw, na.rm=T);
            out$plot.ydn[[i]]   = out$plot.y[[i]] - plot.error;
            out$plot.yup[[i]]   = out$plot.y[[i]] + plot.error;
         }else if (measure == "median"){
            out$plot.ydn[[i]] = sapply(positions, bwAverage, rawdata$plot.y[[i]], measure="medianQuantile", bw=bw, probs=0.25, na.rm=T);
            out$plot.yup[[i]] = sapply(positions, bwAverage, rawdata$plot.y[[i]], measure="medianQuantile", bw=bw, probs=0.75, na.rm=T);
         }else {
            stop(paste("Error: confidence interval plotting not supported for measure '", measure, "'.\n"));
         }
      }
      
      # Exponentiate the geometric means
      if (measure == "geomean"){
         out$plot.y[[i]]   = exp(out$plot.y[[i]])-1;
         if(interval){
            out$plot.ydn[[i]] = exp(out$plot.ydn[[i]])-1;
            out$plot.yup[[i]] = exp(out$plot.yup[[i]])-1;
         }
      }
      
      # Carry over the xmax and xmin from the rawdata object
      out$xmin        = rawdata$xmin;
      out$xmax        = rawdata$xmax;
      
      # Find ymax and ymin for averaged dataset
      tmp.max       = max( c(out$plot.y[[i]], out$plot.yup[[i]], out$plot.ydn[[i]]), na.rm=T);
      tmp.min       = min( c(out$plot.y[[i]], out$plot.yup[[i]], out$plot.ydn[[i]]), na.rm=T);
      if(is.null(out$ymax)){
         out$ymax = tmp.max;
      }
      else{
         out$ymax = ifelse(tmp.max > out$ymax, tmp.max, out$ymax);
      }
      if(is.null(out$ymin)){
         out$ymin = tmp.min;
      }
      else{
         out$ymin = ifelse(tmp.min < out$ymin, tmp.min, out$ymin);
      }
   }
   return(out);
}

# Find extreme values
msExtrema <- function(x, span=3) {
  # find local maxima
  index1 <- peaks(x, span=span, strict=FALSE)

  # find local minima
  index2 <- peaks(-x, span=span, strict=FALSE)

  # remove the interior of plateaus
  index.max <- index1 & !index2
  index.min <- index2 & !index1

  # construct output
  list(index.max=index.max, index.min=index.min)
}


########
# MAIN #
########

# Check arguments
if(!( opt$measure %in% c('mean','geomean','median','sum')   )) stop("Unknown measure, must be either 'mean', 'geomean', 'median' or 'sum'");
if(!( opt$fit %in% c('lines', 'points', 'spline') )) stop("Unknown fit type, must be either 'lines', 'points' or 'spline'");
interval    = ifelse(is.null(opt$interval), FALSE, TRUE);
if(interval & (opt$fit == 'spline')) stop("Confidence intervals can only be plotted for lines or points");

# Read the raw data files
filenames   = unique(unlist(strsplit(opt$input, split=",")));
data.raw    = readRawData(filenames, datacol=opt$datacol);

# Set xlim here
if (is.null(opt$xlim)){
   xlim = c(data.raw$xmin, data.raw$xmax);
} else{
   xlim = as.numeric(unlist(strsplit(opt$xlim, split=",")));
   if (length(xlim) != 2) stop("Error: xlim must have two values 'min,max'");
   xlim = sort(xlim);
}

# Smooth the matrix files
data.smooth = smoothData(data.raw, measure=opt$measure, interval=interval, bw=opt$bw);

# Optionally invert the sign of plotted data
if(!is.null(opt$invert)){
   invert = as.numeric(unlist(strsplit(opt$invert, split=",")));
   if(!all(invert %in% c(1,0))) stop("Error: invert argument must be comma-separated list of 1s and 0s");
   if(length(filenames)!=length(invert)) stop("Error: invert list must match the number of datasets to plot");
   for (i in 1:length(invert)){
      if(invert[i]){ 
         data.raw$plot.y[[i]]    = -data.raw$plot.y[[i]];
         data.smooth$plot.y[[i]] = -data.smooth$plot.y[[i]];
         if(interval){
            tmp.ydn = data.smooth$plot.ydn[[i]];
            data.smooth$plot.ydn[[i]] = -data.smooth$plot.yup[[i]];
            data.smooth$plot.yup[[i]] = -tmp.ydn;
         }
         data.raw$ymin    = min( c(data.raw$plot.y[[i]], data.raw$ymin), na.rm=T);
         data.smooth$ymin = min( c(data.smooth$ymin, data.smooth$plot.y[[i]], data.smooth$plot.yup[[i]], data.smooth$plot.ydn[[i]]), na.rm=T);
      }
   }
}

# Set ylim
plot.raw = ifelse(is.null(opt$raw), FALSE, TRUE);
if (!is.null(opt$ylim)){
   ylim = as.numeric(unlist(strsplit(opt$ylim, split=",")));
   if (length(ylim) != 2) stop("Error: ylim must have two values 'min,max'");
   ylim = sort(ylim);
}else{
   if(plot.raw){
      ylim=c( min(c(data.raw$ymin,data.smooth$ymin), na.rm=T), max(c(data.raw$ymax,data.smooth$ymax), na.rm=T) );
   }else{
      ylim=c(data.smooth$ymin, data.smooth$ymax);
   }
   if(!is.null(opt$zeroy)){
      if (ylim[2]<0) ylim[2] = 0;
      if (ylim[1]>0) ylim[1] = 0;
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

# Plot raw data first to make sure the lines & polygons are overplotted
if(plot.raw){
   for(i in 1:length(filenames)){
      plot.x = rep(data.raw$plot.x[[i]], nrow(data.raw$plot.y[[i]]));
      points(plot.x, as.numeric(data.raw$plot.y[[i]]), pch=".", col=colors.tr[i]);
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
      if(length(ranges) %% 2 != 0) stop("Assertion thrown: ranges!=even in polygon plotting\n");
      ranges = matrix(ranges,ncol=2,byrow=T);
      
      # Plot the ranges
      for(index in 1:nrow(ranges)){
         rs = ranges[index,1];
         re = ranges[index,2];
         if ( re-rs > 0 ){
            pol.x = c(data.smooth$plot.x[[i]][rs:re], rev(data.smooth$plot.x[[i]][rs:re]));
            pol.y = c(data.smooth$plot.ydn[[i]][rs:re],rev(data.smooth$plot.yup[[i]][rs:re]));
            polygon(pol.x, pol.y, col=colors.tr[i], border=NA);
         }
      }
   }
}

# Plot the matrices
for(i in 1:length(filenames)){
   extremes.x = data.smooth$plot.x[[i]];
   extremes.y = data.smooth$plot.y[[i]];
   if(opt$fit == 'points'){
      points(data.smooth$plot.x[[i]],data.smooth$plot.y[[i]], pch=20, col=colors.op[i]);
   }
   if(opt$fit == 'lines' ){
      lines(data.smooth$plot.x[[i]],data.smooth$plot.y[[i]], pch=20, col=colors.op[i]);
   }
   if(opt$fit == 'spline'){
      spline.out  = spline(data.smooth$plot.x[[i]],data.smooth$plot.y[[i]], n=10*length(data.smooth$plot.x[[i]]));
      lines(spline.out, col=colors.op[i]);
      extremes.x = spline.out$x;
      extremes.y = spline.out$y;
   }
   
   # Find local minima and maxima
   if (!is.null(opt$extremes) | !is.null(opt$extremevals)){
      extremes.span   = ifelse('spline' %in% opt$fit, 5, 7);
      extremes        = msExtrema(extremes.y, span=extremes.span);
      extremes.min.xy = cbind(extremes.x[extremes$index.min], extremes.y[extremes$index.min]);
      extremes.max.xy = cbind(extremes.x[extremes$index.max], extremes.y[extremes$index.max]);
      
      if (!is.null(opt$extremes)){
         points(extremes.min.xy, pch=4, col=colors.op[i]);
         points(extremes.max.xy, pch=3, col=colors.op[i]);
      }
      if (!is.null(opt$extremevals)){
         write.table(extremes.min.xy, file=paste(filenames[i], ".extremes.min", sep=""), row.names=F, col.names=F, quote=F, sep="\t");
         write.table(extremes.max.xy, file=paste(filenames[i], ".extremes.max", sep=""), row.names=F, col.names=F, quote=F, sep="\t");
      }
   }
}

# Print the legend
if (opt$plotformat == "png"){
   par(mar = c(5, 0, 4, 2) + 0.1, cex=0.8);
}else{
   par(mar = c(5, 0, 4, 2) + 0.1, cex=0.5);
}
plot(1,1, type="n", xlim=c(0,10), ylim=c(0,10), axes=F, ann = FALSE);
legend("topleft", legend=c(filenames), pch=20, col=colors.op);
text(0,0.3, paste("measure: ", opt$measure), pos=4);

# Close the plot
dev = dev.off();
