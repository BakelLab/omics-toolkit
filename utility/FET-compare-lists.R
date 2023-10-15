#!/usr/bin/env Rscript

# 22.07.2016 
# Zuleyma Peralta <zyperalta@gmail.com>

#############
# ARGUMENTS #
#############

#  Col 1: the long flag name. A multi-character string.
#  Col 2: short flag alias of Column 1. A single-character string.
#  Col 3: 0=no argument, 1=required argument, 2=optional argument.
#  Col 4: Possible values: logical, integer, double, complex, character.
#  Col 5: A brief description of the purpose of the option.
suppressMessages(library(getopt));
args = matrix(c('listA'      , 'a', 1, "character", "Gene set A",
                'listB'      , 'b', 2, "character", "Gene set B",
                'background' , 'x', 2, "character", "Background genes",
                'maxlogp'    , 'm', 2, "numeric",   "Set upper threshold to -logP score (applied to plot only). Default: none",
                'output'     , 'o', 2, "character", "Output file",
                'oddsratio'  , 'r', 0, "logical",   "Use odds ratio as plot label",
                'help'       , 'h', 0, "logical",   "Brief help message"
               ), ncol=5,byrow=T);
opt  = getopt(args);

# Set defaults
if ( is.null(opt$output) )  { opt$output        = "set-overlap.txt"  }

# Help message
if ( is.null(opt$listA) || is.null(opt$listB)  || is.null(opt$background) ) {
   self=sub("^.*/", "", commandArgs()[4], perl=T);
   cat("\n", getopt(args, command=self, usage=T), "\n");
   q(status=1);
}

#################
### LIBRARIES ###
#################

library(ggplot2)
library(parallel)
library(foreach)
library(iterators)


###############
## FUNCTIONS ##
###############

# These functions were provided by the authors of MEGENA.

# Read Genesets
read.geneSet <- function(geneSet.file) {
   gene.list <- readLines(geneSet.file)
   gene.list <- strsplit(gene.list,"\t")
   names(gene.list) <- sapply(gene.list,function(x) x[1])
   gene.list <- lapply(gene.list,function(x) x[2:length(x)])
   return(gene.list)
}

# Output Genesets
output.geneSet.file <- function(geneSet,outputfname) {
   if (!is.list(geneSet)) stop("geneSet is not a list.")
   if (is.null(names(geneSet))) stop("names(geneSet) is not defined properly.")

   sink(outputfname)
   cat(paste(paste(names(geneSet),"\t",sapply(geneSet,function(x) paste(x,collapse = "\t")),sep = ""),collapse = "\n"))
   sink()

   return(0)
}

# selected pairs by ij
make.paired.Tables <- function(geneSets1,geneSets2,ij,background) {
   #mem1 <- do.call(cbind,lapply(geneSets1,function(x,y) {vec <- rep(0,length(y));vec[which(y %in% x)] <- 1;return(vec)},y = background))
   #mem2 <- do.call(cbind,lapply(geneSets2,function(x,y) {vec <- rep(0,length(y));vec[which(y %in% x)] <- 1;return(vec)},y = background))

   pairwise.tables <- mapply(FUN = function(s1,s2,z) {
                                         v1 <- rep(0,length(z));v1[which(z %in% s1)] <- 1; 
										 v2 <- rep(0,length(z));v2[which(z %in% s2)] <- 1;
										 # n11,n12,n21,n22
										 as.table(matrix(c(sum(abs(v1 - 1) * abs(v2 - 1)),sum(abs(v2 - 1) * v1),sum(abs(v1 - 1) * v2),sum(v1 * v2)),nrow = 2))
                                },s1 = geneSets1[ij[,1]],s2 = geneSets2[ij[,2]],MoreArgs = list(z = background),SIMPLIFY = FALSE)
   names(pairwise.tables) <- paste(names(geneSets1)[ij[,1]],"~vs~",names(geneSets2)[ij[,2]],sep = "")
   return(pairwise.tables)
}

# Make pairwise tables
make.Pairwise.Tables <- function(geneSets1,geneSets2,background) {
   mem1 <- do.call(cbind,lapply(geneSets1,function(x,y) {vec <- rep(0,length(y));vec[which(y %in% x)] <- 1;return(vec)},y = background))
   mem2 <- do.call(cbind,lapply(geneSets2,function(x,y) {vec <- rep(0,length(y));vec[which(y %in% x)] <- 1;return(vec)},y = background))

   d <- t(mem1) %*% mem2;
   b <- abs(t(mem1) %*% (mem2-1))
   c <- abs(t(mem1-1) %*% (mem2))
   a <- t(mem1-1) %*% (mem2-1);

   ij <- do.call(rbind,lapply(1:length(geneSets1),function(i,j) cbind(rep(i,length(j)),j),j = 1:length(geneSets2)))

   pairwise.tables <- lapply(1:nrow(ij),function(i,ij,a,b,c,d) as.table(matrix(c(a[ij[i,1],ij[i,2]],b[ij[i,1],ij[i,2]],c[ij[i,1],ij[i,2]],d[ij[i,1],ij[i,2]]),nrow = 2)),ij = ij,a = a,b = b,c = c,d = d)
   names(pairwise.tables) <- paste(names(geneSets1)[ij[,1]],"~vs~",names(geneSets2)[ij[,2]],sep = "")
   return(pairwise.tables)
}

# Do FisherExactTest
do.FisherExactTest <- function(table.count,N_bg = NULL) {
   if (is.null(N_bg)) N_bg = sum(rowSums(table.count))

   out <- fisher.test(x = table.count,or = 1,alternative = "greater")
   odds.ratio <- out$estimate
   p.value <- out$p.value;
   geneSet1.count <- rowSums(table.count)[2]
   geneSet2.count <- colSums(table.count)[2]
   expected.count <- geneSet1.count/N_bg * geneSet2.count
   overlap.count <- table.count[2,2];
   fold.change <- overlap.count/expected.count

   out <- c(N_bg,geneSet1.count,geneSet2.count,expected.count,overlap.count,fold.change,odds.ratio,p.value)
   names(out) <- c("Background","set1_size","set2_size","expected.overlap","actual.overlap","enrichment.foldchange","odds.ratio","FET_pvalue")
   return(out)
}

# Perform an AllPairs FET
perform.AllPairs.FET <- function(geneSets1,geneSets2,background,adjust.FET.pvalue = T,do.multicore = F,n.cores = NULL) {
   pairwise.tables <- make.Pairwise.Tables(geneSets1,geneSets2,background)

   if (do.multicore) {
      #registerDoMC(n.cores)
      set.parallel.backend(n.cores)
      fact <- do.call(c,lapply(1:n.cores,function(i,dn) rep(i,dn),dn = ceiling(length(pairwise.tables)/n.cores)))
      fact <- factor(fact[1:length(pairwise.tables)])
      split.tables <- lapply(split(1:length(pairwise.tables),fact),function(i,obj) obj[i],obj = pairwise.tables);rm(fact)
      output <- foreach (tbl = split.tables,.combine = 'c') %dopar% {
         out <- lapply(tbl,do.FisherExactTest) 
			return(out)
      }
   } else {
      output <- lapply(pairwise.tables,do.FisherExactTest)
   }

   # collect outputs
   output <- data.frame(set1_Name = gsub("~vs~(.*)$","",names(output)),set2_Name = gsub("^(.*)~vs~","",names(output)),as.data.frame(do.call(rbind,output)))
   int.element <- mapply(FUN = function(x,y) paste(intersect(x,y),collapse = ","),x = geneSets1[as.character(output[[1]])],y = geneSets2[as.character(output[[2]])],SIMPLIFY = F)
   int.element <- do.call(c,int.element)
   if (adjust.FET.pvalue) output <- cbind.data.frame(output,data.frame(corrected.FET.pvalue = p.adjust(output$FET_pvalue,"bonferroni"),intersects = int.element))
   return(output)
}

# Perform an ijPairs FET
perform.ijPairs.FET <- function(geneSets1,geneSets2,ij,background,adjust.FET.pvalue = T,do.multicore = F,n.cores = NULL) {
   pairwise.tables <- make.paired.Tables(geneSets1,geneSets2,ij,background)

   if (do.multicore) {
      #registerDoMC(n.cores)
      set.parallel.backend(n.cores)
      fact <- do.call(c,lapply(1:n.cores,function(i,dn) rep(i,dn),dn = ceiling(length(pairwise.tables)/n.cores)))
      fact <- factor(fact[1:length(pairwise.tables)])
      split.tables <- lapply(split(1:length(pairwise.tables),fact),function(i,obj) obj[i],obj = pairwise.tables);rm(fact)
      output <- foreach (tbl = split.tables,.combine = 'c') %dopar% {
         out <- lapply(tbl,do.FisherExactTest) 
			return(out)
      }
   } else {
      output <- lapply(pairwise.tables,do.FisherExactTest)
   }

   # collect outputs
   output <- data.frame(set1_Name = gsub("~vs~(.*)$","",names(output)),set2_Name = gsub("^(.*)~vs~","",names(output)),as.data.frame(do.call(rbind,output)))
   int.element <- mapply(FUN = function(x,y) paste(intersect(x,y),collapse = ","),x = geneSets1[as.character(output[[1]])],y = geneSets2[as.character(output[[2]])],SIMPLIFY = F)
   int.element <- do.call(c,int.element)
   if (adjust.FET.pvalue) output <- cbind.data.frame(output,data.frame(corrected.FET.pvalue = p.adjust(output$FET_pvalue,"bonferroni"),intersects = int.element))
   return(output)
}

# Read MSigDB
read.MSigDB <- function(MSigDB.file) {
   gene.list <- readLines(MSigDB.file)
   gene.list <- strsplit(gene.list,"\t")

   names(gene.list) <- sapply(gene.list,function(x) x[1])
   gene.list <- lapply(gene.list,function(x) x[3:length(x)])
   return(gene.list)
}

# Perform FET output
perform.FET.output <- function(geneSets1,geneSets2,background,outputfname,adjust.FET.pvalue = T,pvalue.cutoff = 0.05,do.multicore = F,n.cores = NULL) {
   if (is.null(names(geneSets1))) {
      cat("#### No names were tagged to geneSets1 #####\nProviding names to geneSets1...\n")
      names(geneSets1) <- paste("A",1:length(geneSets1),sep = "")
   }
 
   if (is.null(names(geneSets2))) {
      cat("#### No names were tagged to geneSets2 #####\nProviding names to geneSets1...\n")
      names(geneSets2) <- paste("B",1:length(geneSets2),sep = "")
   }
   cat("###### Performing Fisher Exact Test for over-representation\n")
   cat(paste("- Number of geneSets1:",length(geneSets1),"\n","- Number of geneSets2:",length(geneSets2),"\n",sep = ""))

   FET.table <- perform.AllPairs.FET(geneSets1 = geneSets1,geneSets2 = geneSets2,background = background,adjust.FET.pvalue = adjust.FET.pvalue,do.multicore = do.multicore,n.cores = n.cores)
 
   # output gene sets
   cat("###### Outputting files...\n- Output gene sets\n")
   geneSet.files <- paste(outputfname,"_geneSets",c(1,2),".txt",sep = "")
   cat(paste(geneSet.files,"\n",sep = ""))
   output.status <- output.geneSet.file(geneSet = geneSets1,outputfname = geneSet.files[1])
   output.status <- output.geneSet.file(geneSet = geneSets2,outputfname = geneSet.files[2])
   cat("- Output FET output table\n")
   write.table(FET.table,file = paste(outputfname,"_FET-Table.txt",sep = ""),sep = "\t",row.names = F,col.names = T,quote = F)
   cat("- Output significant results\n")
   sig.table <- FET.table[FET.table$corrected.FET.pvalue < pvalue.cutoff,];
   write.table(sig.table,file = paste(outputfname,"_SignificantFET-Table.txt",sep = ""),sep = "\t",row.names = F,col.names = T,quote = F)

   return(0) 
}

# Run MSigDB, replaces run.MSigDB.Enrichment() function 
run.MSigDB <- function(module.geneSets,GeneSymbol,MSigDB.files,saveto,
                       annot.table = NULL,id.col = NULL,symbol.col = NULL,
                       do.multicore = F,n.cores = NULL) {
   dir.create(saveto)

   if (is.null(MSigDB.files)) {
      data(MSigDB.geneSets)
   } else {
      names(MSigDB.files) <- gsub("\\.gmt$","",gsub("(.*)/","",MSigDB.files))
      MSigDB.geneSets <- lapply(MSigDB.files,read.MSigDB)
   }
 
   if (is.null(GeneSymbol)) GeneSymbol <- Reduce("union",module.geneSets);

   # if annotation is provided, project MSigDB gene symbols to array
   if (!is.null(annot.table)) {
      gene.vec <- as.character(annot.table[[symbol.col]]);
      names(gene.vec) <- as.character(annot.table[[id.col]])
      for (i in 1:length(MSigDB.geneSets)) {
         MSigDB.geneSets[[i]] <- lapply(MSigDB.geneSets[[i]],function(x,y) {out <- names(y)[match(x,y)];out <- out[!is.na(out)];return(out)},y = gene.vec)
         MSigDB.geneSets[[i]] <- MSigDB.geneSets[[i]][sapply(MSigDB.geneSets[[i]],length) > 0]
      }
   }
 
   # load up MSigDB gene Sets
   for (i in 1:length(MSigDB.geneSets)) {
      output.status <- perform.FET.output(geneSets1 = module.geneSets,geneSets2 = MSigDB.geneSets[[i]],background = GeneSymbol,outputfname = paste(saveto,names(MSigDB.geneSets)[i],sep = "/"),do.multicore = do.multicore,n.cores = n.cores)
   }

   return(0)
}


##### FET table formatting functions

# Extract SigTerms
extract.SigTerms <- function(FET.file,pvalue.cutoff = 0.05,fold.cutoff = 2,pval.colname,fold.colname = "enrichment.foldchange") {
   FET.table <- read.delim(file = FET.file,sep = "\t",header = T)
   pvalue.vec <- FET.table[[which(names(FET.table) == pval.colname)]]
   fold.vec <- FET.table[[which(names(FET.table) == fold.colname)]]

   ii <- which(pvalue.vec < pvalue.cutoff & fold.vec > fold.cutoff)

   sig.module <- factor(FET.table[[1]][ii])
   sig.term <- as.character(FET.table[[2]][ii])
   sig.pval <- pvalue.vec[ii]
   sig.fold <- fold.vec[ii]

   ## order by pvalues
   ii <- order(sig.pval);
   sig.module <- sig.module[ii];sig.term <- sig.term[ii];sig.pval <- sig.pval[ii];sig.fold <- sig.fold[ii]
   sig.str <- paste(sig.term,signif(sig.pval,4),signif(sig.fold,4),sep = "/");
   split.str <- split(sig.str,sig.module)

   output <- data.frame(Module = names(split.str),Significant.Term = sapply(split.str,function(x) paste(x,collapse = ",")))
   return(output)
}

# Format FET output
format.FET.output <- function(FET.folder,pvalue.cutoff = 0.05,fold.cutoff = 2) {
   cat(paste("##### formatting FET outputs from folder:",FET.folder,"\n",sep = ""))
   FET.table.files <- list.files(path = FET.folder,pattern = "_FET-Table.txt$",full.names = T);
   cat("Processing the outputs from the following files:\n");
   cat(paste(paste(FET.table.files,collapse = "\n"),"\n",sep = ""));

   names(FET.table.files) <- gsub("_FET-Table.txt","",gsub("(.*)/","",FET.table.files))

   SigTerms <- lapply(FET.table.files,extract.SigTerms,pval.colname = "FET_pvalue",pvalue.cutoff = pvalue.cutoff,fold.cutoff = fold.cutoff)
   AllPair.SigTerms <- lapply(FET.table.files,extract.SigTerms,pval.colname = "corrected.FET.pvalue",pvalue.cutoff = pvalue.cutoff,fold.cutoff = fold.cutoff)

   module.ticks <- sort(union(Reduce("union",lapply(SigTerms,function(x) as.character(x[[1]]))),Reduce("union",lapply(AllPair.SigTerms,function(x) as.character(x[[1]])))))

   output.matrix <- matrix(NA,nrow = length(module.ticks),ncol = length(FET.table.files)*2)

   sig.i <- seq(1,ncol(output.matrix)-1,2)
   allpair.i <- seq(2,ncol(output.matrix),2)

   for (i in 1:length(FET.table.files)) {
      ii <- match(as.character(SigTerms[[i]][[1]]),module.ticks)
      if (length(ii) > 0) output.matrix[ii,sig.i[i]] <- as.character(SigTerms[[i]][[2]])
      rm(ii)

      ii <- match(as.character(AllPair.SigTerms[[i]][[1]]),module.ticks);
      if (length(ii) > 0) output.matrix[ii,allpair.i[i]] <- as.character(AllPair.SigTerms[[i]][[2]])
      rm(ii)
   }
 
   colname.vec <- rep("",ncol(output.matrix))
   colname.vec[sig.i] <- paste("SigTerms.",names(FET.table.files),sep = "")
   colname.vec[allpair.i] <- paste("AllPair.SigTerms.",names(FET.table.files),sep = "")

   colnames(output.matrix) <- colname.vec

   output.matrix <- data.frame(Module = module.ticks,as.data.frame(output.matrix))

   return(output.matrix)
}


##############
#### MAIN ####
##############

# Read genesets
set.A = read.geneSet(opt$listA);
set.B = read.geneSet(opt$listB);

# Set up background genes
bg.genes = readLines(opt$background);

# Do geneset testing and write results to table
cat( paste("Calculating enrichments\n", sep="") );
output.table = perform.AllPairs.FET(geneSets1 = set.A, geneSets2 = set.B, background = bg.genes, adjust.FET.pvalue = T,do.multicore = F,n.cores = NULL);
output.table$corrected.FET.pvalue = p.adjust(output.table$FET_pvalue,"BH");
write.table(output.table,file = opt$output ,sep = "\t",row.names = F,col.names = T,quote = F);

output.table$log10fdr = -log10(output.table$corrected.FET.pvalue);

pd      = output.table[,c("set1_Name","set2_Name","log10fdr","odds.ratio")];
if (is.null(opt$oddsratio)){
   pd$text = signif(pd$log10fdr,2);
   pd$text[ pd$log10fdr < -log10(0.05) ] = "";
   
   # Restrict max score for plotting if requested
   if ( !is.null(opt$maxlogp) ){
     pd$log10fdr[pd$log10fdr > opt$maxlogp] = opt$maxlogp
   }
} else {
   pd$text = signif(pd$odds.ratio,2);
   pd$text[ pd$log10fdr < -log10(0.05) ] = "";
}

pdf( paste(opt$output, ".pdf", sep="") )
g2 = ggplot(pd, aes(x=set2_Name,y=set1_Name, label=text)) +
        geom_tile(aes(fill=log10fdr),color="grey60") + scale_fill_gradient(low = "white", high = "red","-log10FDR\n") + 
        geom_text(size=3, color="black") + ylab("") + xlab("") + theme(axis.text.x = element_text(angle=30, hjust=1)) + 
        ggtitle("Enrichments")
print(g2)
dev.off()

