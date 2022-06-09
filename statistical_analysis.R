###### Statistical Analysis Part ######
#!/usr/bin/env Rscript
### load input variables ###
args = commandArgs(trailingOnly=TRUE)

for(v in args)
{
  vTmp <- strsplit(v,"=")[[1]]
  assign(vTmp[[1]],vTmp[[2]])
}

# Libraries
library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(parallel)
library(PRROC)
library(stringr)

## 1) Binning the genome in windows of length w
w <- binLength # binLength parameter from outside
genesBed <- read.table(genesbed, sep = "\t") # bed_genome parameter from outside
colnames(genesBed) <- c("chr","start","end","name","score","strand")

genesBinsList <- mclapply(1:nrow(genesBed),function(k) 
{
  i <- genesBed[k,]
  
  lTmp <- (i$end - i$start)+1
  nTmp <- floor(lTmp/w)
  
  swTmp <- lTmp - (nTmp*w)
  
  if(lTmp<=w)
  {
    grangeTmp <- GRanges(seqnames=i$chr,ranges=IRanges(i$start,i$end),strand=i$strand)
    names(grangeTmp) <- paste0(i$name,"-B1")
    grangeTmp
  }else{
    if(as.character(i$strand)=="+")
    {
      if(swTmp!=0)
      {
        breaksTmp <- c(1,seq(from=swTmp,to=lTmp,by=w))+(i$start-1)
      }else{
        breaksTmp <- seq(from=swTmp,to=lTmp,by=w)+(i$start-1)
      }
      grangeTmp <- GRanges(seqnames=i$chr,ranges=IRanges(breaksTmp[-length(breaksTmp)],breaksTmp[-1]-1),strand=i$strand)
      names(grangeTmp) <- paste0(i$name,"-B",seq_along(breaksTmp[-1]-1))
      grangeTmp
    }else{
      if(swTmp!=0)
      {
        breaksTmp <- c(seq(from=1,to=(lTmp-swTmp+1),by=w),lTmp)+(i$start-1)
      }else{
        breaksTmp <- seq(from=1,to=(lTmp-swTmp+1),by=w)+(i$start-1)
      }
      
      grangeTmp <- GRanges(seqnames=i$chr,ranges=IRanges(breaksTmp[-length(breaksTmp)],breaksTmp[-1]-1),strand=i$strand)
      names(grangeTmp) <- paste0(i$name,"-B",seq_along(breaksTmp[-1]-1))
      grangeTmp
    }
  }
},mc.cores=6)

genesBins <- unlist(as(genesBinsList,"GRangesList"))
genesBins <- genesBins[which(width(genesBins) == w), ] # Remove bins without length equal to w

## 2) Load all the output files from the output directory
files <- list.files(bed_folder, full.names = TRUE) # output_directory parameter from outside

listmax <- paste0(c("dena", "epinanoErr", "epinanoSvm", "nanodoc"), collapse = "|") # for these tools we need to maximize the filtering paramenter when there are more than 1 in a bin
listmin <- paste0(c("differ", "drummer", "yanocomp", "nanocompore", "eligos", "xpore", "tomboComparison"), collapse = "|") # for these tools we need to minimize the filtering parameter

## 3) Convert .bed output file of each tool + peaks file into a Granges object
# Peaks
peaks_bed <- read.table(peaks, header = FALSE, sep = "\t") # MiCLIP and MAZTER-seq peaks file parameter from outside
colnames(peaks_bed) <- c("chr", "start", "end", "desc", "score", "strand")
peaks_granges <- makeGRangesFromDataFrame(peaks_bed)

# Build matrix of zeros
hitsMatrix <- matrix(c(0), nrow = length(genesBins), ncol = length(files) + 1)
colnames(hitsMatrix) <- c("Gold_standard", basename(files))
row.names(hitsMatrix) <- c(1:nrow(hitsMatrix))

# Overlap between peaks of gold_standard and genome binned
peaks_overlap <- as.matrix(findOverlaps(query = peaks_granges, subject = genesBins, minoverlap = 1, type = "any"))
hitsMatrix[peaks_overlap[, 2], "Gold_standard"] <- 1

matrix_nanom6A <- hitsMatrix[,"Gold_standard"]
names_nanom6A <- c()

listF1score <- list()
listPRcurves <- list()
noise <- 0.000001
threshold_default <- c(0.1, 0.05, 0.01, 0.05, 0.01, 0.001, 0.1, 0.5, 0.05, 0.02, 0.05)
names(threshold_default) <- c("dena_output.bed", "differr_output.bed", "drummer_output.bed", "yanocomp_output.bed", "nanocompore_output.bed", "eligos_output.bed", 
                                  "epinanoErr_output.bed", "epinanoSvm_output.bed", "xpore_output.bed", "nanodoc_output.bed", "tomboComparison_output.bed")

for (y in files) {
  x <- basename(y)
  # Extraction of bed file + convertion to granges
  bed_file <- read.table(y, header = T, sep = "\t")
  granges <- makeGRangesFromDataFrame(bed_file, keep.extra.columns = T)
  # Overlap between m6A detected site of each tool and genome binned
  overlap <- as.matrix(findOverlaps(query = granges, subject = genesBins, minoverlap = 1, type = "any"))
  if (grepl(x, pattern = paste0(c("dena","drummer","differr","yanocomp","nanocompore","eligos","epinanoErr",
                                  "epinanoSvm","xpore","nanodoc","tomboComparison"), collapse = "|"))) {
    # Add column of filtering parameter
    filtering_parameter <- bed_file[overlap[,"queryHits"] , 6]
    overlap_w_parameter <- cbind(overlap, filtering_parameter)
    default <- unname(threshold_default[grep(x, pattern = paste0(c("dena","drummer","differr","yanocomp","nanocompore","eligos","epinanoErr","epinanoSvm",
                                                                                                                                "xpore","nanodoc","tomboComparison"), collapse = "|"), value = T)])
    # Recognize from the name of the tool if we need to keep the maximum or minimum value (when there are more hits in a single bin)
    if(grepl(x, pattern = listmax)) { 
      score <- sapply(split(overlap_w_parameter[,3], overlap_w_parameter[,2]), max) + noise
      default_thr <- default
    } else {
      score <- -1*sapply(split(overlap_w_parameter[,3], overlap_w_parameter[,2]), min) - noise
      default_thr <- - default
    }
    hitsMatrix[, x] <- rep((min(score) - noise), length(genesBins))
    hitsMatrix[names(score), x] <- score
    positive <- hitsMatrix[which(hitsMatrix[,"Gold_standard"] == 1), x]
    negative <- hitsMatrix[which(hitsMatrix[,"Gold_standard"] == 0), x]
    #par(mfrow = c(2, 1))
    #pdf(file = paste0(x,"_scores_distribution.pdf"))
    #hist(positive, main = paste0(x, " - Scores for positive peaks"))
    #hist(negative, main = paste0(x, " - Scores for negative peaks"))
    #dev.off()
    pr <- pr.curve(scores.class0 = unname(positive), scores.class1 = unname(negative), curve=T)
    pdf(file = paste0(x,"_PRcurve.pdf"), width = 8, height = 8)
    plot(pr, main = paste0(x, " Precision-Recall curve"))
    dev.off()
    listPRcurves[[x]] <- pr
    # Plot "manual" PR curve
    thresholds <- c(seq(from = min(score), to = max(score), length.out = 10000), default_thr)
    for (t in 1:length(thresholds)) {
      thr <- thresholds[t]
      TP <- length(which(positive >= thr))
      FN <- length(which(positive < thr))
      FP <- length(which(negative >= thr))
      TN <- length(which(negative < thr))
      recall <- c(recall, TP/(TP + FN))
      precision <- c(precision, TP/(TP + FP))
    }
    names(recall) <- thresholds
    names(precision) <- thresholds
    # F1 score
    F1score <- 2*(precision[as.character(default_thr)]*recall[as.character(default_thr)])/(precision[as.character(default_thr)]+recall[as.character(default_thr)])
    listF1score[[x]] <- F1score 
    #pdf(paste0(x, "_PRcurve_manual.pdf"))
    #plot(recall, precision, main = paste0(x, " - PR manual"), type = "l", xlim = c(0, 1), ylim = c(0, 1))
    #dev.off()
    }
  } 
  else if (grepl(x, pattern = "mines")) {
    hitsMatrix[overlap[, 2], x] <- 1
    TP <- 0
    for (y in 1:nrow(hitsMatrix)){
      if ((hitsMatrix[y, "Gold_standard"] == 1) && (hitsMatrix[y, x] == 1)){
        TP <- TP + 1
      }
    }
    # Recall + Precision
    totPositiveGS <- count(hitsMatrix[, "Gold_standard"] == 1)
    totPositiveTool <- count(hitsMatrix[, x] == 1)
    recall <- TP/totPositiveGS
    precision <- TP/totPositiveTool
    # F1 score
    F1score <- 2*(precision*recall)/(precision+recall)
    names(F1score) <- "default"
    listF1score[[x]] <- F1score
  } 
  else if (grepl(x, pattern = "nanom6a")) {
    hitsMatrix[overlap[, 2], x] <- 1
    matrix_nanom6A <- cbind(matrix_nanom6A, hitsMatrix[, x])
    names_nanom6A <- c(names_nanom6A, x)
    if (grepl(x, pattern = "0.5")){
      TP <- 0
      for (y in 1:nrow(hitsMatrix)){
        if ((hitsMatrix[y, "Gold_standard"] == 1) && (hitsMatrix[y, x] == 1)){
          TP <- TP + 1
        }
      }
      # Recall + Precision
      totPositiveGS <- count(hitsMatrix[, "Gold_standard"] == 1)
      totPositiveTool <- count(hitsMatrix[, x] == 1)
      recall <- TP/totPositiveGS
      precision <- TP/totPositiveTool
      # F1 score
      F1score <- 2*(precision*recall)/(precision+recall)
      names(F1score) <- "default"
      listF1score[[x]] <- F1score
    }
  }
}

### Code for PR curve nanom6A which has multiple files each run with a different threshold
names_nanom6A <- str_extract(names_nanom6A, "[0-9]\\.[0-9]*")
colnames(matrix_nanom6A) <- c("GS", names_nanom6A)

max_thr <- c()
short <- matrix_nanom6A[,2:ncol(matrix_nanom6A)]

for (row in 1:nrow(short)) {
 v <- c()
 if (sum(short[row,]) != 0){
   for (col in 1:ncol(short)) {
     if(short[row,col] == 1){
       v <- c(v, as.numeric(colnames(short)[col]))
     }
   }
   max <- max(v)
   max_thr <- c(max_thr, max)
 }
 else{
   max_thr <- c(max_thr, 0)
 }
}

new_matrix <- cbind(matrix_nanom6A, max_thr)
posit <- new_matrix[which(new_matrix[, "GS"] == 1), "max_thr"]
negat <- new_matrix[which(new_matrix[, "GS"] == 0), "max_thr"]

pr <- pr.curve(posit, negat, curve = T)
pdf(file = paste0(resultsFolder, "nanom6A_PRcurve.pdf"), width = 8, height = 8)
plot(pr, main = "nanom6A Precision-Recall curve")
dev.off()
listPRcurves[["nanom6A_output.bed"]] <- pr
#par(mfrow = c(2, 1))
#pdf(file = paste0(resultsFolder, "nanom6A_scores_distribution.pdf"))
#hist(posit, 100,main = "nanom6A - Scores for positive peaks")
#hist(negat, 100,main = "nanom6A - Scores for negative peaks")
#dev.off()

### Save list of F1 scores
capture.output(listF1score, file = "F1score_list.csv")

### Plot all the Precision-Recall curves together
col <- c(7,8,420,153,31,100,33,47,53,62,400,454,28,10)
pdf(file = "Summary PR curves.pdf", width = 8, height = 8)
for (x in 1:length(listPRcurves)) {
  if (x == 1){
    plot(listPRcurves[[x]], color = colors()[col[x]], main = "Summary PR curves")
  }
  else{
    plot(listPRcurves[[x]], add = T, color = colors()[col[x]], main = "Summary PR curves")
  }
}
legend("bottomright", legend = names(listPRcurves), col = colors()[col[1:length(listPRcurves)]], lty=1:1, cex=0.8, bg = "lightblue" )
dev.off()
