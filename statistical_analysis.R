#!/usr/bin/env Rscript
### load input variables ###
args = commandArgs(trailingOnly=TRUE)

for(v in args)
{
  vTmp <- strsplit(v,"=")[[1]]
  assign(vTmp[[1]],vTmp[[2]])
}

#avoid scientific notation
options(scipen = 100)

# Libraries
library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(parallel)
library(PRROC)
library(stringr)
library(Biostrings)
library(pheatmap)

## 1) Load all the output files from the output directory
files <- list.files(bed_folder, full.names = TRUE, pattern = "\\.bed") # output_directory parameter from outside

listmax <- paste0(c("DENA", "EpiNano-Error", "EpiNano-SVM", "NanoDoc", "m6Anet"), collapse = "|") # for these tools we need to maximize the filtering paramenter when there are more than 1 in a bin
listmin <- paste0(c("DiffErr", "DRUMMER", "Yanocomp", "Nanocompore", "ELIGOS", "xPore", "Tombo"), collapse = "|") # for these tools we need to minimize the filtering parameter

threshold_default <- c(0.1, 0.05, 0.05, 0.05, 0.01, 1, 0.1, 0.5, 0.05, 0.02, 0.05, 0.9)
names(threshold_default) <- c("DENA", "DiffErr", "DRUMMER", "Yanocomp", "Nanocompore", "ELIGOS", 
                              "EpiNano-Error", "EpiNano-SVM", "xPore", "NanoDoc", "Tombo", "m6Anet")

chrs <- readDNAStringSet(genomefile, format="fasta")
RRACH_plus <- GRanges(vmatchPattern(pattern = "RRACH", subject = chrs, fixed = "subject"), strand = "+")
RRACH_minus <- GRanges(vmatchPattern(pattern = "DGTYY", subject = chrs, fixed = "subject"), strand = "-")
RRACH <- c(RRACH_plus, RRACH_minus)
RRACH_bed <- cbind(as.data.frame(seqnames(RRACH)), start(RRACH), end(RRACH), as.data.frame(strand(RRACH)))
colnames(RRACH_bed) <- c("chr", "start", "end", "strand")
#write.table(RRACH_bed, file = "RRACH_coords.bed", sep = "\t", quote = F, row.names = F)

## 2) Binning the genome in windows of length w
w <- as.numeric(binLength) # binLength parameter from outside
genesBed <- read.table(genesbed, sep = "\t") # genesbed parameter from outside
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
}, mc.cores=as.numeric(mccores))

genesBins <- unlist(as(genesBinsList,"GRangesList"))
genesBins <- genesBins[which(width(genesBins) == w), ] # Remove bins without length equal to w

###### Find bins with RRACH motifs
RRACH_granges <- unique(sort(makeGRangesFromDataFrame(RRACH_bed)))
# Peaks
peaks_bed <- read.table(peaks, header = TRUE, sep = "\t") # gold-standard peaks file parameter from outside
colnames(peaks_bed) <- c("chr", "start", "end", "desc", "score", "strand")
peaks_granges <- unique(sort(makeGRangesFromDataFrame(peaks_bed)))
###### Find bins with RRACH motifs
RRACH_overlap <- findOverlaps(query = RRACH_granges, subject = genesBins, minoverlap = min(5, w), type = "any")
genesBins_RRACH <- unique(sort(genesBins[subjectHits(RRACH_overlap)]))
# Peaks
RRACH_overlap_genesBins_peaks <-  findOverlaps(query = genesBins_RRACH, subject = peaks_granges, minoverlap = 1, type = "any")
peaks_RRACH_granges <- genesBins_RRACH[sort(unique(queryHits(RRACH_overlap_genesBins_peaks)))]

##### Find bins with high coverage

if (!exists("highcov_bed_file")) {
  highcov_GRanges <- NULL
} else {
  highcov_bed <- read.table(highcov_bed_file, sep = "\t", header = FALSE)
  highcov_GRanges <- GRanges(seqnames = highcov_bed[, 1],
                             ranges = IRanges(start = highcov_bed[, 2], end = highcov_bed[, 3]),
                             strand = highcov_bed[, 6])
  names(highcov_GRanges) <- highcov_bed[, 4]
  highcov_overlap <- findOverlaps(query = highcov_GRanges, subject = genesBins, minoverlap = 1, type = "any")
  genesBins_highcov <- unique(sort(genesBins[subjectHits(highcov_overlap)]))
  overlap_genesBins_highcov_peaks <-  findOverlaps(query = genesBins_highcov, subject = peaks_granges, minoverlap = 1, type = "any")
  peaks_highcov_granges <- genesBins_highcov[sort(unique(queryHits(overlap_genesBins_highcov_peaks)))]
  
  ###### Find high coverage bins with RRACH motifs
  highcov_RRACH_overlap <- findOverlaps(query = RRACH_granges, subject = genesBins_highcov, minoverlap = min(5, w), type = "any")
  genesBins_highcov_RRACH <- unique(sort(genesBins_highcov[subjectHits(highcov_RRACH_overlap)]))
  genesBins_highcov_RRACH_peaks_overlap <-  findOverlaps(query = genesBins_highcov_RRACH, subject = peaks_granges, minoverlap = 1, type = "any")
  peaks_highcov_RRACH_granges <- genesBins_highcov_RRACH[sort(unique(queryHits(genesBins_highcov_RRACH_peaks_overlap)))]
}

Run_statistical_analysis <- function(genesBins_par, peaks_par, files_par, notes = "", w) {
  tools <- basename(files_par)
  tools[grep(pattern = "dena", x = tools)] <- "DENA"
  tools[grep(pattern = "drummer", x = tools)] <- "DRUMMER"
  tools[grep(pattern = "differr", x = tools)] <- "DiffErr"
  tools[grep(pattern = "eligos", x = tools)] <- "ELIGOS"
  tools[grep(pattern = "epinanoErr", x = tools)] <- "EpiNano-Error"
  tools[grep(pattern = "epinanoSvm", x = tools)] <- "EpiNano-SVM"
  tools[grep(pattern = "m6anet", x = tools)] <- "m6Anet"
  tools[grep(pattern = "mines", x = tools)] <- "MINES"
  tools[grep(pattern = "nanocompore", x = tools)] <- "Nanocompore"
  tools[grep(pattern = "nanom6a", x = tools)] <- gsub(x = gsub(pattern = "nanom6a_output\\.bed", replacement = "Nanom6A", x = tools[grep(pattern = "nanom6a", x = tools)]), pattern = "\\.tsv\\.bed", replacement = "")
  tools[grep(pattern = "tomboComparison", x = tools)] <- "Tombo"
  tools[grep(pattern = "xpore", x = tools)] <- "xPore"
  tools[grep(pattern = "yanocomp", x = tools)] <- "Yanocomp"
  tools[grep(pattern = "nanodoc", x = tools)] <- "NanoDoc"
  # Build matrix of zeros
  hitsMatrix <- matrix(c(0), nrow = length(genesBins_par), ncol = length(files_par) + 1)
  colnames(hitsMatrix) <- c("Reference_set", tools)
  
  row.names(hitsMatrix) <- c(1:nrow(hitsMatrix))
  
  peaks_overlap <- findOverlaps(query = peaks_par, subject = genesBins_par, minoverlap = 1, type = "any")
  hitsMatrix[unique(subjectHits(peaks_overlap)), "Reference_set"] <- 1
  overlapMatrix <- hitsMatrix
  ind_nodef <- setdiff(grep(x = colnames(overlapMatrix), pattern = "Nanom6A"), grep(x = colnames(overlapMatrix), pattern = "Nanom6A_ratio\\.0\\.5"))
  if (length(ind_nodef) > 0) {
    overlapMatrix <- overlapMatrix[, -ind_nodef]
  }
  
  matrix_nanom6A <- hitsMatrix[,"Reference_set"]
  names_nanom6A <- c()
  
  Performances <- list()
  listPRcurves <- list()
  
  for (y in files) {
    ind <- which(files == y)
    x <- tools[ind]
    
    cat(sprintf("Processing file: %s\n", y))
    recall <- c()
    precision <- c()
    # Extraction of bed file + conversion to granges
    bed_file <- read.table(y, header = T, sep = "\t")
    granges <- makeGRangesFromDataFrame(bed_file, keep.extra.columns = T)
    # Overlap between m6A detected site of each tool and genome binned
    overlap <- as.matrix(findOverlaps(query = granges, subject = genesBins_par, minoverlap = 1, type = "any"))
    if (grepl(x, pattern = paste0(c("DENA","DRUMMER","DiffErr","Yanocomp","Nanocompore","ELIGOS","EpiNano-Error",
                                    "EpiNano-SVM","xPore","NanoDoc","Tombo", "m6Anet"), collapse = "|"))) {
      # Add column of filtering parameter
      filtering_parameter <- bed_file[overlap[,"queryHits"] , 6]
      overlap_w_parameter <- cbind(overlap, filtering_parameter)
      default <- unname(threshold_default[grep(x, pattern = paste0(c("DENA","DRUMMER","DiffErr","Yanocomp","Nanocompore","ELIGOS","EpiNano-Error","EpiNano-SVM",
                                                                     "xPore","NanoDoc","Tombo", "m6Anet"), collapse = "|"), value = T)])
      # Recognize from the name of the tool if we need to keep the maximum or minimum value (when there are more hits in a single bin)
      if(grepl(x, pattern = listmax)) { 
        score <- sapply(split(overlap_w_parameter[,3], overlap_w_parameter[,2]), max)
        default_thr <- default
        hitsMatrix[, x] <- rep(0, length(genesBins_par))
        overlapMatrix[, x] <- rep(0, length(genesBins_par))
      } else {
        score <- -1*sapply(split(overlap_w_parameter[,3], overlap_w_parameter[,2]), min)
        default_thr <- - default
        hitsMatrix[, x] <- rep(-1, length(genesBins_par))
        overlapMatrix[, x] <- rep(0, length(genesBins_par))
      }
      
      #assign score to hitsMatrix
      hitsMatrix[names(score), x] <- score
      
      #assign 0 or 1 value to overlapMatrix for undetected/detected peaks at default values
      pred_pos_def <- which(hitsMatrix[, x] >= default_thr)
      pred_neg_def <- which(hitsMatrix[, x] < default_thr)
      overlapMatrix[pred_pos_def, x] <- 1
      overlapMatrix[pred_neg_def, x] <- 0
      
      positive <- hitsMatrix[which(hitsMatrix[,"Reference_set"] == 1), x]
      negative <- hitsMatrix[which(hitsMatrix[,"Reference_set"] == 0), x]
      #par(mfrow = c(2, 1))
      #pdf(file = paste0(x,"_scores_distribution.pdf"))
      #hist(positive, main = paste0(x, " - Scores for positive peaks"))
      #hist(negative, main = paste0(x, " - Scores for negative peaks"))
      #dev.off()
      if (length(negative) > 0) {
        pr <- pr.curve(scores.class0 = unname(positive), scores.class1 = unname(negative), curve=T, rand.compute=TRUE)
        save(pr, file = paste0(resultsFolder, "/", x,"_PRcurve", notes, "_window_", w, "bp.Rdata"))
        pdf(file = paste0(resultsFolder, "/", x,"_PRcurve", notes, "_window_", w, "bp.pdf"), width = 8, height = 8)
        plot(pr, main = paste0(x, " Precision-Recall curve"), rand.plot=TRUE)
        dev.off()
        listPRcurves[[x]] <- pr
      } else {
        cat(sprintf("All genome bins include peaks, skipping PR curve plotting for file %s\n", x))
      }
      # Plot "manual" PR curve
      thresholds <- c(seq(from = min(score), to = max(score), length.out = 100), default_thr)
      recall <- seq(from = 0, to = 0, length.out = 101)
      precision <- seq(from = 0, to = 0, length.out = 101)
      F1_score <- seq(from = 0, to = 0, length.out = 101)
      for (t in 1:length(thresholds)) {
        thr <- thresholds[t]
        TP <- length(which(positive >= thr))
        FN <- length(which(positive < thr))
        FP <- length(which(negative >= thr))
        TN <- length(which(negative < thr))
        recall[t] <- TP/(TP + FN)
        precision[t] <- TP/(TP + FP)
        F1_score[t] <- 2*recall[t]*precision[t]/(recall[t] + precision[t])
      }
      Performances[[x]] <- data.frame(tool = x, threshold = thresholds, recall = recall, precision = precision, F1_score = F1_score)      
      #pdf(paste0(x, "_PRcurve_manual.pdf"))
      #plot(recall, precision, main = paste0(x, " - PR manual"), type = "l", xlim = c(0, 1), ylim = c(0, 1))
      #dev.off()
    } 
    else if (grepl(x, pattern = "MINES")) {
      hitsMatrix[overlap[, 2], x] <- 1
      overlapMatrix[overlap[, 2], x] <- 1
      TP <- 0
      for (y in 1:nrow(hitsMatrix)){
        if ((hitsMatrix[y, "Reference_set"] == 1) && (hitsMatrix[y, x] == 1)){
          TP <- TP + 1
        }
      }
      # Recall + Precision
      totPositiveGS <- length(which(hitsMatrix[, "Reference_set"] == 1))
      totPositiveTool <- length(which(hitsMatrix[, x] == 1))
      recall <- TP/totPositiveGS
      precision <- TP/totPositiveTool
      # F1 score
      F1_score <- 2*(precision*recall)/(precision+recall)
      Performances[[x]]	<- data.frame(tool = x, threshold = thresholds, recall = recall, precision = precision, F1_score = F1_score)
    } 
    else if (grepl(x, pattern = "Nanom6A")) {
      hitsMatrix[overlap[, 2], x] <- 1
      matrix_nanom6A <- cbind(matrix_nanom6A, hitsMatrix[, x])
      names_nanom6A <- c(names_nanom6A, x)
      if (grepl(x, pattern = "0\\.5")){
        overlapMatrix[overlap[, 2], x] <- 1
        TP <- 0
        for (y in 1:nrow(hitsMatrix)){
          if ((hitsMatrix[y, "Reference_set"] == 1) && (hitsMatrix[y, x] == 1)){
            TP <- TP + 1
          }
        }
        # Recall + Precision
        totPositiveGS <- length(which(hitsMatrix[, "Reference_set"] == 1))
        totPositiveTool <- length(which(hitsMatrix[, x] == 1))
        recall <- TP/totPositiveGS
        precision <- TP/totPositiveTool
        # F1 score
        F1_score <- 2*(precision*recall)/(precision+recall)
        Performances[[x]] <- data.frame(tool = x, threshold = "0.5", recall = recall, precision = precision, F1_score = F1_score)
      }
    }
  }
  
  if (length(grep(x = tools, pattern = "Nanom6A")) > 0) {
    ### Code for PR curve nanom6A which has multiple files each run with a different threshold
    names_nanom6A <- str_extract(names_nanom6A, "[0-9]\\.[0-9]*")
    colnames(matrix_nanom6A) <- c("RS", names_nanom6A)
    
    max_thr <- c()
    short <- matrix_nanom6A[,2:ncol(matrix_nanom6A)]
    if (length(which(short == 0)) > 0) {
      max_thr <- vector(length = nrow(short), mode = "numeric")
      names(max_thr) <- rownames(short)
      
      short_nozero <- short[names(which(apply(short, 1, function(x) any(as.logical(x)) != 0 ))), ]
      tmp <- apply(short_nozero, 1, function(x) {names(which(x == 1))[length(which(x == 1))]})
      max_thr[names(tmp)] <- as.numeric(tmp)
      max_thr <- unname(max_thr)
      
      new_matrix <- cbind(matrix_nanom6A, max_thr)
      posit <- new_matrix[which(new_matrix[, "RS"] == 1), "max_thr"]
      negat <- new_matrix[which(new_matrix[, "RS"] == 0), "max_thr"]
      
      pr <- pr.curve(posit, negat, curve = T, rand.compute=TRUE)
      pdf(file = paste0(resultsFolder, "/Nanom6A_PRcurve", notes, "_window_", w, "bp.pdf"), width = 8, height = 8)
      plot(pr, main = "Nanom6A Precision-Recall curve", rand.plot=TRUE)
      dev.off()
      listPRcurves[["Nanom6A"]] <- pr      
      #par(mfrow = c(2, 1))
      #pdf(file = paste0(resultsFolder, "nanom6A_scores_distribution.pdf"))
      #hist(posit, 100,main = "nanom6A - Scores for positive peaks")
      #hist(negat, 100,main = "nanom6A - Scores for negative peaks")
      #dev.off()
    } else {
      cat("All genome bins include peaks, skipping PR curve plotting for tool nanom6a\n")
    }
  }
  
  sink(paste0(resultsFolder, "/Performances", notes, "_window_", w, "bp.tsv"))
  print(Performances)
  sink()
  
  if (length(negative) > 0) {
    ### Plot all the Precision-Recall curves together
    col <- c(7,8,420,153,31,100,33,47,53,62,400,454,28,10)
    
    pdf(file = paste0(resultsFolder, "/Summary_PR_curves", notes, "_window_", w, "bp.pdf"), width = 8, height = 8)
    for (x in 1:length(listPRcurves)) {
      if (x == 1){
        plot(listPRcurves[[x]], color = colors()[col[x]], main = "Summary PR curves", rand.plot=TRUE)
      }
      else{
        plot(listPRcurves[[x]], add = T, color = colors()[col[x]], main = "Summary PR curves")
      }
    }
    legend("bottomright", legend = names(listPRcurves), col = colors()[col[1:length(listPRcurves)]], lty=1:1, cex=0.8, bg = "lightblue" )
    dev.off()
  } else {
    cat("All genome bins include peaks, skipping summary PR curve plotting\n")
  }
  
  colnames(overlapMatrix) <- gsub(pattern = "_output.*", replacement = "", x = colnames(overlapMatrix))
  
  data_ovlp <- matrix(data = 0, nrow = dim(overlapMatrix)[2] - 1, ncol = dim(overlapMatrix)[2] - 1)
  tools <- setdiff(colnames(overlapMatrix), "Reference_set")

  colnames(data_ovlp) <- tools
  rownames(data_ovlp) <- tools
  
  for (i in tools) {
    for (j in tools) {
      data_ovlp[i, j] <- length(intersect(which(overlapMatrix[, i] == 1), which(overlapMatrix[, j] == 1)))/length(which(overlapMatrix[, i] == 1))
    }
  }
  pdf(paste0(resultsFolder, "/Tools_overlap_default_par", notes, "_window_", w, "bp.pdf"))
  pheatmap(data_ovlp, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = TRUE, fontsize = 15, display_numbers = FALSE, color = colorRampPalette(c("white", "red"))(30))
  dev.off()
  return(list(hitsMatrix, overlapMatrix, Performances, data_ovlp))
}

results <- Run_statistical_analysis(genesBins_par = genesBins, peaks_par = peaks_granges, files_par = files, notes = "", w = w)
hitsMatrix <- results[[1]]
overlapMatrix <- results[[2]]
Performances <- results[[3]]
data_ovlp <- results[[4]]
save(results, file = paste0(resultsFolder, "/Results_window_", w, "bp.rda"))

results_RRACH <- Run_statistical_analysis(genesBins_par = genesBins_RRACH, peaks_par = peaks_RRACH_granges, files_par = files, notes = "_RRACH", w = w)
hitsMatrix_RRACH <- results_RRACH[[1]]
overlapMatrix <- results_RRACH[[2]]
Performances_RRACH <- results_RRACH[[3]]
data_ovlp_RRACH <- results_RRACH[[4]]
save(results_RRACH, file = paste0(resultsFolder, "/Results_window_", w, "bp_RRACH.rda"))

if (!is.null(highcov_GRanges)) {
  results_highcov <- Run_statistical_analysis(genesBins_par = genesBins_highcov, peaks_par = peaks_highcov_granges, files_par = files, notes = "_highcov", w = w)
  hitsMatrix_highcov <- results_highcov[[1]]
  overlapMatrix_highcov <- results_highcov[[2]]
  Performances_highcov <- results_highcov[[3]]
  data_ovlp_highcov <- results_highcov[[4]]
  save(results_highcov, file = paste0(resultsFolder, "/Results_window_", w, "bp_highcov.rda"))
  
  results_highcov_RRACH <- Run_statistical_analysis(genesBins_par = genesBins_highcov_RRACH, peaks_par = peaks_highcov_RRACH_granges, files_par = files, notes = "_highcov_RRACH", w = w)
  hitsMatrix_highcov_RRACH <- results_highcov_RRACH[[1]]
  overlapMatrix_highcov_RRACH <- results_highcov_RRACH[[2]]
  Performances_highcov_RRACH <- results_highcov_RRACH[[3]]
  data_ovlp_highcov_RRACH <- results_highcov_RRACH[[4]]
  save(results_highcov_RRACH, file = paste0(resultsFolder, "/Results_window_", w, "bp_highcov_RRACH.rda"))
}
