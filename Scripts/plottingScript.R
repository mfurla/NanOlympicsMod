#############################
#packages and global parameters

#load packages
library("GenomicFeatures")
library("Guitar")
library("GenomicRanges")
library("GenomicAlignments")
library("parallel")
library("stringr")
library("Biostrings")
library("BSgenome")
library("ggplot2")
library("memes")
library("RColorBrewer")
library("ggpattern")
library("RColorBrewer")
options(scipen = 1000)

#choose dataset
#dataset <- "oligos"
#dataset <- "yeast"
#dataset <- "mESC"
dataset <- "HEK293T"

#proper_bed_flag <- 0 if the file is in format: chr\tstart\tend\tstrand\tstatus\tscore
#proper_bed_flag <- 1 if the file is in format: chr\tstart\tend\tname\tscore\tstrand
proper_bed_flag <- 0

#############################
#hits number and metagene plots

#resize each hit for metagene plot
tol <- 100
#downsample if there are more than max_filt_sites hits
#(required because potentially very time consuming)
max_filt_sites <- 100000

if (length(grep(x = dataset, pattern = "yeast")) > 0) {
  #read gtf and extract genes
  gtf_file <- "/path/to/sk1_cds_final.V1.1.gtf"
  input_dir <- "/path/to/output_bed_files/"
  gold_standard_peaks_file <- "/path/to/Schwartz_m6a_matzer_to_sk1_sorted.bed"
} else if (length(grep(x = dataset, pattern = "mESC")) > 0) {
  #read gtf and extract genes
  gtf_file <- "/path/to/Mus_musculus.GRCm38.102.gtf"
  input_dir <- "/path/to/output_bed_files/"
  gold_standard_peaks_file <- "/path/to/mESC_miCLIP2_OR_GLORI_intersection.bed"
} else if (length(grep(x = dataset, pattern = "oligos")) > 0) {
  gtf_file <- "/path/to/oligos.V1.1.gtf"
  input_dir <- "/path/to/output_bed_files/"
  gold_standard_peaks_file <- "/path/to/Oligos_peaks.bed"
} else if (length(grep(x = dataset, pattern = "HEK293T")) > 0) {
  gtf_file <- "/path/to/Homo_sapiens.GRCh38.109.gtf"
  input_dir <- "/path/to/output_bed_files/"
  gold_standard_peaks_file <- "/path/to/GLORI_reference_set_chr1.bed"
}

#create txdb from GTF file
txdb <- makeTxDbFromGFF(file = gtf_file, format = "auto")

#import bed files
bed_files <- list.files(path = input_dir, pattern = "\\.bed", full.names = TRUE)
bed_files <- bed_files[-setdiff(grep(pattern = "nanom6a", x = bed_files), grep(pattern = "nanom6a_output\\.bed_ratio\\.0\\.5", x = bed_files))]

#tools with confidence parameter that needs to be maximised
listmax <- paste0(c("DENA", "EpiNano-Error", "EpiNano-SVM", "NanoDoc", "m6Anet"), collapse = "|") # for these tools we need to maximize the filtering parameter
#tools with confidence parameter that needs to be minimised
listmin <- paste0(c("DiffErr", "DRUMMER", "Yanocomp", "Nanocompore", "ELIGOS", "xPore", "Tombo"), collapse = "|") # for these tools we need to minimize the filtering parameter

#default parameter threshold
threshold_default <- c(0.1, 0.05, 0.01, 0.05, 0.01, 0.0001, 0.1, 0.5, 0.05, 0.02, 0.05, 0.9)
names(threshold_default) <- c("DENA", "DiffErr", "DRUMMER", "Yanocomp", "Nanocompore", "ELIGOS", 
                              "EpiNano-Error", "EpiNano-SVM", "xPore", "NanoDoc", "Tombo", "m6Anet")

#tools' names formatting
tools <- unique(gsub(pattern = "_output.*", replacement = "", x = basename(bed_files)))
tools[grep(pattern = "dena", x = tools)] <- "DENA"
tools[grep(pattern = "drummer", x = tools)] <- "DRUMMER"
tools[grep(pattern = "differr", x = tools)] <- "DiffErr"
tools[grep(pattern = "eligos", x = tools)] <- "ELIGOS"
tools[grep(pattern = "epinanoErr", x = tools)] <- "EpiNano-Error"
tools[grep(pattern = "epinanoSvm", x = tools)] <- "EpiNano-SVM"
tools[grep(pattern = "m6anet", x = tools)] <- "m6Anet"
tools[grep(pattern = "mines", x = tools)] <- "MINES"
tools[grep(pattern = "nanocompore", x = tools)] <- "Nanocompore"
tools[grep(pattern = "nanom6a", x = tools)] <- "Nanom6A"
tools[grep(pattern = "tomboComparison", x = tools)] <- "Tombo"
tools[grep(pattern = "xpore", x = tools)] <- "xPore"
tools[grep(pattern = "yanocomp", x = tools)] <- "Yanocomp"
tools[grep(pattern = "nanodoc", x = tools)] <- "NanoDoc"
tools_par <- data.frame(tool = tools, par_optim = ifelse(grepl(pattern = listmax, x = tools), "max", "min"), thr_default = threshold_default[tools])
rownames(tools_par) <- tools_par$tool

#read reference-set peaks file
gs <- read.table(file = gold_standard_peaks_file, header = TRUE, sep = "\t")
colnames(gs) <- c("chr", "start", "end", "id", "score", "strand")
if (length(grep(x = dataset, pattern = "mESC")) > 0) {
  gs[, 5] <- unlist(lapply(strsplit(gs[, 5], ","), function(x) mean(as.numeric(x))))
  gs[, 5] <- as.numeric(sprintf("%.5f", gs[, 5]))
}

#downsample reference set peaks and write file in bed format
max_filt_sites_gs <- 10000
if (dim(gs)[1] > max_filt_sites_gs) {
  cat(sprintf("Downsampling %d peaks for bed file: %s\n", max_filt_sites_gs, basename(gold_standard_peaks_file)))
  idx <- sample(x = 1:dim(gs)[1], size = max_filt_sites_gs)
  gs <- gs[idx, ]
}
tmp_file_gs <- paste0("/tmp/", dataset, "_gold_standard_peaks.bed")
write.table(x = gs, file = tmp_file_gs, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

plot_type <- c("tx", "mrna", "tx", "mrna")
names(plot_type) <- c("yeast", "mESC", "oligos", "HEK293T")

#plot reference set metagene
pdf(paste0("/path/to/", dataset, "_gold_standard_peaks_metagene.pdf"))
gs_file_list <- as.list(tmp_file_gs)
gs_file_list_plat <- c("refSet_m6A", "refSet_MAZTER-seq_m6A-seq", "refSet_miCLIP2_GLORI", "refSet_GLORI")
names(gs_file_list_plat) <- c("oligos", "yeast", "mESC", "HEK293T")
names(gs_file_list) <- gs_file_list_plat[dataset]

GuitarPlot(txTxdb = txdb , stBedFiles = gs_file_list, stGroupName = names(gs_file_list), pltTxType = plot_type[dataset])
dev.off()

#import genes as GRanges
genes_txdb <- GenomicFeatures::genes(txdb)
#import reference set as GRanges
gs_GRanges <- GRanges(seqnames = gs[, 1],
                      ranges = IRanges(start = gs[, 2], end = gs[, 3]),
                      strand = gs[, 6])
gs_hits <- findOverlaps(query = gs_GRanges, subject = genes_txdb)
gs_genes <- unique(names(genes_txdb[subjectHits(gs_hits)]))
length(gs_genes)

#plot the number of hits for each tool and corresponding confidence parameter threshold
pdf(paste0("/path/to/", dataset, "_peaks_summary.pdf"))
par(mfrow = c(2, 2))
tmp_file_list <- c()
for (i in 1:length(bed_files)) {
  bed_file_curr <- bed_files[i]
  #tool_curr <- gsub(pattern = "_output.*", replacement = "", x = basename(bed_file_curr))
  tool_curr <- tools[i]
  cat(sprintf("Processing bed file: %s\n", basename(bed_files[i])))
  bed_curr <- read.table(file = bed_file_curr, header = TRUE, sep = "\t")
  if (proper_bed_flag == 1) {
    colnames(bed_curr) <- c("Chr", "Start", "End", "Id", "Score", "Strand")
    peaks_size <- bed_curr$End - bed_curr$Start
    plot_score_flag <- length(which(!is.na(bed_curr[, 5])) > 0) > 0
  } else {
    peaks_size <- bed_curr[, 3] - bed_curr[, 2]
    plot_score_flag <- dim(bed_curr)[2] > 5
  }
  if (plot_score_flag) {
    if (proper_bed_flag == 1) {
      score <- bed_curr[, 5]
      score_name <- colnames(bed_curr)[5] 
    } else {
      score <- bed_curr[, 6]
      score_name <- colnames(bed_curr)[6]
    }
    num_tot_peaks <- length(score)
    if (grepl(x = tools_par[tool_curr, "par_optim"], pattern = "max")) {
      ind_filtered_peaks <- which(score > tools_par[tool_curr, "thr_default"])
    } else {
      ind_filtered_peaks <- which(score < tools_par[tool_curr, "thr_default"])
    }
    num_filt_peaks <- length(ind_filtered_peaks)
    hist(peaks_size, main = paste0("Tool: ", tool_curr, "\nPeaks size (nt)"), breaks = 100,  cex.main = 0.7)
    hist(score, main = paste0("Tool: ", tool_curr, "\n", score_name, " distribution\nNum. tot. peaks = ", num_tot_peaks, "\nNum. filt. peaks = ", num_filt_peaks), breaks = 100, cex.main = 0.7, xlim = c(0, 1))
    col_curr <- ifelse(grepl(x = tools_par[tool_curr, "par_optim"], pattern = "max"), 2, 3)
    abline(v = tools_par[tool_curr, "thr_default"], col = col_curr, lty = 3, lwd = 3)
    bed_curr_filt <- bed_curr[ind_filtered_peaks, ]
  } else {
    num_filt_peaks <- length(peaks_size)
    hist(peaks_size, main = paste0("Tool: ", tool_curr, "\nPeaks size (nt)\nNum. filt. peaks = ", num_filt_peaks), breaks = 100, cex.main = 0.7)
    bed_curr_filt <- cbind(bed_curr, ".")
    plot.new()
  }
  chr_names <- levels(seqnames(genes_txdb))
  ind_chr_ok <- which(bed_curr_filt[, 1] %in% chr_names)
  bed_curr_filt <- bed_curr_filt[ind_chr_ok, ]
  tmp_file <- paste0("/tmp/", dataset, "_", tool_curr, "_peaks.bed")
  if (dim(bed_curr_filt)[1] > max_filt_sites) {
    cat(sprintf("Downsampling %d peaks for bed file: %s\n", max_filt_sites, basename(bed_files[i])))
    idx <- sample(x = 1:dim(bed_curr_filt)[1], size = max_filt_sites)
    bed_curr_filt <- bed_curr_filt[idx, ]
  }
  if (proper_bed_flag == 1) {
    write.table(x = cbind(bed_curr_filt[, 1], bed_curr_filt[, 2] - tol, bed_curr_filt[, 3] + tol, ".", bed_curr_filt[, 5], bed_curr_filt[, 4]), file = tmp_file, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    tool_GRanges <- GRanges(seqnames = bed_curr_filt[, 1],
                            ranges = IRanges(start = bed_curr_filt[, 2], end = bed_curr_filt[, 3]),
                            strand = bed_curr_filt[, 6])
  } else {
    write.table(x = cbind(bed_curr_filt[, 1], bed_curr_filt[, 2] - tol, bed_curr_filt[, 3] + tol, ".", bed_curr_filt[, 6], bed_curr_filt[, 4]), file = tmp_file, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    tool_GRanges <- GRanges(seqnames = bed_curr_filt[, 1],
                            ranges = IRanges(start = bed_curr_filt[, 2], end = bed_curr_filt[, 3]),
                            strand = bed_curr_filt[, 4])
  }
  tool_hits <- findOverlaps(query = tool_GRanges, subject = genes_txdb)
  tool_genes <- unique(names(genes_txdb[subjectHits(tool_hits)]))
  
  num_genes_tool <- length(tool_genes)
  recall_tool <- length(intersect(tool_genes, gs_genes))/length(gs_genes)
  precision_tool <- length(intersect(tool_genes, gs_genes))/length(tool_genes)
  F1_tool <- 2*recall_tool*precision_tool/(precision_tool + recall_tool)
  cat(sprintf("Detected m6A in %d genes\nRecall=%.3f Precision=%.3f F1score=%.3f\n\n", num_genes_tool, recall_tool, precision_tool, F1_tool))
  
  tmp_file_list <- c(tmp_file_list, tmp_file)
}
dev.off()

tmp_file_list <- as.list(tmp_file_list)
names(tmp_file_list) <- tools

#plot m6A hits metagene with Guitar package
pdf(paste0("/path/to/", dataset, "_peaks_metagene_noRS.pdf"))
p1 <- list()
p2 <- list()
p3 <- list()

for (i in 1:length(bed_files)) {
  tool_curr <- tools[i]
  cat(sprintf("Processing bed file: %s\n", basename(bed_files[i])))
  tmp_file_list_curr <- as.list(tmp_file_list[[i]])
  p1[[i]] <- GuitarPlot(txTxdb = txdb , stBedFiles = tmp_file_list_curr, stGroupName = tool_curr, pltTxType = plot_type[dataset])
}
txpromoterLength_val = c(1, 1000, 1000, 1000)
txtailLength_val = c(1, 1000, 1000, 1000)
names(txpromoterLength_val) <- c("oligos", "yeast", "mESC", "HEK293T")
names(txtailLength_val) <- c("oligos", "yeast", "mESC", "HEK293T")

p2[[1]] <- GuitarPlot(txTxdb = txdb , stBedFiles = c(gs_file_list, tmp_file_list), stGroupName = c(names(gs_file_list), names(tmp_file_list)), pltTxType = plot_type[dataset], enableCI = FALSE, txpromoterLength = unname(txpromoterLength_val[dataset]), txtailLength = unname(txtailLength_val[dataset]))
p3[[1]] <- GuitarPlot(txTxdb = txdb , stBedFiles = tmp_file_list, stGroupName = names(tmp_file_list), pltTxType = plot_type[dataset], enableCI = FALSE, txpromoterLength = unname(txpromoterLength_val[dataset]), txtailLength = unname(txtailLength_val[dataset]))
print(p1)
print(p2)
print(p3)
dev.off()

#plot m6A hits metagene with Rbase using data from Guitar package
data_split <- split(p3[[1]]$data[c("x", "density")], p3[[1]]$data$group)
colors_MP <- c(brewer.pal(n = 12, name = 'Paired'))
colors_MP[11] <- "gold"
tools <- c("Nanocompore", "Tombo", "xPore", "Yanocomp", "ELIGOS", "EpiNano-Error", "DiffErr", "DRUMMER", "m6Anet", "DENA", "Nanom6A", "EpiNano-SVM", "MINES")
colors_MP <- c(colors_MP, "black")
names(colors_MP) <- tools

maxval <- max(unlist(lapply(data_split, '[[', 2)))

pdf(paste0("/path/to/", dataset, "_metagene_plot.pdf"))
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
for (i in match(intersect(tools, names(data_split)), names(data_split))) {
  if (i == match(intersect(tools, names(data_split)), names(data_split))[1]) {
    plot(data_split[[names(data_split)[i]]][, 1], data_split[[names(data_split)[i]]][, 2], type = 'l', lty = 1, lwd = 4, col = colors_MP[names(data_split[i])], xlab = "Metagene", ylab = "Density", xlim = c(0, 1), ylim = c(0, maxval), cex.main = 0.7)
  } else {
    points(data_split[[names(data_split)[i]]][, 1], data_split[[names(data_split)[i]]][, 2], type = "l", lty = 1, lwd = 4, col = colors_MP[names(data_split[i])])
  }
  legend("topright", inset=c(-0.2, 0), legend = names(data_split)[match(intersect(tools, names(data_split)), names(data_split))], col = c(colors_MP[names(data_split)[match(intersect(tools, names(data_split)), names(data_split))]], col = colors_MP[match(intersect(tools, names(data_split)), names(data_split))]), lty = 1, cex = 0.5, lwd = 4)
}
dev.off()

#############################
#plot barplots with the number of hits for each tool

num_filt_peaks_oligos <- vector(mode = "numeric", length = length(tools))
names(num_filt_peaks_oligos) <- tools
num_filt_peaks_oligos["DiffErr"] <- 10275 #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_oligos["DRUMMER"] <- 757  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_oligos["ELIGOS"] <- 3972  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_oligos["EpiNano-Error"] <- 305  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_oligos["EpiNano-SVM"] <- 135  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_oligos["m6Anet"] <- 67  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_oligos["MINES"] <- 38  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_oligos["Nanocompore"] <- 4742   #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_oligos["Nanom6A"] <- 130  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_oligos["Tombo"] <- 9989  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_oligos["xPore"] <- 3785  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_oligos["Yanocomp"] <- 6243  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.

num_filt_peaks_yeast <- vector(mode = "numeric", length = length(tools))
names(num_filt_peaks_yeast) <- tools
num_filt_peaks_yeast["DiffErr"] <- 1421  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_yeast["ELIGOS"] <- 391  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_yeast["m6Anet"] <- 171  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_yeast["MINES"] <- 3623  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_yeast["Nanocompore"] <- 2011  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_yeast["Nanom6A"] <- 64656  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_yeast["Tombo"] <- 1080239  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_yeast["xPore"] <- 72760  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_yeast["Yanocomp"] <- 3888  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.

num_filt_peaks_mESC <- vector(mode = "numeric", length = length(tools))
names(num_filt_peaks_mESC) <- tools
num_filt_peaks_mESC["DENA"] <- 144336  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_mESC["DiffErr"] <- 2178  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_mESC["ELIGOS"] <- 1665  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_mESC["EpiNano-Error"] <- 70639  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_mESC["m6Anet"] <- 8218  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_mESC["MINES"] <- 38131  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_mESC["Nanocompore"] <-  4348  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_mESC["Nanom6A"] <- 40236  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_mESC["Tombo"] <- 792700  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_mESC["xPore"] <- 168813  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_mESC["Yanocomp"] <- 47851  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.

num_filt_peaks_HEK293T <- vector(mode = "numeric", length = length(tools))
names(num_filt_peaks_HEK293T) <- tools
num_filt_peaks_HEK293T["DENA"] <- 45281  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_HEK293T["DiffErr"] <- 235104  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_HEK293T["ELIGOS"] <- 1114  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_HEK293T["EpiNano-SVM"] <- 36403  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_HEK293T["EpiNano-Error"] <- NA  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_HEK293T["m6Anet"] <- 2627  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_HEK293T["MINES"] <- 17772  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_HEK293T["Nanocompore"] <-  756  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_HEK293T["Nanom6A"] <- 25931  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_HEK293T["Tombo"] <- 109387  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_HEK293T["xPore"] <- 72127  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
num_filt_peaks_HEK293T["Yanocomp"] <- 4700  #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.

num_hits_df <- data.frame(Tool = c(names(num_filt_peaks_oligos), names(num_filt_peaks_yeast), names(num_filt_peaks_mESC), names(num_filt_peaks_HEK293T)), Dataset = c(rep("Oligos", length(num_filt_peaks_oligos)), rep("Yeast", length(num_filt_peaks_yeast)), rep("mESC", length(num_filt_peaks_mESC)), rep("HEK293T", length(num_filt_peaks_HEK293T))), Num_hits = c(num_filt_peaks_oligos, num_filt_peaks_yeast, num_filt_peaks_mESC, num_filt_peaks_HEK293T))
num_hits_df <- num_hits_df[which(num_hits_df$Num_hits != 0), ]
num_hits_df_oligos <- num_hits_df[which(num_hits_df$Dataset == "Oligos"), ]
num_hits_df_oligos <- num_hits_df_oligos[sort(num_hits_df_oligos$Num_hits, decreasing = TRUE, index.return = TRUE)$ix, ]
num_hits_df_oligos$Tool <- factor(num_hits_df_oligos$Tool, levels = as.character(num_hits_df_oligos$Tool))

num_hits_df_yeast <- num_hits_df[which(num_hits_df$Dataset == "Yeast"), ]
num_hits_df_yeast <- num_hits_df_yeast[sort(num_hits_df_yeast$Num_hits, decreasing = TRUE, index.return = TRUE)$ix, ]
num_hits_df_yeast$Tool <- factor(num_hits_df_yeast$Tool, levels = as.character(num_hits_df_yeast$Tool))

num_hits_df_mESC <- num_hits_df[which(num_hits_df$Dataset == "mESC"), ]
num_hits_df_mESC <- num_hits_df_mESC[sort(num_hits_df_mESC$Num_hits, decreasing = TRUE, index.return = TRUE)$ix, ]
num_hits_df_mESC$Tool <- factor(num_hits_df_mESC$Tool, levels = as.character(num_hits_df_mESC$Tool))

num_hits_df_HEK293T <- num_hits_df[which(num_hits_df$Dataset == "HEK293T"), ]
num_hits_df_HEK293T <- num_hits_df_HEK293T[sort(num_hits_df_HEK293T$Num_hits, decreasing = TRUE, index.return = TRUE)$ix, ]
num_hits_df_HEK293T$Tool <- factor(num_hits_df_HEK293T$Tool, levels = as.character(num_hits_df_HEK293T$Tool))

ggplot(num_hits_df_oligos, aes(x = Tool, y = Num_hits, fill = Dataset, group = Dataset)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x="Tool", y = "Num. hits") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        panel.background = element_blank()) + 
  scale_y_log10()
ggsave("path/to/Num_hits_default_oligos.pdf")

ggplot(num_hits_df_yeast, aes(x = Tool, y = Num_hits, fill = Dataset, group = Dataset)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x="Tool", y = "Num. hits") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        panel.background = element_blank()) + 
  scale_y_log10()
ggsave("/path/to/Num_hits_default_yeast.pdf")

ggplot(num_hits_df_mESC, aes(x = Tool, y = Num_hits, fill = Dataset, group = Dataset)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x="Tool", y = "Num. hits") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        panel.background = element_blank()) + 
  scale_y_log10()
ggsave("/path/to/Num_hits_default_mESC.pdf")

ggplot(num_hits_df_HEK293T, aes(x = Tool, y = Num_hits, fill = Dataset, group = Dataset)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x="Tool", y = "Num. hits") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        panel.background = element_blank()) + 
  scale_y_log10()
ggsave("/path/to/Num_hits_default_HEK293T.pdf")

#############
#RRACH accuracy for mESC dataset

#read gtf and extract genes
gtf_file <- "/path/to/Mus_musculus.GRCm38.102.gtf"
input_dir <- "/path/to/output_bed_files/"
gold_standard_peaks_file <- "/path/to/mESC_miCLIP2_OR_GLORI_intersection.bed"
genome_file <- "/path/to/GRCm38.p6.genome.fasta"
genesbed <- "/path/to/gencode.vM23.genes.bed6"
binLength <- 50
mccores <- 1

#import gtf file as txdb
txdb <- makeTxDbFromGFF(file = gtf_file, format = "auto")

#read reference-set bed file
gs <- read.table(file = gold_standard_peaks_file, header = FALSE, sep = "\t")
gs[, 5] <- unlist(lapply(strsplit(gs[, 5], ","), function(x) mean(as.numeric(x))))
colnames(gs) <- c("chr", "start", "end", "name", "score", "strand")

#convert reference-set to GRanges
gs_GRanges <- GRanges(seqnames = gs[, 1],
                      ranges = IRanges(start = gs[, 3], end = gs[, 3]),
                      strand = gs[, 6])

#load results at default value
load("/path/to/Results_window_50bp.rda")
hitsMatrix <- results[[2]]

#read genome file and find RRACH motifs
chrs <- readDNAStringSet(genome_file, format="fasta")

RRACH_plus <- GRanges(vmatchPattern(pattern = "RRACH", subject = chrs, fixed = "subject"), strand = "+")
RRACH_minus <- GRanges(vmatchPattern(pattern = "DGTYY", subject = chrs, fixed = "subject"), strand = "-")

RRACH <- c(RRACH_plus, RRACH_minus)
RRACH_seq <- as.character(getSeq(chrs, as(RRACH, "GRanges")))
mcols(RRACH) <- data.frame(motif = RRACH_seq)

#find overlaps between reference-set and RRACH motifs
hits_gs_RRACH <- findOverlaps(query = gs_GRanges, subject = RRACH)
gs_RRACH <- gs_GRanges[queryHits(hits_gs_RRACH)]
mcols(gs_RRACH) <- data.frame(motif = mcols(RRACH)$motif[subjectHits(hits_gs_RRACH)])

#bin the genome in windows of length w
#discard the 5' most window in case it is smaller than w
w <- as.numeric(binLength)
genesBed <- read.table(genesbed, sep = "\t")
colnames(genesBed) <- c("chr", "start", "end", "name", "score", "strand")

gc()

#genes binning
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

#genesBins object has row names allowing to match the gold standard peaks and the hitsMatrix bins
genesBins <- unlist(as(genesBinsList,"GRangesList"))
# Remove bins without length equal to w 
genesBins <- genesBins[which(width(genesBins) == w), ]
save(genesBins, file = "/path/to/genesBins.Rdata")

#find overlaps between reference-set peaks in RRACH motifs and genes bins
gs_RRACH_geneBins_hits <- findOverlaps(query = gs_RRACH, subject = genesBins)

#remove bins including more than one peak
ind_single_match <- which(isUnique(subjectHits(gs_RRACH_geneBins_hits)))
genesBins_filtered <- genesBins[subjectHits(gs_RRACH_geneBins_hits)[ind_single_match]]
gs_RRACH_filtered <- gs_RRACH[queryHits(gs_RRACH_geneBins_hits)[ind_single_match]]

#split genesBins by RRACH motif where the peak occurs
genesBins_filtered_split <- split(names(genesBins_filtered), mcols(gs_RRACH_filtered)$motif)

#initialize performances df
#evaluate performances for each tool and for each RRACH motif variant
performances <- data.frame()
for (i in 1:length(genesBins_filtered_split)) {
  ind_curr <- which(names(genesBins) %in% genesBins_filtered_split[[i]])
  motif_curr <- names(genesBins_filtered_split)[i]
  cat(sprintf("RRACH motif: %s\nNum. bins: %d\n", motif_curr, length(ind_curr)))
  hitsMatrix_subset <- hitsMatrix[ind_curr, ]
  for (j in 2:ncol(hitsMatrix_subset)) {
    tool_curr <- colnames(hitsMatrix_subset)[j]
    TP <- length(intersect(which(hitsMatrix_subset[, 1] == 1), which(hitsMatrix_subset[, j] == 1)))
    FP <- length(intersect(which(hitsMatrix_subset[, 1] == 0), which(hitsMatrix_subset[, j] == 1)))
    FN <- length(intersect(which(hitsMatrix_subset[, 1] == 1), which(hitsMatrix_subset[, j] == 0)))
    recall_curr <- TP/(TP + FN)
    #precision_curr <- TP/(TP + FP)
    #F1_score_curr <- 2*precision_curr*recall_curr/(precision_curr + recall_curr)
    cat(sprintf("Tool: %s\nRecall = %.2f\n\n", tool_curr, recall_curr))
    performances_curr <- data.frame(tool = tool_curr, motif = paste(motif_curr, " - ", length(ind_curr), " bins"), num_bins = length(ind_curr), accuracy = recall_curr)
    performances <- rbind(performances, performances_curr)
  }
  cat(sprintf("\n\n"))
}

#save to file
save(performances, file = "/path/to/RRACH_accuracy.Rdata")
performances$tool[grep(pattern = "Nanom6A", x = performances$tool)] <- "Nanom6A"

#plot barplot with accuracy split by RRACH motif variant
ggplot(performances, aes(x = tool, y = accuracy, fill = motif)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle("Accuracy split by RRACH motif") + 
  labs(x="Tool", y = "Accuracy") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        panel.background = element_blank())
ggsave("/path/to/RRACH_accuracy.pdf")

###############################
#False Positive, True Positive, False Negative bins motif enrichment analysis

#add sequence to each GRange 
genesBins_wseq <- genesBins
mcols(genesBins_wseq) <- cbind(mcols(genesBins), sequence = unname(as.character(getSeq(chrs, genesBins))))

control_seq <- genesBins_wseq[sample(1:length(genesBins), 1000)]

FP_bins_all <- list()
FP_bins_all[[1]] <- control_seq
TP_bins_all <- list()
TP_bins_all[[1]] <- control_seq
FN_bins_all <- list()
FN_bins_all[[1]] <- control_seq
for (i in 2:ncol(hitsMatrix)) {
  tool_curr <- colnames(hitsMatrix)[i]
  FP_bins_ind <- setdiff(which(hitsMatrix[, i] == 1), which(hitsMatrix[, 1] == 1))
  genesBins_FP_curr <- genesBins_wseq[FP_bins_ind]
  FP_bins_all[[i]] <- genesBins_FP_curr
  TP_bins_ind <- intersect(which(hitsMatrix[, i] == 1), which(hitsMatrix[, 1] == 1))
  genesBins_TP_curr <- genesBins_wseq[TP_bins_ind]
  TP_bins_all[[i]] <- genesBins_TP_curr
  FN_bins_ind <- intersect(which(hitsMatrix[, i] != 1), which(hitsMatrix[, 1] == 1))
  genesBins_FN_curr <- genesBins_wseq[FN_bins_ind]
  FN_bins_all[[i]] <- genesBins_FN_curr
}

names(FP_bins_all) <- c("control_seq", colnames(hitsMatrix)[2:dim(hitsMatrix)[2]])
save(FP_bins_all, file = "/path/to/FP_bins.RData")
names(TP_bins_all) <- c("control_seq", colnames(hitsMatrix)[2:dim(hitsMatrix)[2]])
save(TP_bins_all, file = "/path/to/TP_bins.RData")
names(FN_bins_all) <- c("control_seq", colnames(hitsMatrix)[2:dim(hitsMatrix)[2]])
save(FN_bins_all, file = "/path/to/FN_bins.RData")

#run XSTREME from MEME-suite online

################
#plot performances at default parameters in terms of recall, precision and F1 score

#oligos
#list tools
tools_oligos <- c("DiffErr", "DRUMMER", "ELIGOS", "EpiNano-Error", "EpiNano-SVM", "m6Anet", "MINES", "Nanocompore", "Nanom6A", "Tombo", "xPore", "Yanocomp")
#bin tools in classes
tools_oligos_class <- vector(mode = "character", length = length(tools_oligos))
names(tools_oligos_class) <- tools_oligos
tools_oligos_class[c("Nanocompore", "Tombo", "xPore", "Yanocomp")] <- "Transcriptome_Comparison"
tools_oligos_class[c("ELIGOS", "EpiNano-Error", "DiffErr", "DRUMMER")] <- "Genome_Comparison"
tools_oligos_class[c("m6Anet", "DENA", "Nanom6A", "MINES")] <- "Transcriptome_ML"
tools_oligos_class[c("EpiNano-SVM")] <- "Genome_ML"
tools_oligos_class <- tools_oligos_class[intersect(names(tools_oligos_class), tools_oligos)]
#performances at default parameters
recall_default_oligos <- c(0.9779650, 0.3817239, 0.8334413, 0.0900842515, 0.0635126377, 0.0401814647, 0.02333117, 0.9060272, 0.08360337, 0.9773169, 0.8813999, 0.9455606)
precision_default_oligos <- c(0.7526185, 0.7739816, 0.8460526, 0.9788732, 1, 0.9253731, 0.9473684, 0.7862767, 1, 0.7532468, 0.7758129, 0.7925041)
F1_score_default_oligos <- 2*recall_default_oligos*precision_default_oligos/(recall_default_oligos + precision_default_oligos)
names(recall_default_oligos) <- tools_oligos
names(precision_default_oligos) <- tools_oligos
names(F1_score_default_oligos) <- tools_oligos
val_F1_sorted <- unlist(lapply(split(F1_score_default_oligos, tools_oligos_class), function(x) sort(x, decreasing = TRUE)))
idx_order <- match(val_F1_sorted, F1_score_default_oligos)
tools_oligos <- tools_oligos[idx_order]
tools_oligos_class <- tools_oligos_class[idx_order]
recall_default_oligos <- recall_default_oligos[idx_order]
precision_default_oligos <- precision_default_oligos[idx_order]
F1_score_default_oligos <- F1_score_default_oligos[idx_order]
performances_default_oligos <- data.frame(tool = rep(tools_oligos, 3), tool_class = rep(tools_oligos_class, 3), metric = rep(c("Recall", "Precision", "F1 score"), each = length(tools_oligos)), metric_value = c(recall_default_oligos, precision_default_oligos, F1_score_default_oligos))
performances_default_oligos$tool <- factor(performances_default_oligos$tool, levels = unique(performances_default_oligos$tool))

#plot barplot with results
col_tool_class <- tools_oligos_class
col_tool_class[which(tools_oligos_class == "Genome_Comparison")] <- "red"
col_tool_class[which(tools_oligos_class == "Genome_ML")] <- "purple"
col_tool_class[which(tools_oligos_class == "Transcriptome_Comparison")] <- "blue" 
col_tool_class[which(tools_oligos_class == "Transcriptome_ML")] <- "green"

ggplot(performances_default_oligos, aes(x = tool, y = metric_value, fill = metric)) +
  geom_bar(stat="identity", width = 0.8, position=position_dodge(), colour = "black") +
  labs(x="Tool", y = "Metric value") +
  scale_fill_manual(values = c("white", "grey", "black")) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, colour = col_tool_class),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        panel.background = element_blank())
ggsave("/path/to/Metrics_default_oligos.pdf", device = "pdf", height = 5, width = 10)

#yeast
#list tools
tools_yeast <- c("DiffErr", "ELIGOS", "m6Anet", "MINES", "Nanocompore", "Nanom6A", "Tombo", "xPore", "Yanocomp")
#bin tools in classes
tools_yeast_class <- vector(mode = "character", length = length(tools_yeast))
names(tools_yeast_class) <- tools_yeast
tools_yeast_class[c("Nanocompore", "Tombo", "xPore", "Yanocomp")] <- "Transcriptome_Comparison"
tools_yeast_class[c("ELIGOS", "EpiNano-Error", "DiffErr", "DRUMMER")] <- "Genome_Comparison"
tools_yeast_class[c("m6Anet", "DENA", "Nanom6A", "MINES")] <- "Transcriptome_ML"
tools_yeast_class[c("EpiNano-SVM")] <- "Genome_ML"
tools_yeast_class <- tools_yeast_class[intersect(names(tools_yeast_class), tools_yeast)]
#performances at default parameters
recall_default_yeast <- c(0.0853406745, 0.086717137, 0.015141087, 0.1403992, 0.0949759119, 0.6834136, 0.878183070, 0.4969029594, 0.134205093) #values stored in associated Performances_window_50bp.tsv file; last line for each tool.
precision_default_yeast <- c(0.302439024, 0.41176471, 0.12865497, 0.05633803, 0.09458533, 0.01943363, 0.01116292, 0.01663058, 0.07959184) #values stored in associated Performances_window_50bp.tsv file; last line for each tool.
F1_score_default_yeast <- 2*recall_default_yeast*precision_default_yeast/(recall_default_yeast + precision_default_yeast)
names(recall_default_yeast) <- tools_yeast
names(precision_default_yeast) <- tools_yeast
names(F1_score_default_yeast) <- tools_yeast
val_F1_sorted <- unlist(lapply(split(F1_score_default_yeast, tools_yeast_class), function(x) sort(x, decreasing = TRUE)))
idx_order <- match(val_F1_sorted, F1_score_default_yeast)
tools_yeast <- tools_yeast[idx_order]
tools_yeast_class <- tools_yeast_class[idx_order]
recall_default_yeast <- recall_default_yeast[idx_order]
precision_default_yeast <- precision_default_yeast[idx_order]
F1_score_default_yeast <- F1_score_default_yeast[idx_order]
performances_default_yeast <- data.frame(tool = rep(tools_yeast, 3), tool_class = rep(tools_yeast_class, 3), metric = rep(c("Recall", "Precision", "F1 score"), each = length(tools_yeast)), metric_value = c(recall_default_yeast, precision_default_yeast, F1_score_default_yeast))
performances_default_yeast$tool <- factor(performances_default_yeast$tool, levels = unique(performances_default_yeast$tool))
#plot barplots with results
col_tool_class <- tools_yeast_class
col_tool_class[which(tools_yeast_class == "Genome_Comparison")] <- "red"
col_tool_class[which(tools_yeast_class == "Genome_ML")] <- "purple"
col_tool_class[which(tools_yeast_class == "Transcriptome_Comparison")] <- "blue" 
col_tool_class[which(tools_yeast_class == "Transcriptome_ML")] <- "green"

ggplot(performances_default_yeast, aes(x = tool, y = metric_value, fill = metric)) +
  geom_bar(stat="identity", width = 0.8, position=position_dodge(), colour = "black") +
  labs(x="Tool", y = "Metric value") +
  scale_fill_manual(values = c("white", "grey", "black")) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, colour = col_tool_class),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        panel.background = element_blank())
ggsave("/path/to/Metrics_default_yeast.pdf", height = 5, width = 10)

#mESC
#list tools
tools_mESC <- c("DENA", "DiffErr", "ELIGOS", "EpiNano-Error", "m6Anet", "MINES", "Nanocompore", "Nanom6A", "Tombo", "xPore", "Yanocomp")
#bin tools in classes
tools_mESC_class <- vector(mode = "character", length = length(tools_mESC))
names(tools_mESC_class) <- tools_mESC
tools_mESC_class[c("Nanocompore", "Tombo", "xPore", "Yanocomp")] <- "Transcriptome_Comparison"
tools_mESC_class[c("ELIGOS", "EpiNano-Error", "DiffErr", "DRUMMER")] <- "Genome_Comparison"
tools_mESC_class[c("m6Anet", "DENA", "Nanom6A", "MINES")] <- "Transcriptome_ML"
tools_mESC_class[c("EpiNano-SVM")] <- "Genome_ML"
tools_mESC_class <- tools_mESC_class[intersect(names(tools_mESC_class), tools_mESC)]
#performances at default parameters
#pick one reference set at a time

datasets_mESC <- c("mESC", "mESC_RRACH", "mESC_miCLIP2_AND_GLORI", "mESC_highcov", "mESC_miCLIP2_AND_GLORI_highcov", "mESC_miCLIP2_OR_GLORI_7DRACHfilt")

for (i in 1:length(datasets_mESC)) {
  if (datasets_mESC[i] == "mESC") {
    #GLORI_intersection OR miCLIP2
    recall_default_mESC <- c(0.32569593, 0.007244825, 0.01313347609, 0.07458957887, 0.06164644302, 0.1529146, 0.00963597430, 0.1500476, 0.0566024268, 0.145063050, 0.0758743754) #values stored in associated Performances_window_50bp.tsv file; last line for each tool.
    precision_default_mESC <- c(0.2503131, 0.513490725, 0.9192340, 0.2499601, 0.8217571, 0.4750009, 0.2682119, 0.3702519, 0.1546713, 0.1721733, 0.219401445) #values stored in associated Performances_window_50bp.tsv file; last line for each tool.
    #recall_default_mESC <- c(0.29052940, 0.005312328, 0.0102582891, 0.0527569152, 0.0452463821, 0.1168712, 0.0075105331, 0.1137571, 0.0426818099, 0.1104597912, 0.0518409965) #chr1 only
    #precision_default_mESC <- c(0.2616298, 0.5370370370, 0.9333333, 0.2167043, 0.8288591, 0.487395, 0.3129771, 0.3642229, 0.1628232, 0.1657504, 0.2012802276) #chr1 only
  } else if (datasets_mESC[i] == "mESC_RRACH") {
    #mESC GLORI_intersection OR miCLIP2 RRACH
    recall_default_mESC <- c(0.37549521, 0.007894586, 0.01419590056, 0.07669231211, 0.06756330022, 0.1754033, 0.01040649940, 0.1727048, 0.0569558477, 0.1462077281, 0.0777688465) #values stored in associated Performances_window_50bp.tsv file; last line for each tool.
    precision_default_mESC <- c(0.2650724, 0.674019608, 0.9518768, 0.4131931, 0.8617722, 0.4929407, 0.4205336, 0.3800019, 0.2682712, 0.2940786, 0.370841889) #values stored in associated Performances_window_50bp.tsv file; last line for each tool.
    recall_default_mESC_RRACH <- recall_default_mESC
    precision_default_mESC_RRACH <- precision_default_mESC
  } else if (datasets_mESC[i] == "mESC_miCLIP2_AND_GLORI") {
    #mESC GLORI_intersection AND miCLIP2
    recall_default_mESC <- c(0.46591890, 0.02015844, 0.03702251157, 0.14307004471, 0.14738410856, 0.2645698, 0.02494313279, 0.2697466, 0.0963997176, 0.230998510, 0.1421287944) #values stored in associated Performances_window_50bp.tsv file; last line for each tool.
    precision_default_mESC <- c(0.05430857, 0.2166947723, 0.3930058, 0.07271568, 0.29797019, 0.1246443, 0.10529801, 0.1009511, 0.03995189, 0.04158195, 0.0623323013) #values stored in associated Performances_window_50bp.tsv file; last line for each tool.
    recall_default_mESC_miCLIP2_AND_GLORI <- recall_default_mESC
    precision_default_mESC_miCLIP2_AND_GLORI <- precision_default_mESC
  } else if (datasets_mESC[i] == "mESC_highcov") {
    #mESC GLORI_intersection OR miCLIP2 high-coverage
    recall_default_mESC <- c(0.54088265, 0.0646596, 0.0864408713, 0.2811038887, 0.1585129433, 0.2881743, 0.0857566427, 0.532786, 0.412817881, 0.570418520, 0.3610445889) #values stored in associated Performances_window_50bp.tsv file; last line for each tool.
    precision_default_mESC <- c(0.3417147, 0.5344015, 0.9369592, 0.3361516, 0.8574954, 0.5785256, 0.2618384, 0.3338335, 0.1579131, 0.1830156, 0.2445165) #values stored in associated Performances_window_50bp.tsv file; last line for each tool.
    recall_default_mESC_highcov <- recall_default_mESC
    precision_default_mESC_highcov <- precision_default_mESC
  } else if (datasets_mESC[i] == "mESC_miCLIP2_AND_GLORI_highcov") {
    #mESC GLORI_intersection AND miCLIP2 high-coverage
    recall_default_mESC <- c(0.64602477, 0.09988014, 0.1458250100, 0.3695565322, 0.2740711147, 0.3979225, 0.1182580903, 0.5972833, 0.408709549, 0.587295246, 0.419896125) #values stored in associated Performances_window_50bp.tsv file; last line for each tool.
    precision_default_mESC <- c(0.11649856, 0.23562677, 0.4511743, 0.12614210, 0.42319556, 0.228022, 0.10306407, 0.1068239, 0.04462572, 0.05378508, 0.08117084) #values stored in associated Performances_window_50bp.tsv file; last line for each tool.
    recall_default_mESC_miCLIP2_AND_GLORI_highcov <- recall_default_mESC
    precision_default_mESC_miCLIP2_AND_GLORI_highcov <- precision_default_mESC
  } else if (datasets_mESC[i] == "mESC_miCLIP2_OR_GLORI_7DRACHfilt") {
    #mESC GLORI_intersection OR miCLIP2 7DRACH filtered
    recall_default_mESC <- c(0.32510178, 0.005406127, 0.01208035774, 0.06901154642, 0.06707601949, 0.1664553, 0.00744176734, 0.1421611, 0.0469865848, 0.1290796236, 0.0689114330) #values stored in associated Performances_window_50bp.tsv file; last line for each tool.
    precision_default_mESC <- c(0.4142359, 0.736363636, 0.9501312, 0.4723618, 0.8973214, 0.6133792, 0.5674300, 0.6127733, 0.3301290, 0.3621723, 0.448913043) #values stored in associated Performances_window_50bp.tsv file; last line for each tool.
    recall_default_mESC_7DRACH_filt <- recall_default_mESC
    precision_default_mESC_7DRACH_filt <- precision_default_mESC
  }
  
  F1_score_default_mESC <- 2*recall_default_mESC*precision_default_mESC/(recall_default_mESC + precision_default_mESC)
  names(recall_default_mESC) <- tools_mESC
  names(precision_default_mESC) <- tools_mESC
  names(F1_score_default_mESC) <- tools_mESC
  val_F1_sorted <- unlist(lapply(split(F1_score_default_mESC, tools_mESC_class), function(x) sort(x, decreasing = TRUE)))
  idx_order <- match(val_F1_sorted, F1_score_default_mESC)
  tools_mESC <- tools_mESC[idx_order]
  tools_mESC_class <- tools_mESC_class[idx_order]
  recall_default_mESC <- recall_default_mESC[idx_order]
  precision_default_mESC <- precision_default_mESC[idx_order]
  F1_score_default_mESC <- F1_score_default_mESC[idx_order]
  recall_default_mESC_base <- recall_default_mESC
  precision_default_mESC_base <- precision_default_mESC
  performances_default_mESC <- data.frame(tool = rep(tools_mESC, 3), tool_class = rep(tools_mESC_class, 3), metric = rep(c("Recall", "Precision", "F1 score"), each = length(tools_mESC)), metric_value = c(recall_default_mESC, precision_default_mESC, F1_score_default_mESC))
  performances_default_mESC$tool <- factor(performances_default_mESC$tool, levels = unique(performances_default_mESC$tool))
  #plot barplots with results
  col_tool_class <- tools_mESC_class
  col_tool_class[which(tools_mESC_class == "Genome_Comparison")] <- "red"
  col_tool_class[which(tools_mESC_class == "Genome_ML")] <- "purple"
  col_tool_class[which(tools_mESC_class == "Transcriptome_Comparison")] <- "blue" 
  col_tool_class[which(tools_mESC_class == "Transcriptome_ML")] <- "green"
  
  ggplot(performances_default_mESC, aes(x = tool, y = metric_value, fill = metric)) +
    geom_bar(stat="identity", width = 0.8, position=position_dodge(), colour = "black") +
    labs(x="Tool", y = "Metric value") +
    scale_fill_manual(values = c("white", "grey", "black")) +
    theme(legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, colour = col_tool_class),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          panel.background = element_blank())
  ggsave("/path/to/Metrics_default_mESC.pdf", height = 5, width = 10)
}

#HEK293T
#list tools
tools_HEK293T <- c("DENA", "DiffErr", "ELIGOS", "EpiNano-Error", "EpiNano-SVM", "m6Anet", "MINES", "Nanocompore", "Nanom6A", "Tombo", "xPore", "Yanocomp")
#bin tools in classes
tools_HEK293T_class <- vector(mode = "character", length = length(tools_HEK293T))
names(tools_HEK293T_class) <- tools_HEK293T
tools_HEK293T_class[c("Nanocompore", "Tombo", "xPore", "Yanocomp")] <- "Transcriptome_Comparison"
tools_HEK293T_class[c("ELIGOS", "EpiNano-Error", "DiffErr", "DRUMMER")] <- "Genome_Comparison"
tools_HEK293T_class[c("m6Anet", "DENA", "Nanom6A", "MINES")] <- "Transcriptome_ML"
tools_HEK293T_class[c("EpiNano-SVM")] <- "Genome_ML"
tools_HEK293T_class <- tools_HEK293T_class[intersect(names(tools_HEK293T_class), tools_HEK293T)]

datasets_HEK293T <- c("HEK293T", "HEK293T_RRACH", "HEK293T_highcov")

#performances at default parameters
for (i in 1:length(datasets_HEK293T)) {
  if (datasets_HEK293T[i] == "HEK293T") {
    recall_default_HEK293T <- c(0.40159518, 0.5302662159, 0.03776674414, NA, 0.185400380, 0.09534122121, 0.2430805, 0.03000211282, 0.2829072, 0.1777413902, 0.362243820, 0.0878406930) #values stored in associated Performances_window_50bp.tsv file; last line for each tool.
    precision_default_HEK293T <- c(0.2670812, 0.136455077, 0.9420290, NA, 0.11521796, 0.8861070, 0.4153055, 0.9265905, 0.2465022, 0.1380796, 0.1550987, 0.554518173) #values stored in associated Performances_window_50bp.tsv file; last line for each tool.
    } else if (datasets_HEK293T[i] == "HEK293T_RRACH") {
    recall_default_HEK293T <- c(0.49598163, 0.5383264672, 0.04410076315, NA, 0.216924428, 0.11237928007, 0.2996556, 0.03619909502, 0.3505774, 0.1839670426, 0.375970825, 0.0978591207) #values stored in associated Performances_window_50bp.tsv file; last line for each tool.
    precision_default_HEK293T <- c(0.2788790, 0.23485563, 0.9588840, NA, 0.16386919, 0.9043478, 0.4308604, 0.9520426, 0.2544608, 0.2377586, 0.2637889, 0.72161355) #values stored in associated Performances_window_50bp.tsv file; last line for each tool.
    recall_default_HEK293T_RRACH <- recall_default_HEK293T
    precision_default_HEK293T_RRACH <- precision_default_HEK293T
  } else if (datasets_HEK293T[i] == "HEK293T_highcov") {
    recall_default_HEK293T <- c(0.50836782, 0.8679143423, 0.1157099154, NA, 0.1707755983, 0.1795933057, 0.3491092, 0.0969947814, 0.6206586, 0.5790894367, 0.685801692, 0.2366384740) #values stored in associated Performances_window_50bp.tsv file; last line for each tool.
    precision_default_HEK293T <- c(0.2692271, 0.13601241, 0.9525926, NA, 0.18607843, 0.8966757, 0.4170249, 0.9373913, 0.235684, 0.1469675, 0.1664628, 0.64115066) #values stored in associated Performances_window_50bp.tsv file; last line for each tool.
    recall_default_HEK293T_highcov <- recall_default_HEK293T
    precision_default_HEK293T_highcov <- precision_default_HEK293T
  }
  
  F1_score_default_HEK293T <- 2*recall_default_HEK293T*precision_default_HEK293T/(recall_default_HEK293T + precision_default_HEK293T)
  names(recall_default_HEK293T) <- tools_HEK293T
  names(precision_default_HEK293T) <- tools_HEK293T
  names(F1_score_default_HEK293T) <- tools_HEK293T
  val_F1_sorted <- unlist(lapply(split(F1_score_default_HEK293T, tools_HEK293T_class), function(x) sort(x, decreasing = TRUE)))
  idx_order <- match(val_F1_sorted, F1_score_default_HEK293T)
  tools_HEK293T <- tools_HEK293T[idx_order]
  tools_HEK293T_class <- tools_HEK293T_class[idx_order]
  recall_default_HEK293T <- recall_default_HEK293T[idx_order]
  precision_default_HEK293T <- precision_default_HEK293T[idx_order]
  F1_score_default_HEK293T <- F1_score_default_HEK293T[idx_order]
  performances_default_HEK293T <- data.frame(tool = rep(tools_HEK293T, 3), tool_class = rep(tools_HEK293T_class, 3), metric = rep(c("Recall", "Precision", "F1 score"), each = length(tools_HEK293T)), metric_value = c(recall_default_HEK293T, precision_default_HEK293T, F1_score_default_HEK293T))
  performances_default_HEK293T$tool <- factor(performances_default_HEK293T$tool, levels = unique(performances_default_HEK293T$tool))
  #plot barplots with results
  col_tool_class <- tools_HEK293T_class
  col_tool_class[which(tools_HEK293T_class == "Genome_Comparison")] <- "red"
  col_tool_class[which(tools_HEK293T_class == "Genome_ML")] <- "purple"
  col_tool_class[which(tools_HEK293T_class == "Transcriptome_Comparison")] <- "blue" 
  col_tool_class[which(tools_HEK293T_class == "Transcriptome_ML")] <- "green"
  
  ggplot(performances_default_HEK293T, aes(x = tool, y = metric_value, fill = metric)) +
    geom_bar(stat="identity", width = 0.8, position=position_dodge(), colour = "black") +
    labs(x="Tool", y = "Metric value") +
    scale_fill_manual(values = c("white", "grey", "black")) +
    theme(legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, colour = col_tool_class),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          panel.background = element_blank())
  ggsave("/path/to/Metrics_default_HEK293T.pdf", height = 5, width = 10)
}

##############
#plot precision-recall curves for all the datasets

datasets <- c("oligos", "yeast", "mESC", "mESC_RRACH", "mESC_miCLIP2_AND_GLORI", "mESC_highcov", "mESC_miCLIP2_AND_GLORI_highcov", "mESC_miCLIP2_OR_GLORI_7DRACHfilt")
#cycle across datasets
for (d in 1:length(datasets)) {
  dataset <- datasets[d]
  #downsample points in the PR curve if there are too many
  num_points_plots <- c(2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000)
  names(num_points_plots) <- datasets
  num_points_plot <- num_points_plots[dataset]
  #list PR curves files
  if (dataset == "oligos") {
    PR_curves_files <- list.files(path = "/path/to/oligos/output_statistical_5/", pattern = "PRcurve_window_5bp\\.Rdata", full.names = TRUE)
  } else if (dataset == "yeast") {
    PR_curves_files <- list.files(path = "/path/to/yeast/output_statistical_50/", pattern = "PRcurve_window_50bp\\.Rdata", full.names = TRUE)
  } else if (dataset == "mESC") {
    PR_curves_files <- list.files(path = "/path/to/mESC/output_statistical_50_mESC_miCLIP2_OR_GLORI_intersection/", pattern = "PRcurve_window_50bp\\.Rdata", full.names = TRUE)
  } else if (dataset == "mESC_RRACH") {
    PR_curves_files <- list.files(path = "/path/to/mESC/output_statistical_50_mESC_miCLIP2_OR_GLORI_intersection/", pattern = "PRcurve_RRACH_window_50bp\\.Rdata", full.names = TRUE)
  } else if (dataset == "mESC_miCLIP2_AND_GLORI") {
    PR_curves_files <- list.files(path = "/path/to/mESC/output_statistical_50_mESC_miCLIP2_AND_GLORI_intersection/", pattern = "PRcurve_window_50bp\\.Rdata", full.names = TRUE)
  } else if (dataset == "mESC_miCLIP2_AND_GLORI_highcov") {
    PR_curves_files <- list.files(path = "/path/to/mESC/output_statistical_50_mESC_miCLIP2_AND_GLORI_intersection/", pattern = "PRcurve_highcov_window_50bp\\.Rdata", full.names = TRUE)
  } else if (dataset == "mESC_highcov") {
    PR_curves_files <- list.files(path = "/path/to/mESC/output_statistical_50_mESC_miCLIP2_OR_GLORI_intersection/", pattern = "PRcurve_highcov_window_50bp\\.Rdata", full.names = TRUE)
  } else if (dataset == "mESC_miCLIP2_OR_GLORI_7DRACHfilt") {
    PR_curves_files <- list.files(path = "/path/to/mESC/output_statistical_50_mESC_miCLIP2_OR_GLORI_intersection_7DRACH_filtered/", pattern = "7DRACH_filt_window_50bp\\.Rdata", full.names = TRUE)
  } else if (dataset == "HEK293T") {
    PR_curves_files <- list.files(path = "/Users/simonemaestri/Projects/nanopore_benchmark/HEK293T/output_statistical_50/", pattern = "PRcurve_window_50bp\\.Rdata", full.names = TRUE)
  } else if (dataset == "HEK293T_RRACH") {
    PR_curves_files <- list.files(path = "/Users/simonemaestri/Projects/nanopore_benchmark/HEK293T/output_statistical_50/", pattern = "PRcurve_RRACH_window_50bp\\.Rdata", full.names = TRUE)
  } else if (dataset == "HEK293T_highcov") {
    PR_curves_files <- list.files(path = "/Users/simonemaestri/Projects/nanopore_benchmark/HEK293T/output_statistical_50/", pattern = "PRcurve_highcov_window_50bp\\.Rdata", full.names = TRUE)
  } else {
    stop("Wrong dataset specified!\n")
  }
  #load PR curves
  PR_curves <- list()
  for (i in 1:length(PR_curves_files)) {
    load(PR_curves_files[i])
    PR_curves[[i]] <- pr
  }
  names(PR_curves) <- gsub(x = basename(PR_curves_files), pattern = "_PRcurve.*", replacement = "")
  #set graphical parameters
  lty_PR <- vector(mode = "numeric", length = length(PR_curves))
  names(lty_PR) <- names(PR_curves)
  tools_tx_comp <- c("Nanocompore", "Tombo", "xPore", "Yanocomp")
  tools_gen_comp <- c("ELIGOS", "EpiNano-Error", "DiffErr", "DRUMMER")
  tools_tx_ml <- c("m6Anet", "DENA", "Nanom6A")
  tools_gen_ml <- c("EpiNano-SVM")
  lty_PR[tools_tx_comp] <- 2 #transcriptome and multi-sample
  lty_PR[tools_gen_comp] <- 1 #genome and multi-sample
  lty_PR[tools_tx_ml] <- 3 #transcriptome and ML-based
  lty_PR[tools_gen_ml] <- 4 #genome and ML-based
  lwd_PR <- vector(mode = "numeric", length = length(PR_curves))
  names(lwd_PR) <- names(PR_curves)
  lwd_PR[tools_tx_comp] <- 2 #transcriptome and multi-sample
  lwd_PR[tools_gen_comp] <- 2 #genome and multi-sample
  lwd_PR[tools_tx_ml] <- 4 #transcriptome and ML-based
  lwd_PR[tools_gen_ml] <- 2 #genome and ML-based
  
  pch_PR <- vector(mode = "numeric", length = length(PR_curves))
  names(pch_PR) <- names(PR_curves)
  
  pch_PR[tools_tx_comp] <- NA #transcriptome and multi-sample
  pch_PR[tools_gen_comp] <- NA #genome and multi-sample
  pch_PR[tools_tx_ml] <- NA #transcriptome and ML-based
  pch_PR[tools_gen_ml] <- NA #genome and ML-based
  
  colors_PR <- c(brewer.pal(n = 12, name = 'Paired'))
  colors_PR[11] <- "gold"
  colors_PR <- colors_PR[1:(length(tools_gen_comp) + length(tools_tx_comp) + length(tools_tx_ml) + length(tools_gen_ml))]
  names(colors_PR) <- c(tools_tx_comp, tools_gen_comp, tools_tx_ml, tools_gen_ml)
  
  legend_position <- c("topright", "topright", "topright", "topright", "topright", "topright", "bottomright", "topright")
  names(legend_position) <- c("mESC", "mESC_RRACH", "mESC_miCLIP2_AND_GLORI", "mESC_miCLIP2_AND_GLORI_highcov", "mESC_highcov", "mESC_miCLIP2_OR_GLORI_7DRACHfilt", "oligos", "yeast")
  
  #open pdf file and load default recall and precision values for each tool at default conditions
  if (dataset == "oligos") {
    pdf("/path/to/oligos_PR_curves_summary.pdf", width = 6, height = 6)
    recall_default_curr <- recall_default_oligos
    precision_default_curr <- precision_default_oligos
  } else if (dataset == "yeast") {
    pdf("/path/to/yeast_PR_curves_summary.pdf", width = 6, height = 6)
    recall_default_curr <- recall_default_yeast
    precision_default_curr <- precision_default_yeast
  } else if (dataset == "mESC_RRACH") {
    pdf("/path/to/mESC_RRACH_PR_curves_summary.pdf", width = 6, height = 6)
    recall_default_curr <- recall_default_mESC_RRACH
    precision_default_curr <- precision_default_mESC_RRACH
  } else if (dataset == "mESC_miCLIP2_AND_GLORI") {
    pdf("/path/to/mESC_miCLIP2_AND_GLORI_PR_curves_summary.pdf", width = 6, height = 6)
    recall_default_curr <- recall_default_mESC_miCLIP2_AND_GLORI
    precision_default_curr <- precision_default_mESC_miCLIP2_AND_GLORI
  } else if (dataset == "mESC_miCLIP2_AND_GLORI_highcov") {
    pdf("/path/to/mESC_miCLIP2_AND_GLORI_highcov_PR_curves_summary.pdf", width = 6, height = 6)
    recall_default_curr <- recall_default_mESC_miCLIP2_AND_GLORI_highcov
    precision_default_curr <- precision_default_mESC_miCLIP2_AND_GLORI_highcov
  } else if (dataset == "mESC_highcov") {
    pdf("/path/to/mESC_highcov_PR_curves_summary.pdf", width = 6, height = 6)
    recall_default_curr <- recall_default_mESC_highcov
    precision_default_curr <- precision_default_mESC_highcov
  } else if (dataset == "mESC") {
    pdf("/path/to/mESC_PR_curves_summary.pdf", width = 6, height = 6)
    recall_default_curr <- recall_default_mESC_base
    precision_default_curr <- precision_default_mESC_base
  } else if (dataset == "mESC_miCLIP2_OR_GLORI_7DRACHfilt") {
    pdf("/path/to/mESC_PR_curves_summary_7DRACH_filtered.pdf", width = 6, height = 6)
    recall_default_curr <- recall_default_mESC_7DRACH_filtered
    precision_default_curr <- precision_default_mESC_7DRACH_filtered
  } else if (dataset == "HEK293T") {
    pdf("/path/to/HEK293T_PR_curves_summary.pdf", width = 6, height = 6)
    recall_default_curr <- recall_default_HEK293T_base
    precision_default_curr <- precision_default_HEK293T_base
  } else if (dataset == "HEK293T_RRACH") {
    pdf("/path/to/HEK293T_RRACH_PR_curves_summary.pdf", width = 6, height = 6)
    recall_default_curr <- recall_default_HEK293T_RRACH
    precision_default_curr <- precision_default_HEK293T_RRACH
  } else if (dataset == "HEK293T_highcov") {
    pdf("/path/to/HEK293T_highcov_PR_curves_summary.pdf", width = 6, height = 6)
    recall_default_curr <- recall_default_HEK293T_highcov
    precision_default_curr <- precision_default_HEK293T_highcov
  } else {
    stop("Wrong dataset specified!\n")
  }
  #cycle across tools and plot PR curves in the same order among datasets
  for (i in match(intersect(c(tools_gen_comp, tools_tx_comp, tools_tx_ml, tools_gen_ml), names(PR_curves)), names(PR_curves))) {
    if (length(PR_curves[[i]]$curve[, 1]) > num_points_plot) {
      range_x <- range(PR_curves[[i]]$curve[, 1])
      quantiles_x <- seq(from = range_x[1], to = range_x[2], length.out = 101)
      for (k in 1:(length(quantiles_x) - 1)) {
        indexes_range_curr <- intersect(which(PR_curves[[i]]$curve[, 1] >= quantiles_x[k]), which(PR_curves[[i]]$curve[, 1] <= quantiles_x[k + 1]))
        subset_point <- c(sort(resample(x = indexes_range_curr, size = min(length(indexes_range_curr), round(num_points_plots[i]/100)))), subset_point)
      }
      #subset_point <- sort(sample(1:length(PR_curves[[i]]$curve[, 1]), num_points_plot))
    } else {
      subset_point <- 1:length(PR_curves[[i]]$curve[, 1])
    }
    if (!is.na(pch_PR[names(PR_curves)[[i]]])) {
      subset_point <- sort(sample(1:length(PR_curves[[i]]$curve[, 1]), 100))
    }
    if (i == match(intersect(c(tools_gen_comp, tools_tx_comp, tools_tx_ml, tools_gen_ml), names(PR_curves)), names(PR_curves))[1]) {
      plot(seq(from = min(PR_curves[[i]]$rand$curve[, 1], na.rm = TRUE), to = max(PR_curves[[i]]$rand$curve[, 1], na.rm = TRUE), length.out = 100), seq(from = min(PR_curves[[i]]$rand$curve[, 2], na.rm = TRUE), to = max(PR_curves[[i]]$rand$curve[, 2], na.rm = TRUE), length.out = 100), type = 'l', pch = NA, lty = 2, lwd = 2, col = "lightgrey", main = "Summary PR curves", xlab = "Recall", ylab = "Precision", xlim = c(0, 1), ylim = c(0, 1))
      lines(PR_curves[[i]]$curve[subset_point, 1], PR_curves[[i]]$curve[subset_point, 2], lwd = lwd_PR[names(PR_curves)[[i]]], col = colors_PR[names(PR_curves)[[i]]])
      points(recall_default_curr[names(PR_curves)[i]], precision_default_curr[names(PR_curves)[i]], pch = 15, cex = 2, col = colors_PR[names(PR_curves)[[i]]])
    } else {
      lines(PR_curves[[i]]$curve[subset_point, 1], PR_curves[[i]]$curve[subset_point, 2], lty = lty_PR[names(PR_curves)[i]], lwd = lwd_PR[names(PR_curves)[[i]]], col = colors_PR[names(PR_curves)[[i]]])
      points(recall_default_curr[names(PR_curves)[i]], precision_default_curr[names(PR_curves)[i]], pch = 15, cex = 2, col = colors_PR[names(PR_curves)[[i]]])
    }
  }
  legend(legend_position[dataset], legend = c("Random", names(PR_curves)[match(intersect(c(tools_gen_comp, tools_tx_comp, tools_tx_ml, tools_gen_ml), names(PR_curves)), names(PR_curves))]), col = c("lightgrey", colors_PR[names(PR_curves)[match(intersect(c(tools_gen_comp, tools_tx_comp, tools_tx_ml, tools_gen_ml), names(PR_curves)), names(PR_curves))]]), pch = c(NA, pch_PR[match(intersect(c(tools_gen_comp, tools_tx_comp, tools_tx_ml, tools_gen_ml), names(PR_curves)), names(PR_curves))]), lty = c(2, lty_PR[match(intersect(c(tools_gen_comp, tools_tx_comp, tools_tx_ml, tools_gen_ml), names(PR_curves)), names(PR_curves))]), cex = 1, lwd = c(2, lwd_PR[match(intersect(c(tools_gen_comp, tools_tx_comp, tools_tx_ml, tools_gen_ml), names(PR_curves)), names(PR_curves))]))
  dev.off()
}

############
#plot running time

time_h_oligos <- c(310.8, 429.8, 109.8, 478.8, 295.8, 295.8, 295.8, 295.8, 109.8, 109.8, 109.8, 109.8, 30, 90, 1.2, 0.2, 4.1, 150.5, 23.5, 16, 42, 6, 6, 6) # derived from NextflowTower report.
tools_oligos <- c("Tombo", "Nanom6A", "DiffErr", "MINES", "xPore", "Nanocompore", "Yanocomp", "m6Anet", "DRUMMER", "EpiNano-SVM", "EpiNano-Error", "ELIGOS")
tools_oligos <- factor(tools_oligos, levels = names(sort(unlist(lapply(split(time_h_oligos, tools_oligos), sum)), decreasing = TRUE)))
analysis_type_oligos <- factor(rep(c("Preprocessing", "Tool execution"), each = length(tools_oligos)), levels = c("Tool execution", "Preprocessing"))
oligos_running_time <- data.frame(tool = rep(tools_oligos, 2), "analysis_type" = analysis_type_oligos, time_h = time_h_oligos)

ggplot(oligos_running_time, aes(x = tool, y = time_h, fill = analysis_type, group = analysis_type)) +
  geom_bar(stat="identity", width = 0.8, position="stack") + 
  labs(x="Tool", y = "Running time\n(CPU hours)") +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        panel.background = element_blank())
ggsave("/path/to/Running_time_oligos.pdf", height = 5, width = 7)  

time_h_yeast <- c(5814.6, 5814.6, 74.1, 8109.6, 4146.6, 4146.6, 4146.6, 4146.6, 444.6, 2295, 252, 24, 0.2, 120, 1582, 435, 470, 102) # derived from NextflowTower report.
tools_yeast <- c("Tombo", "Nanom6A", "DiffErr", "MINES", "xPore", "Nanocompore", "Yanocomp", "m6Anet", "ELIGOS")
tools_yeast <- factor(tools_yeast, levels = names(sort(unlist(lapply(split(time_h_yeast, tools_yeast), sum)), decreasing = TRUE)))
analysis_type_yeast <- factor(rep(c("Preprocessing", "Tool execution"), each = length(tools_yeast)), levels = c("Tool execution", "Preprocessing"))
yeast_running_time <- data.frame(tool = rep(tools_yeast, 2), analysis_type = analysis_type_yeast, time_h = time_h_yeast)

ggplot(yeast_running_time, aes(x = tool, y = time_h, fill = analysis_type, group = analysis_type)) +
  geom_bar(stat="identity", width = 0.8, position="stack") + 
  labs(x="Tool", y = "Running time\n(CPU hours)") +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        panel.background = element_blank())
ggsave("/path/to/Running_time_yeast.pdf", height = 5, width = 7)  

time_h_mESC <- c(4227.8, 4227.8, 667.8, 4395.8, 4227.8, 4587.8, 4587.8, 4587.8, 4587.8, 667.8, 240, 1180, 10, 0.2, 700, 95, 1560, 1220, 502, 5572) # derived from NextflowTower report.
tools_mESC <- c("Tombo", "Nanom6A", "DiffErr", "MINES", "DENA", "xPore", "Nanocompore", "Yanocomp", "m6Anet", "ELIGOS")
tools_mESC <- factor(tools_mESC, levels = names(sort(unlist(lapply(split(time_h_mESC, tools_mESC), sum)), decreasing = TRUE)))
analysis_type_mESC <- factor(rep(c("Preprocessing", "Tool execution"), each = length(tools_mESC)), levels = c("Tool execution", "Preprocessing"))
mESC_running_time <- data.frame(tool = rep(tools_mESC, 2), analysis_type = analysis_type_mESC, time_h = time_h_mESC)

ggplot(mESC_running_time, aes(x = tool, y = time_h, fill = analysis_type, group = analysis_type)) +
  geom_bar(stat="identity", width = 0.8, position="stack") + 
  labs(x="Tool", y = "Running time\n(CPU hours)") +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        panel.background = element_blank())
ggsave("/path/to/Running_time_mESC.pdf", height = 5, width = 7)

time_h_HEK293T <- c(612.2, 612.2, 94.2, 992.2, 94.2, 94.2, 94.2, 94.2, 94.2, 94.2, 612.2, 20, 740, 20, 2, 140, 870, 21.5, 102, 672, 360, 280)
tools_HEK293T <- c("Tombo", "Nanom6A", "DiffErr", "MINES", "xPore", "Nanocompore", "Yanocomp", "m6Anet", "ELIGOS", "EpiNano-SVM", "DENA")
tools_HEK293T <- factor(tools_HEK293T, levels = names(sort(unlist(lapply(split(time_h_HEK293T, tools_HEK293T), sum)), decreasing = TRUE)))
analysis_type_HEK293T <- factor(rep(c("Preprocessing", "Tool execution"), each = length(tools_HEK293T)), levels = c("Tool execution", "Preprocessing"))
HEK293T_running_time <- data.frame(tool = rep(tools_HEK293T, 2), analysis_type = analysis_type_HEK293T, time_h = time_h_HEK293T)

ggplot(HEK293T_running_time, aes(x = tool, y = time_h, fill = analysis_type, group = analysis_type)) +
  geom_bar(stat="identity", width = 0.8, position="stack") + 
  labs(x="Tool", y = "Running time\n(CPU hours)") +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        panel.background = element_blank())
ggsave("/path/to/Running_time_HEK293T.pdf", height = 5, width = 7)  

############
#plot downsampling analyses

tools_HEK293T <- c("DENA", "DiffErr", "ELIGOS", "EpiNano-Error", "EpiNano-SVM", "m6Anet", "MINES", "Nanocompore", "Nanom6A", "Tombo", "xPore", "Yanocomp")
lty_PR <- vector(mode = "numeric", length = length(tools_HEK293T))
names(lty_PR) <- tools_HEK293T
tools_tx_comp <- c("Nanocompore", "Tombo", "xPore", "Yanocomp", "MINES") #MINES added
tools_gen_comp <- c("ELIGOS", "EpiNano-Error", "DiffErr", "DRUMMER")
tools_tx_ml <- c("m6Anet", "DENA", "Nanom6A")
tools_gen_ml <- c("EpiNano-SVM")

lty_PR[tools_tx_comp] <- 2 #transcriptome and multi-sample
lty_PR[tools_gen_comp] <- 1 #genome and multi-sample
lty_PR[tools_tx_ml] <- 3 #transcriptome and ML-based
lty_PR[tools_gen_ml] <- 4 #genome and ML-based
lwd_PR <- vector(mode = "numeric", length = length(tools_HEK293T))
names(lwd_PR) <- tools_HEK293T
lwd_PR[tools_tx_comp] <- 2 #transcriptome and multi-sample
lwd_PR[tools_gen_comp] <- 2 #genome and multi-sample
lwd_PR[tools_tx_ml] <- 4 #transcriptome and ML-based
lwd_PR[tools_gen_ml] <- 2 #genome and ML-based

colors_PR <- c(brewer.pal(n = 12, name = 'Paired'))
colors_PR[11] <- "gold"
colors_PR <- colors_PR[1:(length(tools_gen_comp) - 1 + length(tools_tx_comp) + length(tools_tx_ml) + length(tools_gen_ml))]
colors_PR <- c(colors_PR, "black")
names(colors_PR) <- c(tools_tx_comp[-length(tools_tx_comp)], tools_gen_comp, tools_tx_ml, tools_gen_ml, tools_tx_comp[length(tools_tx_comp)])

#25perc dataset
F1_score_default_HEK293T_25perc <- c(0.3238242, 0.312690676, 0.0183326643, 0.1395275591, 0.207060894, 0.0724111082, 0.2609279, 0.0100779360, 0.1699717, 0.0662697967, 0.201093245, NA) #values stored in associated Performances_window_50bp.tsv file; last line for each tool.
names(F1_score_default_HEK293T_25perc) <- tools_HEK293T
AUC_default_HEK293T_25perc <- c(0.1040578, 0.1144668, 0.05006911, 0.03699215, 0.06709548, 0.1008135, NA, 0.0142202, 0.04065803, 0.02208066, 0.06310136, NA) #values reported in PR_curves.pdf file
names(AUC_default_HEK293T_25perc) <- tools_HEK293T
Num_hits_default_HEK293T_25perc <- c(26029, 1980260, 213, 16822, 21279, 728, 8522, 96, 9071, 23921, 33700, NA) #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
names(Num_hits_default_HEK293T_25perc) <- tools_HEK293T
#50perc dataset
F1_score_default_HEK293T_50perc <- c(0.31356605, 0.2079635140, 0.0387043361, 0.1615737785, 0.156928908, 0.1114799584, 0.2696946, 0.0261471952, 0.207709, 0.1061134569, 0.189981171, 0.0857489682) #values stored in associated Performances_window_50bp.tsv file; last line for each tool.
names(F1_score_default_HEK293T_50perc) <- tools_HEK293T
AUC_default_HEK293T_50perc <- c(0.1230267, 0.1507256, 0.08666965, 0.06148881, 0.06315302, 0.164139, NA, 0.03220917, 0.06857422, 0.05152056, 0.1086919, 0.05695282) #values reported in PR_curves.pdf file
names(AUC_default_HEK293T_50perc) <- tools_HEK293T
Num_hits_default_HEK293T_50perc <- c(35108, 1429357, 555, 28242, 28047, 1463, 12599, 301, 16257, 55814, 51545, 2464) #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
names(Num_hits_default_HEK293T_50perc) <- tools_HEK293T
#75perc dataset
F1_score_default_HEK293T_75perc <- c(0.31694807, 0.2119244188, 0.0555840427, 0.1948570956, 0.151693952, 0.1463107408, 0.2935723, 0.0406243539, 0.2427354, 0.1379579720, 0.206357224, 0.1239657832) #values stored in associated Performances_window_50bp.tsv file; last line for each tool.
names(F1_score_default_HEK293T_75perc) <- tools_HEK293T
AUC_default_HEK293T_75perc <- c(0.1300546, 0.1755889, 0.1053225, 0.08356577, 0.06039177, 0.2076433, NA, 0.04724435, 0.08934331, 0.07465844, 0.1372065, 0.08102062) #values reported in PR_curves.pdf file
names(AUC_default_HEK293T_75perc) <- tools_HEK293T
Num_hits_default_HEK293T_75perc <- c(40731, 762111, 820, 36635, 32576, 2067, 15435, 503, 21677, 83574, 61704, 3756) #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
names(Num_hits_default_HEK293T_75perc) <- tools_HEK293T
#Full dataset
F1_score_default_HEK293T_full <- c(0.32080846, 0.2170547664, 0.0726220101, NA, 0.142116771, 0.1721589012, 0.3066671, 0.0581222819, 0.263453, 0.1554200730, 0.217200589, 0.1516574712) #values stored in associated Performances_window_50bp.tsv file; last line for each tool.
names(F1_score_default_HEK293T_full) <- tools_HEK293T
AUC_default_HEK293T_full <- c(0.1351879, 0.1921865, 0.1254483, NA, 0.05593083, 0.2410183, NA, 0.06591276, 0.1029251, 0.09603134, 0.1614879, 0.09951042) #values reported in PR_curves.pdf file
names(AUC_default_HEK293T_full) <- tools_HEK293T
Num_hits_default_HEK293T_full <- c(45281, 235104, 1114,  NA, 36403, 2627, 17772, 756, 25931, 109387, 72127, 4700) #number of lines in the associated output bed file after filtering at default parameters; previous num_filt_peaks object.
names(Num_hits_default_HEK293T_full) <- tools_HEK293T

#create dataframes
performances_default_HEK293T_F1score_downsampling <- data.frame(tool = rep(tools_HEK293T, 4), downsampling_value = rep(c("0.25", "0.50", "0.75", "1"), each = length(tools_HEK293T)), F1_score = c(F1_score_default_HEK293T_25perc, F1_score_default_HEK293T_50perc, F1_score_default_HEK293T_75perc, F1_score_default_HEK293T_full))
performances_default_HEK293T_AUC_downsampling <- data.frame(tool = rep(tools_HEK293T, 4), downsampling_value = rep(c("0.25", "0.50", "0.75", "1"), each = length(tools_HEK293T)), AUC = c(AUC_default_HEK293T_25perc, AUC_default_HEK293T_50perc, AUC_default_HEK293T_75perc, AUC_default_HEK293T_full))
performances_default_HEK293T_Num_hits_downsampling <- data.frame(tool = rep(tools_HEK293T, 4), downsampling_value = rep(c("0.25", "0.50", "0.75", "1"), each = length(tools_HEK293T)), Num_hits = c(Num_hits_default_HEK293T_25perc, Num_hits_default_HEK293T_50perc, Num_hits_default_HEK293T_75perc, Num_hits_default_HEK293T_full))

#define function to plot downsampling results
plot_downsampling = function(data, legend = FALSE, ylim_range_coeff = c(0.1, 0.1), title, tools, legend_flag = TRUE, ylabel = "", var_name = "", logyaxis_flag = "") {
  pdf(title)
  if (legend_flag == TRUE) {
    par(mar=c(5.1, 4.1, 4.1, 11.1), xpd = TRUE)
  }
  if (logyaxis_flag != "") {
    ylimval <- c(max(0.01, (range(data[, var_name], na.rm = TRUE)[1] - ylim_range_coeff[1]*(range(data[, var_name], na.rm = TRUE)[2] - range(data[, var_name], na.rm = TRUE)[1]))), (range(data[, var_name], na.rm = TRUE)[2] + ylim_range_coeff[2]*(range(data[, var_name], na.rm = TRUE)[2] - range(data[, var_name], na.rm = TRUE)[1])))
    logval <- "y"
  } else {
    ylimval <- c((range(data[, var_name], na.rm = TRUE)[1] - ylim_range_coeff[1]*(range(data[, var_name], na.rm = TRUE)[2] - range(data[, var_name], na.rm = TRUE)[1])), (range(data[, var_name], na.rm = TRUE)[2] + ylim_range_coeff[2]*(range(data[, var_name], na.rm = TRUE)[2] - range(data[, var_name], na.rm = TRUE)[1])))
    logval <- ""
  }
  for (i in seq_along(tools)) {
    ind_tool_curr <- which(data$tool == tools[i])
    if (i == 1) {
      plot(data[ind_tool_curr, "downsampling_value"], data[ind_tool_curr, var_name], lwd = lwd_PR[tools[i]], col = colors_PR[tools[i]], xlim = c(0.25, 1), ylim = ylimval, xlab = "Fraction of dataset", ylab = ylabel, type = 'l', lty = lty_PR[tools[i]], xaxt = "n", log = logval)
      axis(side = 1, at = c(0.25, 0.5, 0.75, 1), labels = c("25%", "50%", "75%", "100%"))
      if (logyaxis_flag != "") {
        options(scipen=0)
      }
      points(data[ind_tool_curr, "downsampling_value"], data[ind_tool_curr, var_name], lwd = lwd_PR[tools[i]], col = colors_PR[tools[i]], pch = 15, cex = 2)
    } else {
      points(data[ind_tool_curr, "downsampling_value"], data[ind_tool_curr, var_name], lwd = lwd_PR[tools[i]], col = colors_PR[tools[i]], lty = lty_PR[tools[i]], type = 'l')
      points(data[ind_tool_curr, "downsampling_value"], data[ind_tool_curr, var_name], lwd = lwd_PR[tools[i]], col = colors_PR[tools[i]], pch = 15, cex = 2)
    }
  }
  if (legend_flag == TRUE) {
    legend("topright", inset=c(-0.5,0), legend = tools, col = colors_PR[tools], lty = lty_PR[tools], cex = 1, lwd = lwd_PR[tools])
  }
  dev.off()
}

#do plots
plot_downsampling(data=performances_default_HEK293T_F1score_downsampling, ylim_range_coeff = c(0.1, 0.1), title="/path/to/HEK293T_downsampling_F1score.pdf", tools = tools_HEK293T, legend_flag = FALSE, ylabel = "F1 score", var_name = "F1_score")
plot_downsampling(data=performances_default_HEK293T_AUC_downsampling, ylim_range_coeff = c(0.1, 0.1), title="/path/to/HEK293T_downsampling_AUC.pdf", tools = tools_HEK293T, legend_flag = TRUE, ylabel = "AUC", var_name = "AUC")
plot_downsampling(data=performances_default_HEK293T_Num_hits_downsampling, ylim_range_coeff = c(0.1, 10), title="/path/to/HEK293T_downsampling_Num_hits.pdf", tools = tools_HEK293T, legend_flag = FALSE, ylabel = "Log10(Num. hits)", var_name = "Num_hits", logyaxis_flag = 1)

#split hits by tool
num_hits_split <- split(performances_default_HEK293T_Num_hits_downsampling, performances_default_HEK293T_Num_hits_downsampling$tool)
#normalise the number of hits within tool
num_hits_norm <- lapply(num_hits_split, function(x) {
  max_hits <- x$Num_hits[which(x$downsampling_value == 1)]
  hits_norm <- x$Num_hits/max_hits
  x$Num_hits_norm <- hits_norm
  return(x)
})

#create dataframe
performances_default_HEK293T_Num_hits_downsampling_norm <- data.frame(tool=unname(unlist(lapply(num_hits_norm, '[[', "tool"))), downsampling_value=as.numeric(unname(unlist(lapply(num_hits_norm, '[[', "downsampling_value")))), num_hits_norm=as.numeric(unname(unlist(lapply(num_hits_norm, '[[', "Num_hits_norm")))))
#do plot
plot_downsampling(data=performances_default_HEK293T_Num_hits_downsampling_norm, ylim_range_coeff = c(0.1, 0.1), title="/path/toHEK293T_downsampling_Num_hits_norm.pdf", tools = tools_HEK293T, legend_flag = FALSE, ylabel = "Num. hits compared to full dataset", var_name = "num_hits_norm", logyaxis_flag = "")
#create dataframe without DiffErr
performances_default_HEK293T_Num_hits_downsampling_norm_noDiffErr <- performances_default_HEK293T_Num_hits_downsampling_norm[-which(performances_default_HEK293T_Num_hits_downsampling_norm$tool == "DiffErr"), ]
#do plot
plot_downsampling(data=performances_default_HEK293T_Num_hits_downsampling_norm_noDiffErr, ylim_range_coeff = c(0.1, 0.1), title="/path/to/HEK293T_downsampling_Num_hits_norm_noDiffErr.pdf", tools = tools_HEK293T, legend_flag = FALSE, ylabel = "Num. hits compared to full dataset", var_name = "num_hits_norm", logyaxis_flag = "")
