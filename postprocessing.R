#!/usr/bin/env Rscript
### load input variables ###
args = commandArgs(trailingOnly=TRUE)

for(v in args)
{
  vTmp <- strsplit(v,"=")[[1]]
  assign(vTmp[[1]],vTmp[[2]])
}

library(IRanges)
library(ensembldb)
library(GenomicRanges)
library(stringr)
library(Biostrings)
library(parallel)

output_processing <- function(tool, path_folder, output_file, filtering_parameter, genome_gtf, genome_bed){
  if (!file.exists(path_folder)){
    message(paste0(tool,"'s output folder doesn't exist."))
  }
  else {
    message(paste0("Processing output file from tool: ", tool,"\n"))
    switch (tool,
            differr = {output_differr <- function(){if (file.exists(paste0(path_folder, "/", "differrOut.bed")) 
                                                        && file.info(paste0(path_folder, "/", "differrOut.bed"))$size != 0){
              data_differr <- read.table(paste0(path_folder, "/", "differrOut.bed"))
              differr <- data_differr[, c(1:3, 6)]
              differr$V2 <- differr$V2 + 1
              differr$V3 <- differr$V3
              differr <- cbind(differr, rep("Mod", nrow(differr)))
              differr <- cbind(differr, 10**(-data_differr$V5)) # Parameter of filtering is FDR not -log10(FDR)
              colnames(differr) <- c("Chr", "Start", "End", "Strand", "Status", "FDR")
              write.table(differr, file = output_file, quote = F, sep = "\t", row.names = F)
            }
              else {message(paste0(tool,"'s output files don't exist."))}
            }
            if (!file.exists(output_file)) {
              output_differr()
            }
            },
            dena = {output_dena <- function(){if (file.exists(paste0(path_folder, "/", "dena_label.tsv"))
                                                  && file.info(paste0(path_folder, "/", "dena_label.tsv"))$size != 0){
              data_dena <- read.table(paste0(path_folder, "/", "dena_label.tsv"))
              dena <- data_dena[, c(1,2,2,6,6)]
              dena$V2 <- dena$V2
              dena$V2.1 <- dena$V2.1
              dena$V6 <- ifelse(!is.nan(dena$V6) & !is.na(dena$V6) & dena$V6 > filtering_parameter, "Mod", "Unmod")
              dena <- dena[which(dena$V6 == "Mod"), ]
              # Creation Edb Database from genome GTF
              EnsDb <- suppressWarnings(suppressMessages(ensDbFromGtf(gtf = genome_gtf)))
              edb <- EnsDb(EnsDb)
              # Lift-over + output bed file
              test_dena <- IRanges(start = dena[,2], end = dena[,3], names = c(dena[,1]))
              
              num_rows_chunk <- 1
              mc.cores <- as.numeric(mccores)
              if (length(test_dena) < num_rows_chunk) {
                test_dena_split <- list(test_dena)
              } else {
                test_dena_split <- split(test_dena, rep(seq(from = 1, to = ceiling(length(test_dena)/num_rows_chunk)), each = num_rows_chunk)[1:length(test_dena)])
              }
              
              tmp1 <- vector(mode = "list", length = length(test_dena_split))
              names(tmp1) <- 1:length(test_dena_split)
              tmp <- vector(mode = "list", length = length(test_dena_split))
              names(tmp) <- 1:length(test_dena_split)
              ind_retry <- 1:length(test_dena_split)
              while(any(unlist(lapply(tmp, is.null)))) {
                cat(sprintf("Starting new iteration for dena; %d sites missing\n", length(which(unlist(lapply(tmp, is.null))))))
                tmp1 <- tmp1[ind_retry]
                tmp1 <- mclapply(test_dena_split[ind_retry], function(x) {
                  tryCatch({
                    coordinate_dena_unlisted <- unlist(transcriptToGenome(x, edb))
                    return(coordinate_dena_unlisted)
                  }, warning = function(w) {
                    print("Warning")
                    return(NULL)
                  }, error = function(e) {
                    print("Error")
                    return(NULL)
                  }
                  )}, mc.cores = mc.cores)
                ind_retry <- names(which(unlist(lapply(tmp1, function(x) is.null(x)))))
                ind_ok <- names(which(unlist(lapply(tmp1, function(x) !is.null(x)))))
                tmp[ind_ok] <- tmp1[ind_ok]
                if (length(ind_retry) > 0) {
                  tmp1 <- tmp1[ind_retry]
                }
              }
              
              coordinate_dena_unlisted <- unlist(as(tmp, "GRangesList"))
              df_dena <- as.data.frame(unname(coordinate_dena_unlisted[,c(0,2,4,5)]))[,c(1:3,5,6,7,8)]
              
              tmp_rownames <- paste0(df_dena[, 6], "_", df_dena[, 7], "_", df_dena[, 5])
              dup_names <- names(which(table(tmp_rownames) > 1))
              ind_dup <- which(tmp_rownames %in% dup_names)
              ind_dup_rm <- ind_dup[which(duplicated(df_dena[ind_dup, c(5,6,7)]))]
              if (length(ind_dup_rm) > 0) {
                df_dena <- df_dena[-ind_dup_rm, ]
              }
              rownames(df_dena) <- paste0(df_dena[, 6], "_", df_dena[, 7], "_", df_dena[, 5])
              rownames(dena) <- paste0(dena[, 2], "_", dena[, 3], "_", dena[, 1])
              
              df_dena$Status <- dena[rownames(df_dena), 4]
              df_dena$Mod.Ratio <- dena[rownames(df_dena), 5]
              df_dena_final <- df_dena[,c(1,2,3,4,8,9)]
              colnames(df_dena_final) <- c("Chr", "Start", "End", "Strand", "Status", "Mod.Ratio")
              write.table(df_dena_final, file = output_file, quote = F, sep = "\t", row.names = F)
            }
              else {message(paste0(tool,"'s output files don't exist."))}
            }
            if (!file.exists(output_file)) {
              output_dena()
            }
            },
            drummer = {output_drummer <- function(){if (length(list.files(path = path_folder, pattern = "multiple_comp.txt", recursive = T)) != 0){
              all_data_drummer <- list.files(path = path_folder, pattern = "multiple_comp.txt", recursive = T)
              data_drummer <- data.frame()
              for (file in all_data_drummer) {
                if (file.info(paste0(path_folder, "/", file))$size != 0){
                  table <- read.table(paste0(path_folder, "/", file), header = TRUE)
                  data_drummer <- rbind(data_drummer, table)
                }
              }
              drummer <- data_drummer[, c("transcript_id", "position", "position", "max.G_padj", "max.G_padj")]
              #drummer <- data_drummer[, c(1,2,2,5,5)] #ok if -m is not set to TRUE in DRUMMER command
              #drummer <- data_drummer[, c(1,2,2,7,7)]
              drummer$position <- drummer$position
              drummer$position.1 <- drummer$position.1
              drummer$max.G_padj <- rep("Mod", length(drummer$max.G_padj))
              drummer$strand <- rep("*", nrow(drummer))
              drummer <- drummer[,c(1,2,3,6,4,5)]
              colnames(drummer) <- c("Chr", "Start", "End", "Strand", "Status", "Padj")
              write.table(drummer, file = output_file, quote = F, sep = "\t", row.names = F)
            }
              else {message(paste0(tool,"'s output files don't exist."))}
            }
            if (!file.exists(output_file)) {
              output_drummer()
            }
            },
            yanocomp = {output_yanocomp <- function(){if (file.exists(paste0(path_folder, "/", "yanocomp_output.bed"))
                                                          && file.info(paste0(path_folder, "/", "yanocomp_output.bed"))$size != 0){
              data_yanocomp <- read.table(paste0(path_folder, "/", "yanocomp_output.bed"))
              yanocomp <- data_yanocomp[, c(1:3,6,9,9)]
              yanocomp$V2 <- yanocomp$V2 + 2
              yanocomp$V3 <- yanocomp$V2
              yanocomp$V9 <- rep("Mod", length(yanocomp$V9))
              colnames(yanocomp) <- c("Chr", "Start", "End", "Strand", "Status", "Score")
              write.table(yanocomp, file = output_file, quote = F, sep = "\t", row.names = F)
            }
              else {message(paste0(tool,"'s output files don't exist."))}
            }
            if (!file.exists(output_file)) {
              output_yanocomp()
            }
            },
            nanocompore = {output_nanocompore <- function(){if (file.exists(paste0(path_folder, "/", "outnanocompore_results.tsv"))
                                                                && file.info(paste0(path_folder, "/", "outnanocompore_results.tsv"))$size != 0){
              data_nanocompore <- read.table(paste0(path_folder, "/", "outnanocompore_results.tsv"), header = TRUE)
              nanocompore <- data_nanocompore[, c(2,3,3,5,7,13,7)]
              nanocompore$genomicPos <- nanocompore$genomicPos + 2
              nanocompore$genomicPos.1 <- nanocompore$genomicPos.1 + 2
              nanocompore$Logit_LOR <- ifelse(nanocompore$Logit_LOR == "NC", NA, nanocompore$Logit_LOR)
              nanocompore$GMM_logit_pvalue <- ifelse(!is.nan(nanocompore$GMM_logit_pvalue) & !is.na(nanocompore$GMM_logit_pvalue) &
                                                       !is.nan(nanocompore$Logit_LOR) & !is.na(nanocompore$Logit_LOR) & 
                                                       nanocompore$GMM_logit_pvalue < filtering_parameter & 
                                                       abs(as.numeric(nanocompore$Logit_LOR)) > 0.5, "Mod", "Unmod")
              nanocompore <- nanocompore[which(nanocompore$GMM_logit_pvalue == "Mod"), ]
              nanocompore <- nanocompore[,c(1,2,3,4,5,7)]
              colnames(nanocompore) <- c("Chr", "Start", "End", "Strand", "Status", "Pvalue")
              write.table(nanocompore, file = output_file, quote = F, sep = "\t", row.names = F)
            }
              else {message(paste0(tool,"'s output files don't exist."))}
            }
            if (!file.exists(output_file)) {
              output_nanocompore()
            }
            },
            eligos = {output_eligos <- function(){if (file.exists(paste0(path_folder, "/", "minimap.sortG.1_vs_minimap.sortG.2_on_genome_combine.txt"))
                                                      && file.info(paste0(path_folder, "/", "minimap.sortG.1_vs_minimap.sortG.2_on_genome_combine.txt"))$size != 0){
              data_eligos <- read.table(paste0(path_folder, "/", "minimap.sortG.1_vs_minimap.sortG.2_on_genome_combine.txt"), header = TRUE, fill = TRUE)
              rows_eligos <- apply(data_eligos, 1, function(x){all(!is.na(x))})
              eligos <- data_eligos[which(rows_eligos),]
              eligos <- data_eligos[, c(1:4,18,16,18)]
              eligos$start_loc <- eligos$start_loc + 1
              eligos$end_loc <- eligos$end_loc
              eligos$adjPval <- ifelse(!is.nan(eligos$adjPval) & !is.na(eligos$adjPval) & 
                                         !is.nan(eligos$oddR) & !is.na(eligos$oddR) &
                                         eligos$adjPval < filtering_parameter & eligos$oddR > 1.2, "Mod", "Unmod")
              eligos <- eligos[which(eligos$adjPval == "Mod"), ]
              eligos <- eligos[,c(1,2,3,4,5,7)]
              colnames(eligos) <- c("Chr", "Start", "End", "Strand", "Status", "Padj")
              write.table(eligos, file = output_file, quote = F, sep = "\t", row.names = F)
            }
              else {message(paste0(tool,"'s output files don't exist."))}
            }
            if (!file.exists(output_file)) {
              output_eligos()
            }
            },
            mines = {output_mines <- function(){if (file.exists(paste0(path_folder, "/", "m6A_output_filename.bed"))
                                                    && file.info(paste0(path_folder, "/", "m6A_output_filename.bed"))$size != 0){
              data_mines <- read.table(paste0(path_folder, "/", "m6A_output_filename.bed")) 
              mines <- data_mines[, c(1:3,7,7)]
              mines$V2 <- mines$V2 + 1
              mines$V3 <- mines$V3
              mines$V7 <- rep("Mod", length(mines$V7))
              # Creation Edb Database from genome GTF
              EnsDb <- suppressWarnings(suppressMessages(ensDbFromGtf(gtf = genome_gtf)))
              edb <- EnsDb(EnsDb)
              # Lift-over + output bed
              test_mines <- IRanges(start = mines[,2], end = mines[,3], names = c(mines[,1]))
              
              num_rows_chunk <- 1
              mc.cores <- as.numeric(mccores)
              if (length(test_mines) < num_rows_chunk) {
                test_mines_split <- list(test_mines)
              } else {
                test_mines_split <- split(test_mines, rep(seq(from = 1, to = ceiling(length(test_mines)/num_rows_chunk)), each = num_rows_chunk)[1:length(test_mines)])
              }
              
              tmp1 <- vector(mode = "list", length = length(test_mines_split))
              names(tmp1) <- 1:length(test_mines_split)
              tmp <- vector(mode = "list", length = length(test_mines_split))
              names(tmp) <- 1:length(test_mines_split)
              ind_retry <- 1:length(test_mines_split)
              while(any(unlist(lapply(tmp, is.null)))) {
                cat(sprintf("Starting new iteration for mines; %d sites missing\n", length(which(unlist(lapply(tmp, is.null))))))
                tmp1 <- tmp1[ind_retry]
                tmp1 <- mclapply(test_mines_split[ind_retry], function(x) {
                  tryCatch({
                    coordinate_mines_unlisted <- unlist(transcriptToGenome(x, edb))
                    return(coordinate_mines_unlisted)
                  }, warning = function(w) {
                    print("Warning")
                    return(NULL)
                  }, error = function(e) {
                    print("Error")
                    return(NULL)
                  }
                  )}, mc.cores = mc.cores)
                ind_retry <- names(which(unlist(lapply(tmp1, function(x) is.null(x)))))
                ind_ok <- names(which(unlist(lapply(tmp1, function(x) !is.null(x)))))
                tmp[ind_ok] <- tmp1[ind_ok]
                if (length(ind_retry) > 0) {
                  tmp1 <- tmp1[ind_retry]
                }
              }
              
              coordinate_mines_unlisted <- unlist(as(tmp, "GRangesList"))
              
              df_mines <- as.data.frame(unname(coordinate_mines_unlisted[,c(0,2,4,5)]))[,c(1:3,5,6,7,8)]
              
              tmp_rownames <- paste0(df_mines[, 6], "_", df_mines[, 7], "_", df_mines[, 5])
              dup_names <- names(which(table(tmp_rownames) > 1))
              ind_dup <- which(tmp_rownames %in% dup_names)
              ind_dup_rm <- ind_dup[which(duplicated(df_mines[ind_dup, c(5,6,7)]))]
              if (length(ind_dup_rm) > 0) {
                df_mines <- df_mines[-ind_dup_rm, ]
              }
              rownames(df_mines) <- paste0(df_mines[, 6], "_", df_mines[, 7], "_", df_mines[, 5])
              rownames(mines) <- paste0(mines[, 2], "_", mines[, 3], "_", mines[, 1])
              df_mines$Status <- mines[rownames(df_mines), 4]
              df_mines$Ratiomod<- mines[rownames(df_mines), 5]
              df_mines_final <- df_mines[,c(1,2,3,4,8)]
              colnames(df_mines_final) <- c("Chr", "Start", "End", "Strand", "Status")
              write.table(df_mines_final, file = output_file, quote = F, sep = "\t", row.names = F)
            }
              else {message(paste0(tool,"'s output files don't exist."))}
            }
            if (!file.exists(output_file)) {
              output_mines()
            }
            },
            epinanoErr = {output_epinano_error <- function(){if (file.exists(paste0(path_folder, "/plus/", "diffErr.delta-sum_err.prediction.csv")) 
                                                                 || file.exists(paste0(path_folder, "/minus/", "diffErr.delta-sum_err.prediction.csv"))){
              if (file.exists(paste0(path_folder, "/plus/", "diffErr.delta-sum_err.prediction.csv")) 
                  && file.exists(paste0(path_folder, "/minus/", "diffErr.delta-sum_err.prediction.csv"))){
                data_epinanoerr_plus <- read.table(paste0(path_folder, "/plus/", "diffErr.delta-sum_err.prediction.csv"), header = TRUE, sep = ",")
                data_epinanoerr_minus <- read.table(paste0(path_folder, "/minus/", "diffErr.delta-sum_err.prediction.csv"), header = TRUE, sep = ",")
                data_epinanoerr <- rbind(data_epinanoerr_plus, data_epinanoerr_minus)
              }
              else if (file.exists(paste0(path_folder, "/plus/", "diffErr.delta-sum_err.prediction.csv"))){
                data_epinanoerr <- read.table(paste0(path_folder, "/plus/", "diffErr.delta-sum_err.prediction.csv"), header = TRUE, sep = ",")
              }
              else {data_epinanoerr <- read.table(paste0(path_folder, "/minus/", "diffErr.delta-sum_err.prediction.csv"), header = TRUE, sep = ",")}
              
              epinanoerr <- data.frame("Chr" = sapply(data_epinanoerr$chr_pos, function(x){return(strsplit(x, split = " ")[[1]][1])}),
                                       "Start" = sapply(data_epinanoerr$chr_pos, function(x){return(as.numeric(strsplit(x, split = " ")[[1]][2]))}),
                                       "End" = sapply(data_epinanoerr$chr_pos, function(x){return(as.numeric(strsplit(x, split = " ")[[1]][2]))}),
                                       "Strand" = sapply(data_epinanoerr$chr_pos, function(x){return(strsplit(x, split = " ")[[1]][4])}),
                                       "Status" = data_epinanoerr$z_score_prediction,
                                       "Delta sum err" = data_epinanoerr$delta_sum_err
              )
              epinanoerr <- epinanoerr[which(epinanoerr$Status == "mod"), ]
              write.table(epinanoerr, file = output_file, quote = F, sep = "\t", row.names = F)
            }
              else {message(paste0(tool,"'s output files don't exist."))}
            }
            if (!file.exists(output_file)) {
              output_epinano_error()
            }
            },
            epinanoSvm = {output_epinano_svm <- function(){if (file.exists(paste0(path_folder, "/", "plus_mod_prediction.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv")) 
                                                               || file.exists(paste0(path_folder, "/", "minus_mod_prediction.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv"))){
              if (file.exists(paste0(path_folder, "/", "plus_mod_prediction.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv")) 
                  && file.exists(paste0(path_folder, "/", "minus_mod_prediction.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv"))){
                data_epinanosvm_plus <- read.table(paste0(path_folder, "/", "plus_mod_prediction.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv"), sep=",") 
                data_epinanosvm_minus <- read.table(paste0(path_folder, "/", "minus_mod_prediction.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv"), sep=",")
                data_epinanosvm <- rbind(data_epinanosvm_plus, data_epinanosvm_minus)
              }
              else if (file.exists(paste0(path_folder, "/", "plus_mod_prediction.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv"))){
                data_epinanosvm <- read.table(paste0(path_folder, "/", "plus_mod_prediction.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv"), sep = ",")
              }
              else {data_epinanosvm <- read.table(paste0(path_folder, "/", "minus_mod_prediction.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv"), sep = ",")}
              
              epinanosvm <- data.frame("Kmer" = data_epinanosvm$V1,
                                       "Chr" = data_epinanosvm$V3,
                                       "Start" = sapply(data_epinanosvm$V2, function(x){return(strsplit(x, split = "\\-")[[1]][1])}),
                                       "End" = sapply(data_epinanosvm$V2, function(x){return(strsplit(x, split = "\\-")[[1]][1])}),
                                       "Strand" = data_epinanosvm$V4,
                                       "Status" = ifelse(data_epinanosvm$V28 > filtering_parameter, "Mod", "Unmod"), 
                                       "ProbM" = data_epinanosvm$V28
              )
              epinanosvm$Start <- as.numeric(epinanosvm$Start) + 2
              epinanosvm$End <- as.numeric(epinanosvm$End) + 2
              epinanosvm <- epinanosvm[which(epinanosvm$Kmer %in% rrach), ]
              epinanosvm <- epinanosvm[which(epinanosvm$Status == "Mod"), ]
              epinanosvm <- epinanosvm[, 2:7]
              write.table(epinanosvm, file = output_file, quote = F, sep = "\t", row.names = F)
            }
              else {message(paste0(tool,"'s output files don't exist."))}
            }
            if (!file.exists(output_file)) {
              output_epinano_svm()
            }
            },
            xpore = {output_xpore <- function(){if (file.exists(paste0(path_folder, "/", "diffmod.table"))
                                                    && file.info(paste0(path_folder, "/", "diffmod.table"))$size != 0){
              data_xpore <- read.table(paste0(path_folder, "/", "diffmod.table"), header = TRUE, sep=",") 
              xpore <- data.frame("GeneID" = data_xpore[,1],
                                  "Start" = data_xpore[,2] + 2,
                                  "End" = data_xpore[,2] + 2,
                                  "Status" = p.adjust(data_xpore[,5], method = "BH"), # Adviced to use FDR instead of pvalue
                                  "FDR" = p.adjust(data_xpore[,5], method = "BH")
              )
              xpore$Status <- ifelse(!is.nan(xpore$Status) & !is.na(xpore$Status) & xpore$Status < filtering_parameter, "Mod", "Unmod")
              xpore <- xpore[which(xpore$Status == "Mod"), ]
              # Add strand and chromosome comparing gene ID with genome bed file
              genome <- read.table(genome_bed, header = FALSE, sep="\t")
              rownames(genome) <- genome[,4]
              xpore$Strand <- genome[xpore$GeneID, "V6"]
              xpore$Chr <- genome[xpore$GeneID, "V1"]
              xpore <- xpore[,c(7,2,3,6,4,5)]
              write.table(xpore, file = output_file, quote = F, sep = "\t", row.names = F)
            }
              else {message(paste0(tool,"'s output files don't exist."))}
            }
            if (!file.exists(output_file)) {
              output_xpore()
            }
            },
            nanodoc = {output_nanodoc <- function(){if (length(list.files(path = path_folder, pattern = "*.txt")) != 0){
              txt_files_ls <- list.files(path = path_folder, pattern="*.txt")
              txt_files_df <- lapply(txt_files_ls, function(x) {if (file.info(paste0(path_folder, "/", x))$size != 0){
                table <- read.table(file = paste0(path_folder, "/", x), sep = "\t")
                table <- cbind(table, x) 
                table
              } 
              })
              if (length(txt_files_df) != 0){
                data_nanodoc <- do.call("rbind", lapply(txt_files_df, as.data.frame))
                nanodoc <- data_nanodoc[,c(1,1,12,12,13)] 
                nanodoc$V1 <- nanodoc$V1
                nanodoc$V1.1 <- nanodoc$V1.1
                nanodoc$V12 <- ifelse(!is.nan(nanodoc$V12) & !is.na(nanodoc$V12) & nanodoc$V12 > filtering_parameter, "Mod", "Unmod")
                nanodoc <- nanodoc[which(nanodoc$V12 == "Mod"), ]
                nanodoc$x <- str_extract(nanodoc$x, "(chr)[0-9]+|[IVX]{1,3}")
                nanodoc$strand <- rep("*", nrow(nanodoc))
                nanodoc <- nanodoc[,c(5,1,2,6,3,4)]
                colnames(nanodoc) <- c("Chr", "Start", "End", "Strand", "Status", "Score")
                if (nrow(nanodoc) > 0) {
                  write.table(nanodoc, file = output_file, quote = F, sep = "\t", row.names = F)
                }
              }
              else {message(paste0(tool,"'s output files don't exist."))}
            }
              else {message(paste0(tool,"'s output files don't exist."))}
            }
            if (!file.exists(output_file)) {
              output_nanodoc()
            }
            },
            nanom6a = {output_nanom6a <- function(){if (length(list.files(path = path_folder, pattern = "ratio.*.tsv")) != 0){
              files <- list.files(path = path_folder, pattern = "ratio.*.tsv")
              for (file in files){
                if (file.info(paste0(path_folder, "/", file))$size != 0){
                  data_nanom6a <- read.table(paste0(path_folder, "/", file), sep="\t", fill = T, col.names = 1:100)
                  nanom6a <- data.frame()
                  for (row in 1:nrow(data_nanom6a)) {
                    for (col in 1:length(data_nanom6a[row, ])) {
                      if (col == 1){
                        chr <- strsplit(data_nanom6a[row, col], split = "\\|")[[1]][2]
                        geneID <- strsplit(data_nanom6a[row, col], split = "\\|")[[1]][1]
                      }
                      else {
                        if(data_nanom6a[row, col] != "" && !is.na(data_nanom6a[row, col])){
                          start <- as.numeric(strsplit(data_nanom6a[row, col], split = "\\|")[[1]][1])
                          end <- as.numeric(strsplit(data_nanom6a[row, col], split = "\\|")[[1]][1])
                          mod_ratio <- strsplit(data_nanom6a[row, col], split = "\\|")[[1]][4]
                          x <- c(chr, geneID, start, end, mod_ratio, mod_ratio)
                          nanom6a <- rbind(nanom6a, x)
                        }
                      }
                    }
                  }
                  # Add strand comparing gene ID with genome bed file
                  genome <- read.table(genome_bed, header = FALSE, sep="\t")
                  rownames(genome) <- genome[,4]
                  nanom6a$Strand <- genome[nanom6a[,2], "V6"]
                  nanom6a <- nanom6a[,c(1,3,4,7)]
                  colnames(nanom6a) <- c("Chr", "Start", "End", "Strand")
                  nanom6a$Status <- rep("Mod", nrow(nanom6a))
                  write.table(nanom6a, file = paste0(output_file, "_", file, ".bed"), quote = F, sep = "\t", row.names = F)
                }
                else {message(paste0(tool,"'s output files is empty."))}
              }
            }
              else {message(paste0(tool,"'s output files don't exist."))}
            }
            if (!file.exists(output_file)) {
              output_nanom6a()
            }
            },
            tomboComparison = {output_tomboComparison <- function(){if (file.exists(paste0(path_folder, "/", "sample.level_samp_comp_detect.statistic.plus.wig"))
                                                                        && file.info(paste0(path_folder, "/", "sample.level_samp_comp_detect.statistic.plus.wig"))$size != 0){
              data_tombo <- read.table(paste0(path_folder, "/", "sample.level_samp_comp_detect.statistic.plus.wig"), fill = T, header = T)
              #tombo <- data.frame()
              tombo <- matrix(data = NA, nrow = dim(data_tombo)[1], ncol = 5)
              counter <- 1
              for (row in 1:nrow(data_tombo)) {
                for (col in 1:2) {
                  if (is.na(as.numeric(data_tombo[row,col]))){
                    if (col == 2){
                      transcriptID <- str_extract(data_tombo[row,col], "[A-Z0-9]{6,}(-[A-Z])?(-mRNA)?")
                    }
                  }
                  else if (col == 1){
                    start <- as.numeric(data_tombo[row, col]) 
                    end <- as.numeric(data_tombo[row, col])
                  }
                  else if (col == 2){
                    pvalue <- data_tombo[row, col]
                    x <- c(transcriptID, start, end, pvalue, pvalue)
                    #tombo <- rbind(tombo, x)
                    tombo[counter, ] <- x
                    counter <- counter + 1 
                  }
                }
              }
              tombo <- as.data.frame(tombo)
              tombo[,4] <- ifelse(!is.nan(tombo[,4]) & !is.na(tombo[,4]) & 10**(-as.numeric(tombo[,4])) < filtering_parameter, "Mod", "Unmod")
              tombo <- tombo[which(tombo[,4] == "Mod"), ]
              # Creation Edb Database from genome GTF
              EnsDb <- suppressWarnings(suppressMessages(ensDbFromGtf(gtf = genome_gtf)))
              edb <- EnsDb(EnsDb)
              # Lift-over + Creation of bed file
              test_tombo <- IRanges(start = as.numeric(tombo[,2]), end = as.numeric(tombo[,3]), names = c(tombo[,1]))
              
              num_rows_chunk <- 1
              mc.cores <- as.numeric(mccores)
              if (length(test_tombo) < num_rows_chunk) {
                test_tombo_split <- list(test_tombo)
              } else {
                test_tombo_split <- split(test_tombo, rep(seq(from = 1, to = ceiling(length(test_tombo)/num_rows_chunk)), each = num_rows_chunk)[1:length(test_tombo)])
              }
              
              tmp1 <- vector(mode = "list", length = length(test_tombo_split))
              names(tmp1) <- 1:length(test_tombo_split)
              tmp <- vector(mode = "list", length = length(test_tombo_split))
              names(tmp) <- 1:length(test_tombo_split)
              ind_retry <- 1:length(test_tombo_split)
              while(any(unlist(lapply(tmp, is.null)))) {
                cat(sprintf("Starting new iteration for Tombo; %d sites missing\n", length(which(unlist(lapply(tmp, is.null))))))
                tmp1 <- tmp1[ind_retry]
                tmp1 <- mclapply(test_tombo_split[ind_retry], function(x) {
                  tryCatch({
                    coordinate_tombo_unlisted <- unlist(transcriptToGenome(x, edb))
                    return(coordinate_tombo_unlisted)
                  }, warning = function(w) {
                    print("Warning")
                    return(NULL)
                  }, error = function(e) {
                    print("Error")
                    return(NULL)
                  }
                  )}, mc.cores = mc.cores)
                
                ind_retry <- names(which(unlist(lapply(tmp1, function(x) is.null(x)))))
                ind_ok <- names(which(unlist(lapply(tmp1, function(x) !is.null(x)))))
                tmp[ind_ok] <- tmp1[ind_ok]
                if (length(ind_retry) > 0) {
                  tmp1 <- tmp1[ind_retry]
                }
              }
              coordinate_tombo_unlisted <- unlist(as(tmp, "GRangesList"))
              df_tombo <- as.data.frame(unname(coordinate_tombo_unlisted[,c(0,2,4,5)]))
              df_tombo <- df_tombo[,c(1:3,5,6,7,8)]
              names_df_tombo <- paste0(df_tombo[, 6], "_", df_tombo[, 7], "_", df_tombo[, 5])
              rownames(tombo) <- paste0(tombo[, 2], "_", tombo[, 3], "_", tombo[, 1])
              df_tombo$Status <- tombo[names_df_tombo, 4]
              df_tombo$Pvalue <- 10**(-as.numeric(tombo[names_df_tombo, 5])) # Parameter of filtering is Pvalue not -log10(Pvalue)
              df_tombo_final <- df_tombo[,c(1,2,3,4,8,9)]
              colnames(df_tombo_final) <- c("Chr", "Start", "End", "Strand", "Status", "Pvalue")
              df_tombo_final$Start <- df_tombo_final$Start
              df_tombo_final$End <- df_tombo_final$End
              write.table(df_tombo_final, file = output_file, quote = F, sep = "\t", row.names = F)
            }
              else {message(paste0(tool,"'s output files don't exist."))}
            }
            if (!file.exists(output_file)) {
              output_tomboComparison()
            }
            },
            m6anet = {output_m6anet <- function(){if (file.exists(paste0(path_folder, "/", "data.result.csv"))
                                                      && file.info(paste0(path_folder, "/", "data.result.csv"))$size != 0){
              data_m6anet <- read.table(paste0(path_folder, "/", "data.result.csv"), header = TRUE, sep=",") 
              m6anet <- data.frame("TranscriptID" = data_m6anet[,1],
                                   "Start" = data_m6anet[,2] + 2,
                                   "End" = data_m6anet[,2] + 2,
                                   "Status" = data_m6anet[,4],
                                   "Prob_mod" = data_m6anet[,4]
              )
              m6anet$Status <- ifelse(!is.nan(m6anet$Status) & !is.na(m6anet$Status) & m6anet$Status > filtering_parameter, "Mod", "Unmod")
              m6anet <- m6anet[which(m6anet$Status == "Mod"), ]
              # Creation Edb Database from genome GTF
              EnsDb <- suppressWarnings(suppressMessages(ensDbFromGtf(gtf = genome_gtf)))
              edb <- EnsDb(EnsDb)
              # Lift-over + output bed
              test_m6anet <- IRanges(start = m6anet[,2], end = m6anet[,3], names = c(m6anet[,1]))
              num_rows_chunk <- 1
              mc.cores <- as.numeric(mccores)
              if (length(test_m6anet) < num_rows_chunk) {
                test_m6anet_split <- list(test_m6anet)
              } else {
                test_m6anet_split <- split(test_m6anet, rep(seq(from = 1, to = ceiling(length(test_m6anet)/num_rows_chunk)), each = num_rows_chunk)[1:length(test_m6anet)])
              }
              
              tmp1 <- vector(mode = "list", length = length(test_m6anet_split))
              names(tmp1) <- 1:length(test_m6anet_split)
              tmp <- vector(mode = "list", length = length(test_m6anet_split))
              names(tmp) <- 1:length(test_m6anet_split)
              ind_retry <- 1:length(test_m6anet_split)
              while(any(unlist(lapply(tmp, is.null)))) {
                cat(sprintf("Starting new iteration for m6Anet; %d sites missing\n", length(which(unlist(lapply(tmp, is.null))))))
                tmp1 <- tmp1[ind_retry]
                tmp1 <- mclapply(test_m6anet_split[ind_retry], function(x) {
                  tryCatch({
                    coordinate_m6anet_unlisted <- unlist(transcriptToGenome(x, edb))
                    return(coordinate_m6anet_unlisted)
                  }, warning = function(w) {
                    print("Warning")
                    return(NULL)
                  }, error = function(e) {
                    print("Error")
                    return(NULL)
                  }
                  )}, mc.cores = mc.cores)
                ind_retry <- names(which(unlist(lapply(tmp1, function(x) is.null(x)))))
                ind_ok <- names(which(unlist(lapply(tmp1, function(x) !is.null(x)))))
                tmp[ind_ok] <- tmp1[ind_ok]
                if (length(ind_retry) > 0) {
                  tmp1 <- tmp1[ind_retry]
                }
              }
              
              coordinate_m6anet_unlisted <- unlist(as(tmp, "GRangesList"))
              #coordinate_m6anet_unlisted <- unlist(transcriptToGenome(test_m6anet, edb))
              df_m6anet <- as.data.frame(unname(coordinate_m6anet_unlisted[,c(0,2,4,5)]))[,c(1:3,5,6,7,8)]
              
              tmp_rownames <- paste0(df_m6anet[, 6], "_", df_m6anet[, 7], "_", df_m6anet[, 5])
              dup_names <- names(which(table(tmp_rownames) > 1))
              ind_dup <- which(tmp_rownames %in% dup_names)
              ind_dup_rm <- ind_dup[which(duplicated(df_m6anet[ind_dup, c(5,6,7)]))]
              if (length(ind_dup_rm) > 0) {
                df_m6anet <- df_m6anet[-ind_dup_rm, ]
              }
              rownames(df_m6anet) <- paste0(df_m6anet[, 6], "_", df_m6anet[, 7], "_", df_m6anet[, 5])
              
              rownames(m6anet) <- paste0(m6anet[, 2], "_", m6anet[, 3], "_", m6anet[, 1])
              
              df_m6anet$Status <- m6anet[rownames(df_m6anet), 4]
              df_m6anet$Prob_Mod <- m6anet[rownames(df_m6anet), 5]
              df_m6anet_final <- df_m6anet[,c(1,2,3,4,8,9)]
              colnames(df_m6anet_final) <- c("Chr", "Start", "End", "Strand", "Status", "Prob_mod")
              write.table(df_m6anet_final, file = output_file, quote = F, sep = "\t", row.names = F)
            }
              else {message(paste0(tool,"'s output files don't exist."))}
            }
            if (!file.exists(output_file)) {
              output_m6anet()
            }
            },
            stop("Enter a valid tool as input!")
    )}
}

# Definining a set of parameters used to filter the results for those tools which give as output all the sites 
rrach <- c("AAACA","AAACT","AAACC","GAACA","GAACT","GAACC","GGACA","GGACT","GGACC","GAACA","GAACT","GAACC")
tools <- c("dena", "drummer", "differr", "yanocomp", "nanocompore", "eligos", "mines"
           , "epinanoErr", "epinanoSvm", "xpore", "nanodoc", "nanom6a", "tomboComparison", "m6anet")

pathTools <- c(pathdena, pathdrummer, pathdifferr, pathyanocomp, pathnanocompore, patheligos
               , pathmines, pathepinanoError, pathepinanoSVM, pathxpore, pathnanodoc, pathnanom6a
               , pathtomboComparison, pathm6anet)

default <- c(0.1, NA, NA, NA, 0.01, 0.001, NA, NA, 0.5, 0.05, 0.02, NA, 0.05, 0.9)
relaxed <- c(0, NA, NA, NA, 1, 1, NA, NA, 0, 1, 0, NA, 1, 0)
value <- rep(threshold, length(default))

names(pathTools) <- names(default) <- names(relaxed) <- names(value) <- tools

parameters_list <- list("default" = unname(default), "relaxed" = unname(relaxed), "value" = unname(value))

if (grepl(x = threshold, pattern = "default")) {
  parameters <- parameters_list$default
} else if (grepl(x = threshold, pattern = "relaxed")) {
  parameters <- parameters_list$relaxed
} else {
  parameters <- as.numeric(parameters_list$value)
}

# Data frame containing the results from all the tools
results_df <- data.frame(row.names = tools,
                         path_folder = unname(pathTools),
                         parameter = parameters)

# Looping through all the elements in the data frame applying the output_processing function
for (x in row.names(results_df)){output_processing(tool = x, 
                                                   path_folder = paste0(path, "/", results_df[x, "path_folder"]), 
                                                   output_file = paste0(resultsFolder, "/", x, "_", "output.bed"), 
                                                   genome_gtf = genomegtf,
                                                   genome_bed = genomebed,
                                                   filtering_parameter = results_df[x, "parameter"])}
