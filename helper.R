library(readr)
library(rstudioapi)
library(data.table) # fread
library(STRINGdb)
library(org.Hs.eg.db) # to get entrez ids
library(svMisc) # to show progress

rm(list=ls())
setwd(dirname(getSourceEditorContext()$path))
dir()
dir.csv <- "CSV/"
dir.raw <- "EXPORT/FINAL/"
dir.deg <- "EXPORT/DEG/"
f.names <- colnames(as.data.frame(read.csv("IMPORT/SERIES.csv", header = TRUE, stringsAsFactors = FALSE)))[2:36]

FormatDf <- function(df) {
  df[is.na(df)] <- ""
  df[df == "-"] <- ""
  df[df == "."] <- ""
  df[df == "Null"] <- ""
  df[df == "--Null"] <- ""
  df[df == "EMPTY"] <- ""
  df[df == "--empty"] <- ""
  df[df == "--unknown"] <- ""
  df[df == "Dye Marker"] <- ""
  return (df)
}

GroupByCol <- function(df, by, method) {
  df <- aggregate(df[, -by], list(toupper(df[, by])), method)
  colnames(df)[1] <- "Symbol"
  df$"F635.Median" <- round(df$"F635.Median", 0)
  df$"TrueVal" <- round(df$"TrueVal", 2)
  df$"SeriesVal" <- round(df$"SeriesVal", 2)
  return(df)
}

GetRawVals <- function(f.names, states) {
  f.raw <- paste(dir.raw, f.names[1], ".csv", sep = "")
  df.result <- as.data.frame(read.csv(f.raw, header = TRUE, stringsAsFactors = FALSE)[, 1])
  colnames(df.result)[1] <- "Symbol"
  
  for(state in states) {
    f.raw <- paste(dir.raw, f.names[state], ".csv", sep = "")
    v.temp <- read.csv(f.raw, header = TRUE, stringsAsFactors = FALSE)[, 2]
    df.result[, f.names[state]] <- v.temp
  }
  return(df.result)
}

df.raw.all <- GetRawVals(f.names, c(1:35))

MapProbes <- function(df.gpl, df.map) {
  
  df.clean.map <- df.map[df.map$Name != "", ]
  df.clean.gpl <- df.gpl[df.gpl$Name != "", ]
  
  df.clean.map <- df.clean.map[order(df.clean.map$Name), ]
  df.clean.gpl <- df.clean.gpl[order(df.clean.gpl$Name), ]
  
  row.names(df.clean.map) <- NULL
  row.names(df.clean.gpl) <- NULL
  
  df.clean.gpl[, "Probe_Name"] <- ""
  
  unmapped <- c()
  for(i in 1:nrow(df.clean.map)) {
    if(df.clean.map[i, "Name"] == df.clean.gpl[i, "Name"])
      df.clean.gpl[i, "Probe_Name"] <- df.clean.map[i, "Probe_Name"]
    else
      unmapped <- cbind(unmapped, i)
  }
  cat(length(unmapped), "probes cannot be mapped!\n")
  return(df.clean.gpl)
}

########################################### deg.R

GetStringMap <- function(df.raw) {
  
  strdb <- STRINGdb$new(species = 9606, version = "10")
  
  cat("Number of proteins:", nrow(strdb$get_proteins()), "\n")
  df.map <- as.data.frame(select(org.Hs.eg.db,
                                 keys = as.vector(df.raw$Symbol),
                                 columns = c("ENTREZID", "ACCNUM", "SYMBOL"),
                                 keytype = "SYMBOL"))
  
  df.map <- strdb$map(df.map,
                      c("ENTREZID", "SYMBOL"),
                      takeFirst = TRUE,
                      removeUnmappedRows = FALSE,
                      quiet = FALSE)
  
  df.map <- df.map[, c("SYMBOL", "STRING_id")]
  df.map <- df.map[!duplicated(df.map$SYMBOL), ]
  df.map <- df.map[!is.na(df.map$STRING_id), ]
  df.map <- df.map[!duplicated(df.map$STRING_id), ]
  
  string.map <- list()
  for(i in 1:nrow(df.map)) {
    string.map[df.map[i, "STRING_id"]] <- df.map[i, "SYMBOL"]
  }
  return(string.map)
}

FilterPlinks <- function(df.plinks, string.map, comb_score_cutoff) {
  df.result <- df.plinks
  cat("Total plinks =", nrow(df.result),"\n")
  
  # Omit insignificant links
  df.result <- df.result[df.result$combined_score >= comb_score_cutoff, ]
  cat("Significant plinks =", nrow(df.result),"\n")
  
  # Omit unmapped proteins
  df.result <- df.result[df.result$protein1 %in% names(string.map), ]
  df.result <- df.result[df.result$protein2 %in% names(string.map), ]
  cat("Mapped significant plinks =", nrow(df.result),"\n")
  return(df.result)
}

GetGlinks <- function(df.plinks, string.map) {
  df.glinks <- df.plinks
  df.glinks$gene1 <- as.vector(as.character(string.map[df.glinks$protein1]))
  df.glinks$gene2 <- as.vector(as.character(string.map[df.glinks$protein2]))
  df.glinks <- df.glinks[, c("gene1", "gene2", "combined_score")]
  return(df.glinks)
}

GetMedian <- function(df) {
  df.result <- as.data.frame(df[, 1])
  colnames(df.result)[1] <- "Symbol"
  df.result[, "MedVal"] <- NA
  for(i in 1:nrow(df))
    df.result[i, "MedVal"] <- median(as.numeric(df[i, 2:ncol(df)]))
  
  return(df.result)
}

GetFC <- function(df.ctrl, df.disease) {
  df.result <- as.data.frame(df.ctrl[, 1])
  colnames(df.result)[1] <- "Symbol"
  df.result[, "LOG2FC"] <- NA
  
  for(i in 1:nrow(df.result))
    df.result[i , "LOG2FC"] <- round(log2(df.disease[i, "MedVal"]/df.ctrl[i, "MedVal"]), 2)
  
  return(df.result)
}

GetDegs <- function(fc_cutoff) {
  fname <- paste(dir.deg, "RA_DEGs_", fc_cutoff, ".csv", sep = "")
  if(file.exists(fname)) {
    cat("fc_cutoff =", fc_cutoff,"\nReading DEGs from existing files...\n")
    list.degs <- list(
      RA   = as.data.frame(read.csv(paste(dir.deg, "RA_DEGs_",   fc_cutoff, ".csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)),
      METS = as.data.frame(read.csv(paste(dir.deg, "METS_DEGs_", fc_cutoff, ".csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)),
      CAD  = as.data.frame(read.csv(paste(dir.deg, "CAD_DEGs_",  fc_cutoff, ".csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)),
      T2D  = as.data.frame(read.csv(paste(dir.deg, "T2D_DEGs_",  fc_cutoff, ".csv", sep = ""), header = TRUE, stringsAsFactors = FALSE))
    )
  }
  else {
    cat("fc_cutoff =", fc_cutoff,"\nGenerating DEGs...\n")
    list.medians <- list (
      CTRL = GetMedian(GetRawVals(f.names, c(1:9))),
      RA   = GetMedian(GetRawVals(f.names, c(10:15))),
      METS = GetMedian(GetRawVals(f.names, c(16:21))),
      CAD  = GetMedian(GetRawVals(f.names, c(22:27))),
      T2D  = GetMedian(GetRawVals(f.names, c(28:35)))
    )
    
    list.FC <- list (
      RA   = GetFC(list.medians$CTRL, list.medians$RA),
      METS = GetFC(list.medians$CTRL, list.medians$METS),
      CAD  = GetFC(list.medians$CTRL, list.medians$CAD),
      T2D  = GetFC(list.medians$CTRL, list.medians$T2D)
    )
    
    list.degs <- list(
      RA   = list.FC[["RA"]][abs(list.FC[["RA"]]$LOG2FC) >= fc_cutoff, ],
      METS = list.FC[["METS"]][abs(list.FC[["METS"]]$LOG2FC) >= fc_cutoff, ],
      CAD  = list.FC[["CAD"]][abs(list.FC[["CAD"]]$LOG2FC) >= fc_cutoff, ],
      T2D  = list.FC[["T2D"]][abs(list.FC[["T2D"]]$LOG2FC) >= fc_cutoff, ]
    )
    ExportDegs(list.degs, fc_cutoff)
    cat("DEG files have been exported for the next time!\n")
  }
  cat("DONE!\n")
  return(list.degs)
}

ExportDegs <- function(list.degs, fc_cutoff) {
  f.out <- paste(dir.deg, "RA_DEGs_", fc_cutoff, ".csv", sep = "")
  write.csv(list.degs$RA, file = f.out, row.names = FALSE)
  f.out <- paste(dir.deg, "METS_DEGs_", fc_cutoff,".csv", sep = "")
  write.csv(list.degs$METS, file = f.out, row.names = FALSE)
  f.out <- paste(dir.deg, "CAD_DEGs_", fc_cutoff, ".csv", sep = "")
  write.csv(list.degs$CAD, file = f.out, row.names = FALSE)
  f.out <- paste(dir.deg, "T2D_DEGs_", fc_cutoff, ".csv", sep = "")
  write.csv(list.degs$T2D, file = f.out, row.names = FALSE)
}

GetMSGenes <- function(df.raw.all, comb_score_cutoff) {
  fname <- paste("EXPORT/MS_GENES_", comb_score_cutoff, ".csv", sep = "")
  if(file.exists(fname)) {
    cat("Reading from existing file:", fname, "\n")
    list.msgenes <- read.csv(fname, stringsAsFactors = FALSE, header = TRUE)$Symbol
  }
  else {
    string.map <- GetStringMap(df.raw.all) # list
    df.plinks <- as.data.frame(fread("IMPORT/PLINKS.tsv", header = TRUE, sep = ' '))[, c(1, 2, 16)]
    df.plinks <- FilterPlinks(df.plinks, string.map, comb_score_cutoff)
    df.glinks <- GetGlinks(df.plinks, string.map)
    list.msgenes <- unique(as.vector(rbind(df.glinks$gene1, df.glinks$gene2))) # mapped sig genes
    write.csv(as.data.frame(cbind(Symbol = list.msgenes)), file = fname, row.names = FALSE)
  }
  return(list.msgenes)
}

PrintDegStats <- function(dflist, fc_cutoff) {
  cat("[GROUP]\t[-T]\t[+T]\t[NROW]\t[MAX]\t[MIN]\n")
  for(i in 1:length(dflist)) {
    df <- dflist[[i]]
    cat(paste(names(dflist)[i], "_", fc_cutoff, sep = ""), "\t")
    cat(nrow(df[df$LOG2FC <= fc_cutoff*(-1), ]), "\t")
    cat(nrow(df[df$LOG2FC >= fc_cutoff, ]), "\t")
    cat(nrow(df), "\t")
    cat(max(df$LOG2FC), "\t")
    cat(min(df$LOG2FC), "\n")
  }
}

########################################### cor.R

GetRawMSDegs <- function(df.raw.all, list.msdegs) {
  # raw vals for mapped sig genes for each group included CTRL
  list.raw.msdegs <- list(
    RA   = df.raw.all[df.raw.all$Symbol %in% list.msdegs$RA$Symbol,   c(1:16)],
    METS = df.raw.all[df.raw.all$Symbol %in% list.msdegs$METS$Symbol, c(1:10,17:22)],
    CAD  = df.raw.all[df.raw.all$Symbol %in% list.msdegs$CAD$Symbol,  c(1:10,23:28)],
    T2D  = df.raw.all[df.raw.all$Symbol %in% list.msdegs$T2D$Symbol,  c(1:10,29:36)]
  )
  
  # set first col as rownames for each group
  for(i in 1:length(list.raw.msdegs)) {
    row.names(list.raw.msdegs[[i]]) <- list.raw.msdegs[[i]]$Symbol
    list.raw.msdegs[[i]] <- list.raw.msdegs[[i]][, -1]
  }
  return(list.raw.msdegs)
}

GetCorMatrices <- function(list.raw.msds) {
  list.cors <- list()
  for(i in 1:length(list.raw.msds)) {
    df.cor <- round(cor(x = t(list.raw.msds[[i]]), y = NULL, method = c("pearson")),2)
    group <- names(list.raw.msds)[i]
    list.cors[[group]] <- df.cor
    fname = paste("EXPORT/COR/", group, ".csv", sep = "")
    write.csv(df.cor, file = fname, row.names = TRUE)
  }
  return(list.cors)
}

GetMSDegLinks <- function(df.raw.all, list.cors, comb_score_cutoff) {
  fname <- paste("EXPORT/MS_GLINKS_", comb_score_cutoff, ".csv", sep = "")
  if(file.exists(fname)) {
    cat("Reading from existing file:", fname, "\n")
    df.msglinks <- as.data.frame(read.csv(fname, stringsAsFactors = FALSE, header = TRUE))
  }
  else {
    string.map <- GetStringMap(df.raw.all) # list
    df.plinks <- as.data.frame(fread("IMPORT/PLINKS.tsv", header = TRUE, sep = ' '))[, c(1, 2, 16)]
    df.plinks <- FilterPlinks(df.plinks, string.map, comb_score_cutoff)
    df.msglinks <- GetGlinks(df.plinks, string.map)
    write.csv(df.msglinks, file = fname, row.names = FALSE)
  }
  rownames(df.msglinks) <- NULL
  
  list.msdeglinks <- list()
  for(i in 1:length(list.cors)) {
    df.temp <- df.msglinks[df.msglinks$gene1 %in% row.names(list.cors[[i]]), 1:2]
    df.temp <- df.temp[df.temp$gene2 %in% row.names(list.cors[[i]]), ]
    list.msdeglinks[[names(list.cors)[i]]] <- df.temp
  }
  return(list.msdeglinks)
}

GetCorEdges <- function(list.msdeglinks, list.cor.matrices) {
  list.pcor <- list()
  for(i in 1:length(list.msdeglinks)) {
    df.pcor <- list.msdeglinks[[i]]
    if(nrow(list.msdeglinks[[i]]) == 0) {
      df.pcor$pearson <- numeric()
      print(dim(df.pcor))
    }
    else {
      df.pcor[, "pearson"] <- NA
      # print(dim(df.pcor))
      df.matrix <- list.cor.matrices[[i]]
      
      for(j in 1:nrow(list.msdeglinks[[i]])) {
        corval <- df.matrix[list.msdeglinks[[i]]$gene1[j], list.msdeglinks[[i]]$gene2[j]]
        df.pcor[j, "pearson"] <- corval
      }
    }
    list.pcor[[names(list.msdeglinks)[i]]] <- df.pcor
  }
  return(list.pcor)
}

GetSigCorEdges <- function(list.cor.edges, p_val) {
  list.sig.cor.edges <- list()
  for(i in 1:length(list.cor.edges)) {
    df <- list.cor.edges[[i]]
    if(names(list.cor.edges)[i] == "T2D") {
      if(p_val == 0.05) coef <- 0.482
      else if(p_val == 0.01) coef <- 0.606
    }
    else {
      if(p_val == 0.05) coef <- 0.514
      else if(p_val == 0.01) coef <- 0.641
    }
    
    list.sig.cor.edges[[names(list.cor.edges)[i]]] <-
      df[abs(df$pearson) >= coef, ]
  }
  return(list.sig.cor.edges)
}

GetSigCorNodes <- function(list.sig.cor.edges) {
  list.sig.cor.genes <- list()
  for(i in 1:length(list.sig.cor.edges)) {
    df <- list.sig.cor.edges[[i]]
    genes <- unique(as.character(rbind(df$gene1, df$gene2)))
    list.sig.cor.genes[[names(list.sig.cor.edges)[i]]] <- genes
  }
  return(list.sig.cor.genes)
}

GetInterNodes <- function(df.list) {
  list.int <- list(
    RA_METS  = Reduce(intersect, list(df.list$RA, df.list$METS)),
    RA_CAD   = Reduce(intersect, list(df.list$RA, df.list$CAD)),
    RA_T2D   = Reduce(intersect, list(df.list$RA, df.list$T2D)),
    METS_T2D = Reduce(intersect, list(df.list$METS, df.list$T2D)),
    METS_CAD = Reduce(intersect, list(df.list$METS, df.list$CAD)),
    CAD_T2D  = Reduce(intersect, list(df.list$CAD, df.list$T2D)),
    
    CAD_RA_T2D   = Reduce(intersect, list(df.list$CAD, df.list$RA, df.list$T2D)),
    CAD_RA_METS  = Reduce(intersect, list(df.list$CAD, df.list$RA, df.list$METS)),
    CAD_T2D_METS = Reduce(intersect, list(df.list$CAD, df.list$T2D, df.list$METS)),
    RA_T2D_METS  = Reduce(intersect, list(df.list$RA, df.list$T2D, df.list$METS)),
    
    RA_T2D_METS_CAD  = Reduce(intersect, list(df.list$RA, df.list$T2D,
                                              df.list$METS, df.list$CAD))
  )
  return(list.int)
}

GetExtendedEdgeList <- function(df) {
  list.ext <- list()
  for(i in 1:nrow(df)) {
    pair <- c(df[i, "gene1"], df[i, "gene2"])
    list.ext <- append(list.ext, list(pair))
    # list.ext <- append(list.ext, list(rev(pair)))
  }
  list.ext <- unique(list.ext)
  return(list.ext)
}

GetInterEdges <- function(df.list) {
  eel_RA   <- GetExtendedEdgeList(df.list$RA)
  eel_METS <- GetExtendedEdgeList(df.list$METS)
  eel_CAD  <- GetExtendedEdgeList(df.list$CAD)
  eel_T2D  <- GetExtendedEdgeList(df.list$T2D)
  
  list.int <- list(
    RA_METS  = Reduce(intersect, list(eel_RA, eel_METS)),
    RA_CAD   = Reduce(intersect, list(eel_RA, eel_CAD)),
    RA_T2D   = Reduce(intersect, list(eel_RA, eel_T2D)),
    METS_T2D = Reduce(intersect, list(eel_METS, eel_T2D)),
    METS_CAD = Reduce(intersect, list(eel_METS, eel_CAD)),
    CAD_T2D  = Reduce(intersect, list(eel_CAD, eel_T2D)),
    
    CAD_RA_T2D   = Reduce(intersect, list(eel_CAD, eel_RA, eel_T2D)),
    CAD_RA_METS  = Reduce(intersect, list(eel_CAD, eel_RA, eel_METS)),
    CAD_T2D_METS = Reduce(intersect, list(eel_CAD, eel_T2D, eel_METS)),
    RA_T2D_METS  = Reduce(intersect, list(eel_RA, eel_T2D, eel_METS)),
    
    RA_T2D_METS_CAD  = Reduce(intersect, list(eel_RA, eel_T2D, eel_METS, eel_CAD))
  )
  return(list.int)
}

ExportEdgeLists <- function(list, p_val, comb_score_cutoff, fc_cutoff) {
  for(i in 1:length(list)) {
    fname = paste("EXPORT/EDGELIST/", names(list)[i], "_",
                  comb_score_cutoff, "_FC", fc_cutoff, "_P",
                  gsub("\\.", "", p_val), ".csv", sep = "")
    write.csv(list[[i]], fname, row.names = FALSE)
  }
}

ExportInterEdges <- function(list, p_val, comb_score_cutoff, fc_cutoff) {
  for(i in 1:length(list)) {
    fname = paste("EXPORT/INTERS/", names(list)[i], "_",
                  comb_score_cutoff, "_FC", fc_cutoff, "_P",
                  gsub("\\.", "", p_val), ".csv", sep = "")
    df <- data.frame(matrix(ncol = 2, nrow = 0))
    for(j in 1:length(list[[i]])) {
      df <- rbind(df, data.frame(gene1 = list[[i]][[j]][1], gene2 = list[[i]][[j]][2]))
    }
    colnames(df) <- c("gene1", "gene2")
    write.csv(df, fname, row.names = FALSE)
  }
}

ExportHubs <- function(list.hubs.ext) {
  df <- cbind(list.hubs.ext[[1]], list.hubs.ext[[2]], list.hubs.ext[[3]], list.hubs.ext[[4]])
  colnames(df) <- c("RA_Gene", "RA_Degree", "METS_Gene", "METS_Degree",
                    "CAD_Gene", "CAD_Degree", "T2D_Gene", "T2D_Degree")
  fname <- paste("EXPORT/HUBS/", "Degree_HUBS.csv", sep = "")
  write.csv(df, fname, row.names = FALSE)
}