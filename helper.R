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

GetDegs <- function(df.fc, cutoff) {
  df.degs <- NULL
  df.degs <- df.fc[abs(df.fc$LOG2FC) >= cutoff, ]
  return(df.degs)
}

PrintStats <- function(dflist) {
  cat("[GROUP]\t[-T]\t[+T]\t[NROW]\t[MAX]\t[MIN]\n")
  for(i in 1:length(dflist)) {
    df <- dflist[[i]]
    cat(paste(names(dflist)[i], "_", cutoff[i], sep = ""), "\t")
    cat(nrow(df[df$LOG2FC <= cutoff[i]*(-1), ]), "\t")
    cat(nrow(df[df$LOG2FC >= cutoff[i], ]), "\t")
    cat(nrow(df), "\t")
    cat(max(df$LOG2FC), "\t")
    cat(min(df$LOG2FC), "\n")
  }
}

GetInters <- function(df.list) {
  list.int <- list(
    RA_METS  = Reduce(intersect, list((df.list$RA)$Symbol, (df.list$METS)$Symbol)),
    RA_CAD   = Reduce(intersect, list((df.list$RA)$Symbol, (df.list$CAD)$Symbol)),
    RA_T2D   = Reduce(intersect, list((df.list$RA)$Symbol, (df.list$T2D)$Symbol)),
    METS_T2D = Reduce(intersect, list((df.list$METS)$Symbol, (df.list$T2D)$Symbol)),
    METS_CAD = Reduce(intersect, list((df.list$METS)$Symbol, (df.list$CAD)$Symbol)),
    CAD_T2D  = Reduce(intersect, list((df.list$CAD)$Symbol, (df.list$T2D)$Symbol)),
    
    CAD_RA_T2D   = Reduce(intersect, list((df.list$CAD)$Symbol, (df.list$RA)$Symbol, (df.list$T2D)$Symbol)),
    CAD_RA_METS  = Reduce(intersect, list((df.list$CAD)$Symbol, (df.list$RA)$Symbol, (df.list$METS)$Symbol)),
    CAD_T2D_METS = Reduce(intersect, list((df.list$CAD)$Symbol, (df.list$T2D)$Symbol, (df.list$METS)$Symbol)),
    RA_T2D_METS  = Reduce(intersect, list((df.list$RA)$Symbol, (df.list$T2D)$Symbol, (df.list$METS)$Symbol))
  )
  return(list.int)
}

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

FilterPlinks <- function(df.plinks, string.map) {
  df.result <- df.plinks
  cat("Total plinks =", nrow(df.result),"\n")
  
  # Omit insignificant links
  df.result <- df.result[df.result$combined_score >= 700, ]
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

GetMSGs <- function(df.raw.all) {
  if(file.exists("EXPORT/MSGENES_14422.csv")) {
    print("Reading from existing file!")
    list.msgs <- read.csv("EXPORT/MSGENES_14422.csv", stringsAsFactors = FALSE, header = TRUE)$Symbol
  }
  else {
    string.map <- GetStringMap(df.raw.all) # list
    df.plinks <- as.data.frame(fread("IMPORT/PLINKS.tsv", header = TRUE, sep = ' '))[, c(1, 2, 16)]
    df.plinks <- FilterPlinks(df.plinks, string.map)
    df.glinks <- GetGlinks(df.plinks, string.map)
    list.msgs <- unique(as.vector(rbind(df.glinks$gene1, df.glinks$gene2))) # mapped sig genes
    write.csv(as.data.frame(cbind(Symbol = list.msgs)), file = "EXPORT/MSGENES_14422.csv", row.names = FALSE)
  }
  return(list.msgs)
}