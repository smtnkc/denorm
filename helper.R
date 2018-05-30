library(readr)
library(rstudioapi)
library(data.table)

rm(list=ls())
setwd(dirname(getSourceEditorContext()$path))
dir()
dir.csv <- "CSV/"
dir.raw <- "EXPORT/FINAL/"
dir.deg <- "EXPORT/DEG/"
f.names <- colnames(as.data.frame(read.csv("SERIES.csv", header = TRUE, stringsAsFactors = FALSE)))[2:36]

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

source("denorm.R")