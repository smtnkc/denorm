rm(list=ls())
setwd(dirname(getSourceEditorContext()$path))
dir()
dir.in <- "EXPORT/FINAL/"
dir.out <- "EXPORT/DEG/"
f.names <- colnames(as.data.frame(read.csv("SERIES.csv", header = TRUE, stringsAsFactors = FALSE)))[2:36]

################################################ MEDIANS

GetRawVals <- function(f.names, states) {
  f.in <- paste(dir.in, f.names[1], ".csv", sep = "")
  df.result <- as.data.frame(read.csv(f.in, header = TRUE, stringsAsFactors = FALSE)[, 1])
  colnames(df.result)[1] <- "Symbol"
  
  for(state in states) {
    f.in <- paste(dir.in, f.names[state], ".csv", sep = "")
    v.temp <- read.csv(f.in, header = TRUE, stringsAsFactors = FALSE)[, 2]
    df.result[, f.names[state]] <- v.temp
  }
  return(df.result)
}

df.ctrl <- GetRawVals(f.names, c(1:9))
df.ra <- GetRawVals(f.names, c(10:15))
df.mets <- GetRawVals(f.names, c(16:21))
df.cad <- GetRawVals(f.names, c(22:27))
df.t2d <- GetRawVals(f.names, c(28:35))

GetMedian <- function(df) {
  df.result <- as.data.frame(df[, 1])
  colnames(df.result)[1] <- "Symbol"
  df.result[, "MedVal"] <- NA
  for(i in 1:nrow(df))
    df.result[i, "MedVal"] <- median(as.numeric(df[i, 2:ncol(df)]))

  return(df.result)
}

df.ctrl <- GetMedian(df.ctrl)
df.ra <- GetMedian(df.ra)
df.mets <- GetMedian(df.mets)
df.cad <- GetMedian(df.cad)
df.t2d <- GetMedian(df.t2d)

################################################ FC

GetFC <- function(df.ctrl, df.disease) {
  df.result <- as.data.frame(df.ctrl[, 1])
  colnames(df.result)[1] <- "Symbol"
  df.result[, "LOG2FC"] <- NA
  
  for(i in 1:nrow(df.result))
    df.result[i , "LOG2FC"] <- round(log2(df.disease[i, "MedVal"]/df.ctrl[i, "MedVal"]), 2)

  return(df.result)
}

fc.ra <- GetFC(df.ctrl, df.ra)
fc.mets <- GetFC(df.ctrl, df.mets)
fc.cad <- GetFC(df.ctrl, df.cad)
fc.t2d <- GetFC(df.ctrl, df.t2d)

################################################ GET DEGS

GetDegs <- function(df.fc, cutoff) {
  df.degs <- NULL
  df.degs <- df.fc[abs(df.fc$LOG2FC) >= cutoff, ]
  return(df.degs)
}

cutoff <- c(2,2,2,2)

degs.ra   <- GetDegs(fc.ra,   cutoff[1])
degs.mets <- GetDegs(fc.mets, cutoff[2])
degs.cad  <- GetDegs(fc.cad,  cutoff[3])
degs.t2d  <- GetDegs(fc.t2d,  cutoff[4])

f.out <- paste(dir.out, "RA_DEGs_", cutoff[1], ".csv", sep = "")
write.csv(degs.ra, file = f.out, row.names = FALSE)
f.out <- paste(dir.out, "METS_DEGs_", cutoff[2],".csv", sep = "")
write.csv(degs.mets, file = f.out, row.names = FALSE)
f.out <- paste(dir.out, "CAD_DEGs_", cutoff[3], ".csv", sep = "")
write.csv(degs.cad, file = f.out, row.names = FALSE)
f.out <- paste(dir.out, "T2D_DEGs_", cutoff[4], ".csv", sep = "")
write.csv(degs.t2d, file = f.out, row.names = FALSE)
