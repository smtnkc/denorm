library(readr)
library(rstudioapi)
library(data.table)
library(marray)
source("https://bioconductor.org/biocLite.R")
library(limma)
library(marray)

rm(list=ls())
setwd(dirname(getSourceEditorContext()$path))
dir()

df.map <- as.data.frame(read.csv("MAP.csv", header = TRUE, stringsAsFactors = FALSE))
df.gpl <- as.data.frame(read.csv("GPL.csv", header = TRUE, stringsAsFactors = FALSE))

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
  #print(colnames(df))
  #v.symbols <- unlist(df[, by])
  df$"F635.Median" <- round(df$"F635.Median", 0)
  df$"TrueVal" <- round(df$"TrueVal", 2)
  df$"SeriesVal" <- round(df$"SeriesVal", 2)
  #row.names(df) <- v.symbols
  return(df)
}

df.map <- FormatDf(df.map)
df.gpl <- FormatDf(df.gpl)

############################################################

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

############################################################

for(i in 10:35) {
  state <- i
  source("denorm.R")
  cat(length(unmapped), "probes cannot be mapped!\n")
  cat("Wrongs for", colnames(df.series)[2], "=", nrow(df.wrongs[df.wrongs$Symbol != "", ]), "\n")
  cat("Final Wrongs for", colnames(df.series)[2], "=", nrow(df.wrongs.final), "\n")
  fout_final <- paste("EXPORT/FINAL/", colnames(df.series)[2], ".csv", sep = "")
  write.csv(df.mapped.grouped, file = fout_final, row.names = FALSE)
  fout_wrong <- paste("EXPORT/WRONG/", colnames(df.series)[2], "_wrong.csv", sep = "")
  write.csv(df.wrongs.final, file = fout_wrong, row.names = FALSE)
  fout_raw <- paste("EXPORT/RAW/", colnames(df.series)[2], "_raw.csv", sep = "")
  write.csv(df.clean.raw[, c(4:9)], file = fout_raw, row.names = FALSE)
}

ctrl4 <- as.data.frame(read.csv("EXPORT/FINAL/GSM577966_CTRL_04.csv", header = TRUE, stringsAsFactors = FALSE))
mets6 <- as.data.frame(read.csv("EXPORT/FINAL/GSM577983_MetS_06.csv", header = TRUE, stringsAsFactors = FALSE))
