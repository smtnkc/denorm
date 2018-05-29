df.series <- as.data.frame(read.csv("SERIES.csv", header = TRUE, stringsAsFactors = FALSE))[, c(1, 1 + state)]
fdir <- paste("CSV/", colnames(df.series)[2], ".csv", sep = "")
df.raw <- as.data.frame(read.table(fdir, sep = ";", header = TRUE, stringsAsFactors = FALSE))
colnames(df.raw)[4] <- "Probe_Name"
df.raw <- FormatDf(df.raw)
total_intensity <- sum(df.raw$F635.Median)

############################################################

df.clean.raw <- df.raw[df.raw$Probe_Name != "", ]
#df.clean.raw <- df.clean.raw[order(df.clean.raw$ID, decreasing = TRUE), ]
df.clean.raw <- df.clean.raw[order(df.clean.raw$Probe_Name), ]
df.clean.gpl <- df.clean.gpl[order(df.clean.gpl$Probe_Name), ]

row.names(df.clean.raw) <- NULL
row.names(df.clean.gpl) <- NULL

df.clean.raw[, "Name"] <- ""
df.clean.raw[, "ID_REF"] <- ""
df.clean.raw[, "Symbol"] <- ""

unmapped <- c()
for(i in 1:nrow(df.clean.raw)) {
  if(df.clean.raw[i, "Probe_Name"] == df.clean.gpl[i, "Probe_Name"]) {
    df.clean.raw[i, "Name"] <- df.clean.gpl[i, "Name"]
    df.clean.raw[i, "ID_REF"] <- df.clean.gpl[i, "ID"]
    df.clean.raw[i, "Symbol"] <- df.clean.gpl[i, "Symbol.v12"]
  }
  else
    unmapped <- cbind(unmapped, i)
}
df.clean.raw$ID_REF <- as.integer(df.clean.raw$ID_REF)

############################################################

df.clean.raw[, "TrueVal"] <- ""
df.clean.raw[, "SeriesVal"] <- ""

for(i in 1:nrow(df.clean.raw)) {
  df.clean.raw[i, "TrueVal"] <- round(df.clean.raw[i, "F635.Median"]*10000/total_intensity, 2)
  df.clean.raw[i, "SeriesVal"] <- df.series[df.series$ID_REF == df.clean.raw[i, "ID_REF"], 2]
}

df.clean.raw$TrueVal <- as.numeric(df.clean.raw$TrueVal)
df.clean.raw$SeriesVal <- as.numeric(df.clean.raw$SeriesVal)

############################################################

df.wrongs <- df.clean.raw[df.clean.raw$TrueVal != df.clean.raw$SeriesVal, ]

############################################################

df.mapped <- df.clean.raw[df.clean.raw$Symbol != "", ]
df.mapped <- df.mapped[, c(6,9,10,11)]
df.mapped <- df.mapped[, c(2,1,3,4)]

df.mapped.grouped <- GroupByCol(df.mapped, by = 1, method = median) # Take average or max?
df.wrongs.final <- df.mapped.grouped[df.mapped.grouped$TrueVal != df.mapped.grouped$SeriesVal, ]
