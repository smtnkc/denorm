source("helper.R")

df.map <- as.data.frame(read.csv("MAP.csv", header = TRUE, stringsAsFactors = FALSE))
df.gpl <- as.data.frame(read.csv("GPL.csv", header = TRUE, stringsAsFactors = FALSE))

df.map <- FormatDf(df.map)
df.gpl <- FormatDf(df.gpl)

df.gpl.mapped <- MapProbes(df.gpl, df.map)

DenormRunner(df.gpl.mapped, c(1:1))
