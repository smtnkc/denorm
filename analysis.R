source("helper.R")

df.string <- as.data.frame(read.csv("STRING_SIG.csv", header = TRUE, stringsAsFactors = FALSE))
df.F635 <- GetRawVals(f.names, c(1:35))
write.csv(df.F635, file = "EXPORT/F635_ALL.csv", row.names = FALSE)

cutoff <- c(3,2,2,2)

################################################ GET DEGS

list.degs <- list(
  RA   = as.data.frame(read.csv(paste(dir.deg, "RA_DEGs_",   cutoff[1], ".csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)),
  METS = as.data.frame(read.csv(paste(dir.deg, "METS_DEGs_", cutoff[2], ".csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)),
  CAD  = as.data.frame(read.csv(paste(dir.deg, "CAD_DEGs_",  cutoff[3], ".csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)),
  T2D  = as.data.frame(read.csv(paste(dir.deg, "T2D_DEGs_",  cutoff[4], ".csv", sep = ""), header = TRUE, stringsAsFactors = FALSE))
)

################################################ MAPPED DEGS

list.mapped.degs <- list(RA   = list.degs[[1]][list.degs[[1]]$Symbol %in% df.string$SYMBOL, ],
                         METS = list.degs[[2]][list.degs[[2]]$Symbol %in% df.string$SYMBOL, ],
                         CAD  = list.degs[[3]][list.degs[[3]]$Symbol %in% df.string$SYMBOL, ],
                         T2D  = list.degs[[4]][list.degs[[4]]$Symbol %in% df.string$SYMBOL, ])

PrintStats(list.degs)
PrintStats(list.mapped.degs)

list.int <- GetInters(list.degs)
list.mapped.int <- GetInters(list.mapped.degs)

