rm(list=ls())
setwd(dirname(getSourceEditorContext()$path))
dir()
dir.in <- "EXPORT/DEG/"
f.names <- colnames(as.data.frame(read.csv("SERIES.csv", header = TRUE, stringsAsFactors = FALSE)))[2:36]
df.string <- as.data.frame(read.csv("STRING_SIG.csv", header = TRUE, stringsAsFactors = FALSE))
df.F635 <- GetRawVals(f.names, c(1:35))
write.csv(df.F635, file = "EXPORT/F635_ALL.csv", row.names = FALSE)

cutoff <- c(3,2,2,2)

################################################ GET DEGS

list.degs <- list(
  RA   = as.data.frame(read.csv(paste(dir.in, "RA_DEGs_",   cutoff[1], ".csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)),
  METS = as.data.frame(read.csv(paste(dir.in, "METS_DEGs_", cutoff[2], ".csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)),
  CAD  = as.data.frame(read.csv(paste(dir.in, "CAD_DEGs_",  cutoff[3], ".csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)),
  T2D  = as.data.frame(read.csv(paste(dir.in, "T2D_DEGs_",  cutoff[4], ".csv", sep = ""), header = TRUE, stringsAsFactors = FALSE))
)

################################################ MAPPED DEGS

list.mapped.degs <- list(RA   = list.degs[[1]][list.degs[[1]]$Symbol %in% df.string$SYMBOL, ],
                         METS = list.degs[[2]][list.degs[[2]]$Symbol %in% df.string$SYMBOL, ],
                         CAD  = list.degs[[3]][list.degs[[3]]$Symbol %in% df.string$SYMBOL, ],
                         T2D  = list.degs[[4]][list.degs[[4]]$Symbol %in% df.string$SYMBOL, ])


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

PrintStats(list.degs)
PrintStats(list.mapped.degs)

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

list.int <- GetInters(list.degs)
list.mapped.int <- GetInters(list.mapped.degs)

