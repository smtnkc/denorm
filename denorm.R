DenormRunner <- function(df.gpl.mapped, range) {
  for(i in range) {
    state <- i
    df.series <- as.data.frame(read.csv("SERIES.csv", header = TRUE, stringsAsFactors = FALSE))[, c(1, 1 + state)]
    f.name <- paste(dir.csv, colnames(df.series)[2], ".csv", sep = "")
    cat("Running for", f.name, "\n")
    df.raw <- as.data.frame(read.table(f.name, sep = ";", header = TRUE, stringsAsFactors = FALSE))
    cat("Raw file size =", nrow(df.raw), "x", ncol(df.raw), "\n")
    colnames(df.raw)[4] <- "Probe_Name"
    df.raw <- FormatDf(df.raw)
    total_intensity <- sum(df.raw$F635.Median)
    cat("Total intensity = ", total_intensity, "\n")

    ############################################################
    
    df.clean.raw <- df.raw[df.raw$Probe_Name != "", ]
    df.clean.raw <- df.clean.raw[order(df.clean.raw$Probe_Name), ]
    df.clean.gpl <- df.gpl.mapped[order(df.gpl.mapped$Probe_Name), ]
    cat("CleanRaw file size =", nrow(df.clean.raw), "x", ncol(df.clean.raw), "\n")
    cat("CleanGpl file size =", nrow(df.clean.gpl), "x", ncol(df.clean.gpl), "\n")
   
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
    cat(length(unmapped), "probes cannot be mapped!")
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

    df.wrongs <- df.clean.raw[(df.clean.raw$TrueVal - df.clean.raw$SeriesVal) > 0.01, ]
    
    ############################################################

    df.mapped <- df.clean.raw[df.clean.raw$Symbol != "", ]
    df.mapped <- df.mapped[, c(6,9,10,11)]
    df.mapped <- df.mapped[, c(2,1,3,4)]

    df.mapped.grouped <- GroupByCol(df.mapped, by = 1, method = median) # Mean, Median, or Max
    df.wrongs.final <- df.mapped.grouped[(df.mapped.grouped$TrueVal - df.mapped.grouped$SeriesVal) > 0.01, ]

    ############################################################

    cat(length(unmapped), "probes cannot be mapped!\n")
    cat("Wrongs for", colnames(df.series)[2], "=", nrow(df.wrongs[df.wrongs$Symbol != "", ]), "\n")
    cat("Wrongs for grouped", colnames(df.series)[2], "=", nrow(df.wrongs.final), "\n")
    fout_final <- paste("EXPORT/FINAL/", colnames(df.series)[2], ".csv", sep = "")
    write.csv(df.mapped.grouped, file = fout_final, row.names = FALSE)
    fout_wrong <- paste("EXPORT/WRONG/", colnames(df.series)[2], "_wrong.csv", sep = "")
    write.csv(df.wrongs.final, file = fout_wrong, row.names = FALSE)
    fout_raw <- paste("EXPORT/RAW/", colnames(df.series)[2], "_raw.csv", sep = "")
    write.csv(df.clean.raw[, c(4:9)], file = fout_raw, row.names = FALSE)
  }
}
