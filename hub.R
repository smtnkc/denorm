list.nodes <- list(  
  RA = as.data.frame(read.csv(file = "EXPORT/EDGELIST/cyto/RA.csv", stringsAsFactors = FALSE), header = TRUE),
  METS = as.data.frame(read.csv(file = "EXPORT/EDGELIST/cyto/METS.csv", stringsAsFactors = FALSE), header = TRUE),
  CAD = as.data.frame(read.csv(file = "EXPORT/EDGELIST/cyto/CAD.csv", stringsAsFactors = FALSE), header = TRUE),
  T2D = as.data.frame(read.csv(file = "EXPORT/EDGELIST/cyto/T2D.csv", stringsAsFactors = FALSE), header = TRUE)
)

list.hubs <- list(
  RA = head(list.nodes$RA[order(-list.nodes$RA$Degree), "name"], 50),
  METS = head(list.nodes$METS[order(-list.nodes$METS$Degree), "name"], 50),
  CAD = head(list.nodes$CAD[order(-list.nodes$CAD$Degree), "name"], 50),
  T2D = head(list.nodes$T2D[order(-list.nodes$T2D$Degree), "name"], 50)
  )

inter.hubs <- GetInterNodes(list.hubs)

list.hubs.ext <- list(
  RA = head(list.nodes$RA[order(-list.nodes$RA$Degree), c("name", "Degree")], 50),
  METS = head(list.nodes$METS[order(-list.nodes$METS$Degree), c("name", "Degree")], 50),
  CAD = head(list.nodes$CAD[order(-list.nodes$CAD$Degree), c("name", "Degree")], 50),
  T2D = head(list.nodes$T2D[order(-list.nodes$T2D$Degree), c("name", "Degree")], 50)
)

ExportHubs(list.hubs.ext)
