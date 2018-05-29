rm(list=ls())
setwd(dirname(getSourceEditorContext()$path))
dir()

################################################ GET DEGS

deg.ra <-  as.data.frame(read.csv("EXPORT/DEG/RA_DEGs_3.csv", header = TRUE, stringsAsFactors = FALSE))
deg.mets <-  as.data.frame(read.csv("EXPORT/DEG/METS_DEGs_2.csv", header = TRUE, stringsAsFactors = FALSE))
deg.cad <-  as.data.frame(read.csv("EXPORT/DEG/CAD_DEGs_2.csv", header = TRUE, stringsAsFactors = FALSE))
deg.t2d <-  as.data.frame(read.csv("EXPORT/DEG/T2D_DEGs_2.csv", header = TRUE, stringsAsFactors = FALSE))

################################################ INTERSECTIONS

common_all <- Reduce(intersect, list(deg.ra$Symbol, deg.mets$Symbol, deg.cad$Symbol, deg.t2d$Symbol))
common_cad_ra <- Reduce(intersect, list(deg.cad$Symbol, deg.ra$Symbol))
common_cad_ra_t2d <- Reduce(intersect, list(deg.cad$Symbol, deg.ra$Symbol, deg.t2d$Symbol))
common_cad_ra_mets <- Reduce(intersect, list(deg.cad$Symbol, deg.ra$Symbol, deg.mets$Symbol))
common_cad_t2d_mets <- Reduce(intersect, list(deg.cad$Symbol, deg.t2d$Symbol, deg.mets$Symbol))
common_ra_mets <- Reduce(intersect, list(deg.mets$Symbol, deg.ra$Symbol))
common_ra_t2d_mets <- Reduce(intersect, list(deg.ra$Symbol, deg.t2d$Symbol, deg.mets$Symbol))
common_ra_t2d <- Reduce(intersect, list(deg.ra$Symbol, deg.t2d$Symbol))
common_mets_t2d <- Reduce(intersect, list(deg.t2d$Symbol, deg.mets$Symbol))
