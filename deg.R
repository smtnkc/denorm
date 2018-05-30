source("helper.R")

################################################ MEDIANS

med.ctrl <- GetMedian(GetRawVals(f.names, c(1:9)))
med.ra   <- GetMedian(GetRawVals(f.names, c(10:15)))
med.mets <- GetMedian(GetRawVals(f.names, c(16:21)))
med.cad  <- GetMedian(GetRawVals(f.names, c(22:27)))
med.t2d  <- GetMedian(GetRawVals(f.names, c(28:35)))

################################################ FC

fc.ra   <- GetFC(med.ctrl, med.ra)
fc.mets <- GetFC(med.ctrl, med.mets)
fc.cad  <- GetFC(med.ctrl, med.cad)
fc.t2d  <- GetFC(med.ctrl, med.t2d)

################################################ GET DEGS

cutoff <- c(2,2,2,2)

degs.ra   <- GetDegs(fc.ra,   cutoff[1])
degs.mets <- GetDegs(fc.mets, cutoff[2])
degs.cad  <- GetDegs(fc.cad,  cutoff[3])
degs.t2d  <- GetDegs(fc.t2d,  cutoff[4])

f.out <- paste(dir.deg, "RA_DEGs_", cutoff[1], ".csv", sep = "")
write.csv(degs.ra, file = f.out, row.names = FALSE)
f.out <- paste(dir.deg, "METS_DEGs_", cutoff[2],".csv", sep = "")
write.csv(degs.mets, file = f.out, row.names = FALSE)
f.out <- paste(dir.deg, "CAD_DEGs_", cutoff[3], ".csv", sep = "")
write.csv(degs.cad, file = f.out, row.names = FALSE)
f.out <- paste(dir.deg, "T2D_DEGs_", cutoff[4], ".csv", sep = "")
write.csv(degs.t2d, file = f.out, row.names = FALSE)
