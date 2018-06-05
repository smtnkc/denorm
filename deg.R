p_val <- 0.05
fc_cutoff <- 1.5
comb_score_cutoff <- 500

list.degs <- GetDegs(fc_cutoff)

################################################ GET MAPPED SIG DEGS (MSDS) AS A LIST

list.msgenes <- GetMSGenes(df.raw.all, comb_score_cutoff)

list.msdegs <- list(   
  RA   = list.degs[[1]][list.degs[[1]]$Symbol %in% list.msgenes, ],
  METS = list.degs[[2]][list.degs[[2]]$Symbol %in% list.msgenes, ],
  CAD  = list.degs[[3]][list.degs[[3]]$Symbol %in% list.msgenes, ],
  T2D  = list.degs[[4]][list.degs[[4]]$Symbol %in% list.msgenes, ]
)

################################################ STATS

# PrintDegStats(list.degs, fc_cutoff)
# PrintDegStats(list.msdegs, fc_cutoff)
