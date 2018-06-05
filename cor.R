list.raw.msdegs <- GetRawMSDegs(df.raw.all, list.msdegs)
list.cor.matrices <- GetCorMatrices(list.raw.msdegs)
list.msdeglinks <- GetMSDegLinks(df.raw.all, list.cor.matrices, comb_score_cutoff)
list.cor.edges <- GetCorEdges(list.msdeglinks, list.cor.matrices)

list.sig.cor.edges <- GetSigCorEdges(list.cor.edges, p_val)
list.sig.cor.genes <- GetSigCorGenes(list.sig.cor.edges)
list.inter.degs <- GetInterDegs(list.sig.cor.genes)