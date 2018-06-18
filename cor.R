list.raw.msdegs <- GetRawMSDegs(df.raw.all, list.msdegs)
list.cor.matrices <- GetCorMatrices(list.raw.msdegs)
list.msdeglinks <- GetMSDegLinks(df.raw.all, list.cor.matrices, comb_score_cutoff)
list.cor.edges <- GetCorEdges(list.msdeglinks, list.cor.matrices)

list.sig.cor.edges <- GetSigCorEdges(list.cor.edges, p_val)
list.sig.cor.nodes <- GetSigCorNodes(list.sig.cor.edges)
list.inter.nodes <- GetInterNodes(list.sig.cor.nodes)
list.inter.edges <- GetInterEdges(list.sig.cor.edges)

ExportEdgeLists(list.sig.cor.edges, p_val, comb_score_cutoff, fc_cutoff)
ExportInterEdges(list.inter.edges, p_val, comb_score_cutoff, fc_cutoff)
