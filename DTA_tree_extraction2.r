DTA_tree_extraction2 = function(dta_tre, mostRecentSamplingDatum)
	{
		dta_tab = data.frame(matrix(nrow=dim(dta_tre$edge)[1], ncol=5))
		colnames(dta_tab) = c("node1","node2","length","startLoc","endLoc")
		dta_tab[,c("node1","node2")] = as.numeric(dta_tre$edge)
		dta_tab[,c("length")] = as.numeric(dta_tre$edge.length)
		for (i in 1:length(dta_tre$annotations))
			{
				annotations = dta_tre$annotations[[i]]
				if (!is.null(dta_tre$annotations[[i]]))
					{
						dta_tab[i,"endLoc"] = annotations$location
					}
			}
		for (i in 1:length(dta_tre$annotations))
			{
				index = which(dta_tab[,"node2"] == dta_tab[i,"node1"])
				if (length(index) > 0)
					{
						dta_tab[i,"startLoc"] = dta_tab[index,"endLoc"]
					}	else	{
						annotations = dta_tre$root.annotation
						dta_tab[i,"startLoc"] = annotations$location
					}
			}
		treeHeight = max(nodeHeights(dta_tre)[,2], na.rm=T)
		nodeYears = matrix(nrow=dim(dta_tab)[1], ncol=2)
		colnames(nodeYears) = c("startYear","endYear")
		for (i in 1:dim(nodeYears)[1])
			{
				nodeHeight1 = nodeheight(dta_tre, dta_tre$edge[i,1])
				nodeHeight2 = nodeheight(dta_tre, dta_tre$edge[i,2])
				nodeYears[i,1] = mostRecentSamplingDatum-(treeHeight-nodeHeight1)
				nodeYears[i,2] = mostRecentSamplingDatum-(treeHeight-nodeHeight2)
			}
		dta_tab = cbind(dta_tab, nodeYears)
		dta_tab = dta_tab[order(dta_tab[,"startYear"],decreasing=F),]
		dta_tab1 = dta_tab[1,]; dta_tab2 = dta_tab[2:dim(dta_tab)[1],]
		dta_tab2 = dta_tab2[order(dta_tab2[,"endYear"],decreasing=F),]
		return(dta_tab)
	}
