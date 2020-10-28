# enrichment
enrichment <- function(x, odds.ratio=FALSE){
	pvals = oddratio = matrix(ncol=ncol(x), nrow=nrow(x))
	if (nrow(x) == 2 & ncol(x) > 2){
		for (i in 1:ncol(x)){
			x1 = cbind(x[,i],rowSums(x[,-i]))
			fet = fisher.test(x1)
			pvals[1,i] = fet$p.value
			oddratio[1,i] = fet$estimate
			fet = fisher.test(x1[c(2,1),])
			pvals[2,i] = fet$p.value
			oddratio[2,i] = fet$estimate
		}
	}
	if (nrow(x) > 2 & ncol(x) == 2){
		for (i in 1:nrow(x)){
			x1 = cbind(x[i,],colSums(x[-i,]))
			fet = fisher.test(x1)
			pvals[i,1] = fet$p.value
			oddratio[i,1] = fet$estimate
			fet = fisher.test(x1[c(2,1),])
			pvals[i,2] = fet$p.value
			oddratio[i,2] = fet$estimate
		}
	}
	if (nrow(x) > 2 & ncol(x) > 2){
		for (i in 1:nrow(x)){
			for (j in 1:ncol(x)){
				x1 = rbind(x[i,], colSums(x[-i,]))
				x1 = cbind(x1[,j], rowSums(x1[,-j]))
				fet = fisher.test(x1)
				pvals[i,j] = fet$p.value
				oddratio[i,j] = fet$estimate
			}
		}
	}
	if (odds.ratio)
		return(oddratio)
	else
		return(pvals)
}
