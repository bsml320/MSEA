MSEA.clust = function(mutations, refseq.length, output ){
	genes_to_test = unique(mutations$RefSeq.ID)
	clust.res = c()	
	for(k in 1:length(genes_to_test)){
		refseq.ID = genes_to_test[[k]]
		exonic = mutations[which(mutations$RefSeq.ID==refseq.ID),]
		gene.length = refseq.length[ refseq.ID ]
		mut_pos = exonic$mut_pos
		
		### sometimes, the annotation has errors
		mut_pos = mut_pos[ mut_pos %in% 1:gene.length ]
		if(length(mut_pos)<4)next
		
		### es.random
		es.random = c()
		for(ii in 1:(gene.length*10)){
			Nm = nrow(exonic)
			mut_pos.pai = sample(1:gene.length, Nm, replace=T)
			
			inc = 1/length(mut_pos.pai)
			table(mut_pos.pai) -> nMut.per.location
			dec = 1/(gene.length-length(nMut.per.location))
			inc.1 = rep(0, gene.length)
			inc.1[as.numeric(names(nMut.per.location))] = inc * nMut.per.location
			dec.1 = rep(-dec, gene.length)
			dec.1[as.numeric(names(nMut.per.location))] = 0
			inc.1 + dec.1 -> ss
			cumsum(ss) -> es.cum
			es.pai = ( max(es.cum) - min(es.cum) )
			
			es.random = c(es.random, es.pai)
		}
		
		### es true
		inc = 1/length(mut_pos)
		table(mut_pos) -> nMut.per.location
		dec = 1/(gene.length-length(nMut.per.location))
		inc.1 = rep(0, gene.length)
		inc.1[as.numeric(names(nMut.per.location))] = inc * nMut.per.location
		dec.1 = rep(-dec, gene.length)
		dec.1[as.numeric(names(nMut.per.location))] = 0
		inc.1 + dec.1 -> ss
		cumsum(ss) -> es.cum
		true.es.cum = es.cum
		es.true = ( max(es.cum) - min(es.cum) )
		
		nes = (es.true-mean(es.random))/sd(es.random); p1 = sum(es.random>=es.true)/(length(es.random) + 1)
		clust.res = rbind(clust.res, c(as.character(exonic$Symbol[1]), refseq.ID, es.true, nes, p1, nrow(exonic), which.min(es.cum), min(es.cum), which.max(es.cum), max(es.cum) ) )
		
		cat(".", sep="")
	}
	colnames(clust.res) = c("symbol", "refseq", "es", "nes", "pvalue", "nMut", "es.min.location", "es.min", "es.max.location", "es.max")
	write.table(cluster.res, file=output, row.names=F, col.names=T, quote=F, sep="\t")
}


MSEA.domain = function(mutations, domain, refseq.length, M1.output=NULL, M2.output=NULL, M3.output=NULL){
	require(MASS)
	genes_to_test = unique(mutations$RefSeq.ID)
	
	M3.res = c()
	all.anno = c()
	M1.res = c()
	M2.res = c()
	
	for(k in 1:length(genes_to_test)){
		refseq.ID = genes_to_test[[k]]
		exonic = mutations[which(mutations$RefSeq.ID==refseq.ID),]
		gene.length = refseq.length[ refseq.ID ]
		mut_pos = exonic$mut_pos
		
		### sometimes, the annotation has errors
		mut_pos = mut_pos[ mut_pos %in% 1:gene.length ]
		if(length(mut_pos)<4)next
		
		### if no domain info or amino acid length does not match, skip
		domain.info = domain[domain$refseq.ID==refseq.ID,]
		domain.info = domain.info[which(gene.length == as.numeric(domain.info[,4])), ]
		if(nrow(domain.info)<1)next
		
		### H0
		table(mut_pos) -> nMut.per.location
		rep(0, gene.length) -> n_mut
		n_mut[ as.numeric(names(nMut.per.location)) ] = nMut.per.location
		glm.nb(n_mut ~ 1) -> h0
		
		########################
		### M1.p
		########################
		domains.X = c()
		this.M1.res = c()
		for(k1 in 1:nrow(domain.info)){
			domain.start  = domain.info$domain.start[k1];
			domain.end    = domain.info$domain.end[k1];
			domain.length = domain.end - domain.start
			if( gene.length - domain.length < 10)next  ### if the domain is so long that it is nearly the same length of the gene itself, not eligible
			
			### keep record of domain regions
			rep(0, gene.length ) -> X
			X[domain.start:domain.end] = 1
			domains.X = rbind(domains.X, X)
			
			if( sum(mut_pos %in% domain.start:domain.end)!=0 ){
				glm.nb(n_mut ~ 1+as.factor(X) ) -> h1
				anova(h0, h1, test="Chisq") -> fit
				M1.res      = rbind(M1.res,      c(domain.info[k1, ], M1.p=fit[2,"Pr(Chi)"]) )
				this.M1.res = rbind(this.M1.res, c(domain.info[k1, ], M1.p=fit[2,"Pr(Chi)"]) )
			}
		}
		if(nrow(this.M1.res)!=0){
			idx = which.min(as.numeric(this.M1.res[,10]))
			M1.p = min(as.numeric(this.M1.res[,10]))
			M1.site = this.M1.res[idx, 9]
		} else {
			M1.p = 1
			M1.site = "NULL"
		}
		if(is.null(this.M1.res))next
		########################
		
		########################
		### M3.p
		########################
		apply(domains.X, 2, function(u)sum(u)!=0) -> domains.X.M3
		if(length(unique(domains.X.M3))==1)next
		
		glm.nb( n_mut ~ 1 + as.factor(domains.X.M3) ) -> h3
		M3.anova = anova(h0, h3, test="Chisq")
		M3.anova[2,"Pr(Chi)"] -> M3.p
		
		########################
		### M2.p
		########################	
		cc = 1; new = c(cc)
		for(nn in 2:length(domains.X.M3) ) {
			if(domains.X.M3[nn] == domains.X.M3[nn-1]){
				new = c(new, cc)
			} else {
				cc = cc+1
				new = c(new, cc)
			}
		}
		new[which(domains.X.M3==0)] = 0
		table(new) -> new.tab
		
		this.M2.res = c()
		M2.anova = list()
		for(k1 in 2:length(new.tab)){
			tag = names(new.tab)[k1]
			domains.X.M2 = ifelse(new==tag, 1, 0)
			
			glm.nb( n_mut ~ 1 + as.factor(domains.X.M2) ) -> h2
			
			anova(h0, h2, test="Chisq") -> fit
			M2.res      = rbind(M2.res,      c(exonic[1,10], gene, paste(min(which(domains.X.M2==1)), " - ", max(which(domains.X.M2==1)), sep=""), fit[2,"Pr(Chi)"] ) )
			this.M2.res = rbind(this.M2.res, c(exonic[1,10], gene, paste(min(which(domains.X.M2==1)), " - ", max(which(domains.X.M2==1)), sep=""), fit[2,"Pr(Chi)"] ) )
		}
		if(nrow(this.M2.res)!=0){
			idx = which.min(as.numeric(this.M2.res[,10]))
			M2.p = min(as.numeric(this.M2.res[,10]))
			M2.site = "M2"
		} else {
			M2.p = 1
			M2.site = "NULL"
		}
		########################
		
		results = c(M1.p, M2.p, M3.p)
		min.site = c(M1.site, "M2", "M3")[which.min(results)]
		M3.res = rbind(M3.res, c(exonic[1,10], gene, results, min.site ) )
	}
	
	write.table(M1.res, file=M1.output, row.names=F, col.names=F, quote=F, sep="\t")
	write.table(M2.res, file=M2.output, row.names=F, col.names=F, quote=F, sep="\t")
	write.table(M3.res, file=M3.output, row.names=F, col.names=F, quote=F, sep="\t");
}
