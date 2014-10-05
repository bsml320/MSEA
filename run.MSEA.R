source("MSEA-clust.R")

###############################################################################

data = unfactor(read.delim("all_non_silent.SNVs.txt"))
sample.info = unfactor(read.delim("/path/to/Supplementary_Table_3a.txt"))

load("/path/to/human.protein.faa.ll.RData")
#load("/scratch/jiap/data/NCBI/refseq/human.protein.faa.ll.RData")

###############################################################################

### use blca as an example
cancer.samples = as.character(sample.info[sample.info$Cancer.Type=="blca",1])
match(data[,42], cancer.samples) -> idx
nonsilent = data[!is.na(idx), c(1:9, 42)]

###############################################################################

mutations <- do.call(rbind, apply(nonsilent, 1, function(u){ 
	strsplit(u[9], split=",")[[1]] -> v
	do.call(rbind, lapply(v, function(w){
		c(u[-9], unlist(strsplit(w, split=":")), col.names=NULL)
	}))
	}))
colnames(mutations)[9:14] = c("TCGA.ID", "Symbol", "RefSeq.ID", "Exon.index", "DNA_change", "amino_acid_change")
mutations = unfactor(as.data.frame(mutations))
mutations = mutations[which(mutations$DNA_change!=""), ]  ### remove occasional errors
mutations = mutations[!is.na(match(mutations$RefSeq.ID, names(refseq.length))), ]  ### remove refseq genes with no gene length

print(paste("# samples: ", length(unique(mutations$TCGA.ID)), sep=""))
print(paste("# mutations: ", nrow(mutations), sep=""))

###################################################################################################
###################################################################################################

tapply(mutations$amino_acid_change, mutations$RefSeq.ID, length) -> nMut.per.gene
names( which(nMut.per.gene>=4) ) -> genes_to_test
print(paste("# input genes: ", length(nMut.per.gene), "; eligible genes to test: ", length(genes_to_test), sep=""))
match(mutations$RefSeq.ID, genes_to_test) -> idx
new.mutations = unfactor(mutations[!is.na(idx), ])
	
### extract mutation positions
as.numeric(unlist( lapply(new.mutations$amino_acid_change, function(u){ 
		if(grepl("_", u)){
			substr(u, 3, regexec("_", u)[[1]][1]-1  ) -> u1
			regexec("[0-9]+", u1) -> m
			regmatches(u1,m) -> a
		} else {
			regexec("[0-9]+", u) -> m
			regmatches(u,m) -> a
		} } ) )) -> mut_pos
new.mutations$mut_pos = mut_pos
	
###################################################################################################
###################################################################################################

MSEA.clust(new.mutations, refseq.length, "output.txt")

###################################################################################################
###################################################################################################

domain = unfactor(read.delim("protein.gbk.regions.symbol.txt"))
MSEA.domain = (new.mutations, domain, refseq.length, M1.output="M1.output.txt", M2.output="M2.output.txt", M3.output="M3.output.txt")
