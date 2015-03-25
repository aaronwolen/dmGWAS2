SNP2Gene.match <-
function(assoc.file, snp2gene.file, boundary=20)
{
	cat("start searching at ", format(Sys.time(), "%H:%M, %b %d %Y"), " ...\n", sep="")
	assoc <- read.table(assoc.file, header=T)
	cat("read in ", nrow(assoc), " SNPs in ", assoc.file, "\n", sep="");
	map.data <- read.table(snp2gene.file, header=T)
	cat("read in ", nrow(map.data), " SNPs in ", snp2gene.file, "\n", sep="")
	colns = names(map.data)
	if(is.element("SNP", colns) & is.element("Dist", colns) & is.element("Gene", colns))
	{
		affy.map2 <- map.data[as.numeric(map.data$Dist)<as.numeric(1000*boundary), ]
		merge(affy.map2, assoc, by="SNP") -> gene.map
		cat("finish searching at ", format(Sys.time(), "%H:%M, %b %d %Y"), " ...\n", sep="")
		gene.map = data.frame(Gene=gene.map$Gene, SNP=gene.map$SNP, P=gene.map$P)
		gene.map = gene.map[gene.map$Gene!="---", ]
		
		#write.table(gene.map, file="gene.map.tmp", row.names=F, col.names=F, quote=F, sep="\t")
		#tmp = read.table("gene.map.tmp", sep=" ")
		#tmp = unique(tmp)
		#write.table(tmp, file="gene.map.tmp", row.names=F, col.names=F, quote=F, sep="\t")
		#gene.map = read.table("gene.map.tmp", sep="\t")
		#names(gene.map) = c("Gene", "SNP", "P")
		#file.remove("gene.map.tmp")
		#rm(tmp)
		
		return(gene.map)
	}
	else
	{
		stop("Please make sure that a header line is provided in the snp2gene.file, including \"SNP\", \"Dist\", and \"Gene\"")
	}
}
#==================================================================================================================#

PCombine <-
function(gene.map, method='NULL')
{
	#gene.map: gene, SNP, P
	if(!is.element(method, c('simes', 'smallest', 'fisher', 'gwbon')) | method=='NULL')
	{
		stop('please choose a method: smallest, simes, fisher, gwbon (gene-wise bonferroni)')
	}
	
	tapply(gene.map$P, gene.map$Gene, c) -> t
	idx = lapply(t, function(x){return(length(x)==0)})
	t.clean = t[!unlist(idx)]
	
	if(method=='simes')
	{
		t2 = lapply(t.clean, function(x){tmp=sort(x); pos=seq(1:length(tmp)); tmp2=tmp*length(tmp)/pos; return(min(tmp2))})
		gene2weight = data.frame(gene=names(t2), weight=unlist(t2))
	}
	if(method=='smallest')
	{
		tapply(gene.map$P, gene.map$Gene, min) -> gene2weight
		data.frame(gene=names(gene2weight), weight=gene2weight) -> gene2weight
		gene2weight = gene2weight[!is.na(gene2weight[,2]), ]
	}
	if(method=='fisher')
	{
		tapply(gene.map$P, gene.map$Gene, function(x){sum(-2*log(x))->chi2; p=1-pchisq(chi2, 2*length(x)); return(p);}) -> gene2weight
		data.frame(gene=names(gene2weight), weight=gene2weight) -> gene2weight
		gene2weight = gene2weight[!is.na(gene2weight[,2]), ]
	}
	if(method=='gwbon')
	{
		tapply(gene.map$P, gene.map$Gene, function(x){tmp2 = p.adjust(x, method='bonferroni'); return(min(tmp2));}) -> gene2weight
		data.frame(gene=names(gene2weight), weight=gene2weight) -> gene2weight
		gene2weight = gene2weight[!is.na(gene2weight[,2]), ]
	}
	return(gene2weight)
}
#==================================================================================================================#

statResults <-
function(res.list, top=500)
{
	if (!require(igraph)) {
		stop('igraph must be pre-installed!\n')
	}	
	sets.length <- c()
	sets.cc <- c()
	sets.as <- c()
	sets.ad <- c()
	genes <- c()
	genesets.name = names(res.list$genesets.clear)
	k=1
	while(TRUE)
	{
		match(res.list$zi.ordered[k,1], genesets.name) -> idx
		genes <- c(genes, res.list$genesets.clear[[idx]])
		genes <- unique(genes)
		sets.length <- rbind(sets.length, length(genes));
		#match(genes, V(res.list$GWPI)$name)-1 -> tmp
		induced.subgraph(res.list$GWPI, genes) -> sub.tmp
		sets.cc <- c(sets.cc, transitivity(sub.tmp))
		sets.as <- c(sets.as, average.path.length(sub.tmp))
		deg <- degree(sub.tmp)
		sets.ad <- c(sets.ad, mean(deg))
		k = k+1
		if(k==top)break();
	}
	par(mfrow=c(2,2))
	plot(sets.length, main='No. of genes in total', xlab='Index of modules', ylab='No.of genes')
	plot(sets.cc, main='Clustering coefficient', xlab='Index of modules', ylab='Clustering coefficient')
	plot(sets.as, main='Average shortest path', xlab='Index of modules', ylab='Average shortest path')
	plot(sets.ad, main='Average degree', xlab='Index of modules', ylab='Average degree')	
}
#==================================================================================================================#

simpleChoose <-
function(res.list, top=0.01, plot=FALSE)
{
	if (!require(igraph)) {
		stop('igraph must be pre-installed!\n')
	}
    
      if (top<1 & top>0) 
        top<-round(nrow(res.list$zi.ordered)*top)      
	
	match(res.list$zi.ordered[1:top,1], names(res.list$genesets.clear)) -> mod.idx
	genes=unique(unlist(res.list$genesets.clear[mod.idx]))
	subG <- induced.subgraph(res.list$GWPI, genes)
	results = list(modules=res.list$genesets.clear[mod.idx], subnetwork=subG)
	x = V(subG)$weight
	nx = (x-min(x))/(max(x)-min(x))
	if(plot)
       tkplot(subG, vertex.label=V(subG)$name, vertex.size=5, vertex.color=gray(nx), vertex.label.dist=1, layout=layout.fruchterman.reingold)
	return(results)
}
#==================================================================================================================#

dualEval <-
function(resfile1, resfile2)
{
	load(resfile1)
	dis.list <- res.list
	rm(res.list)
	load(resfile2)
	eval.list <- res.list
	rm(res.list)
	
	dis.genesets <- dis.list$genesets.clear
	dis.length.null.dis <- dis.list$genesets.length.null.dis
	dis.length.null.stat <- dis.list$genesets.length.null.stat
	eval.length.null.dis <- eval.list$genesets.length.null.dis
	eval.length.null.stat <- eval.list$genesets.length.null.stat
	zi.by.dis <- data.frame(gene=names(dis.genesets), zi.dis=-9, za.dis=-9, nominalP.dis=-9, zi.eval=-9, za.eval=-9, nominalP.eval=-9)
	for(k in 1:length(dis.genesets))
	{
		genes <- dis.genesets[[k]]
		match(genes, dis.list$graph.g.weight[,1]) -> idx
		idx <- idx[!is.na(idx)]
		if(length(idx)!=0){zi.by.dis[k,2] <- sum(dis.list$graph.g.weight[idx, 2])/sqrt(length(idx))}
		tmp <- dis.length.null.stat[[as.character(length(idx))]]
		l.zperm <- dis.length.null.dis[[as.character(length(idx))]]
		zi.by.dis[k,3] = (zi.by.dis[k,2]-tmp[1])/tmp[2]
		zi.by.dis[k,4]=sum(l.zperm>=zi.by.dis[k,2])/length(l.zperm)
		
		#check the weight of these genes in evaluation dataset
		match(genes, eval.list$graph.g.weight[,1]) -> idx
		idx <- idx[!is.na(idx)]
		if(length(idx)!=0){zi.by.dis[k,5] <- sum(eval.list$graph.g.weight[idx, 2])/sqrt(length(idx))}
		
		if(!is.element(length(idx), names(eval.length.null.stat)))
		{
			print(k)
			print(dis.genesets[k])
			next()
		}
		tmp <- eval.length.null.stat[[as.character(length(idx))]]
		l.zperm <- eval.length.null.dis[[as.character(length(idx))]]
		zi.by.dis[k,6] = (zi.by.dis[k,5]-tmp[1])/tmp[2]
		zi.by.dis[k,7]=sum(l.zperm>=zi.by.dis[k,5])/length(l.zperm)
		#cat(k, '.', sep='')
	}
	zod = zi.by.dis[order(zi.by.dis[,3], decreasing=T), ]
	zod.sig = zod[zod[,3]>quantile(zod[,3], probs=.95)&zod[,6]>quantile(zod[,6], probs=.95),]
	zod.res = list(zod=zod, zod.sig=zod.sig)
	return(zod.res)
}
#==================================================================================================================#

moduleChoose <-
function(seed.nodes, res.list, plot=FALSE)
{
	if (!require(igraph)) {
		stop('igraph must be pre-installed!\n')
	}	
	match(seed.nodes, names(res.list$genesets.clear)) -> idx
	module.nodes = unique(unlist(res.list$genesets.clear[idx]))
	subG = induced.subgraph(res.list$GWPI, module.nodes)
	results = list(modules=res.list$genesets.clear[idx], subnetwork=subG)
	if(plot)
	{
		x = V(subG)$weight
		nx = (x-min(x))/(max(x)-min(x))
		tkplot(subG, vertex.label=V(subG)$name, vertex.size=5, vertex.color=gray(nx), vertex.label.dist=1, layout=layout.fruchterman.reingold)
	}
	#return(list(subnetwork = subG, module.list = results))
	return(results)
}
#==================================================================================================================#

zn.permutation <- function(module.list, gene2snp, gene2snp.method="smallest", original.file, permutation.dir)
{
	tmp <- read.table(file=original.file, header=T)
	header = names(tmp)
	if((!is.element("SNP", header)) || (!is.element("P", header)))
	{
		stop("please make sure the head line of original.file contains SNP and P...")
	}
	snp2p <- data.frame(tmp$SNP, tmp$P)
	cat("SNPs in total: ", length(snp2p[,1]), sep="")
	
	match(gene2snp[,2], snp2p[,1]) -> x1;
	x2 <- x1[!is.na(x1)]
	gene2snp2p <- data.frame(Gene=gene2snp[!is.na(x1),1], SNP=gene2snp[!is.na(x1),2], P=snp2p[x2,2]); #only SNPs covered by genes
	print(paste("SNPs covered by genes: ", length(gene2snp2p[,1]), sep=""))
	
	gene2weight = PCombine(gene2snp2p, method=gene2snp.method)
	
	perm.files = dir(permutation.dir)
	
	zn.matrix = matrix(-9, ncol=length(perm.files)+1, nrow=length(module.list))
	rownames(zn.matrix)=names(module.list)
	
	
	for(k in 1:length(module.list))
	{
		genes = module.list[[k]]
		match(genes, gene2weight[,1]) -> idx
		idx = idx[!is.na(idx)]
		zn.matrix[k, length(perm.files)+1]=sum(qnorm(1-as.numeric(gene2weight[idx, 2])))/sqrt(length(idx))
	}
	
	for(l in 1:length(perm.files))
	{
		tmp <- read.table(file=paste(permutation.dir, perm.files[l], sep="/"), header=T)
		header = names(tmp)
		if((!is.element("SNP", header)) || (!is.element("P", header)))
		{
			stop(paste("please make sure the head line contains SNP and P in the file ", permutation.dir, "/", perm.files[l], sep=""))
		}
		snp2p <- data.frame(tmp$SNP, tmp$P)
		
		match(gene2snp[,2], snp2p[,1]) -> x1;
		x2 <- x1[!is.na(x1)]
		gene2snp2p <- data.frame(Gene=gene2snp[!is.na(x1),1], SNP=gene2snp[!is.na(x1),2], P=snp2p[x2,2]); #only SNPs covered by genes
	
		gene2weight = PCombine(gene2snp2p, method=gene2snp.method)
		
		for(k in 1:length(module.list))
		{
			genes = module.list[[k]]
			match(genes, gene2weight[,1]) -> idx
			idx = idx[!is.na(idx)]
			zn.matrix[k, l]=sum(qnorm(1-as.numeric(gene2weight[idx, 2])))/sqrt(length(idx))
		}
		cat(l, ".", sep="")
	}
	
	zn.nominal = matrix(-9, ncol=3, nrow=length(module.list))
	colnames(zn.nominal)=c("Seed", "Zm", "EmpiricalP")
	for(k in 1:length(module.list))
	{
		zn.nominal[k,1] = as.character(names(module.list[k]))
		zn.nominal[k,2] = zn.matrix[k, length(perm.files)+1];
		tmp = zn.matrix[k, 1:length(perm.files)];
		zn.nominal[k,3] = sum(tmp>as.numeric(zn.nominal[k,2]))/length(perm.files)
	}
	zn.nominal=as.data.frame(zn.nominal)
	
	return(list(Zn.EmpiricalP=zn.nominal, Zn.Permutation=zn.matrix))
}
#==================================================================================================================#

integGM <-
function (G, genes, weights, simplify=T){
	if (!require(igraph)) {
		stop('igraph must be pre-installed!\n')
	}
	if (is.element("weight", list.vertex.attributes(G))) {
		cat("Warning: previous G node weight replaced!\n")
	}
	names(weights) <- genes
	genes <- intersect(genes,V(G)$name)
	subG <- induced.subgraph(G,genes)
	V(subG)$weight <- weights[V(subG)$name]
	if (simplify) subG <- simplify(subG)
	return(subG)
}
#==================================================================================================================#

node2treePath <-
function (G, Tnodes, node){
	if (!require(igraph)) {
		stop('igraph must be pre-installed!\n')
	}
	tmp.path <- get.all.shortest.paths(G, node, Tnodes)$res
	tmp.l <- unlist(lapply(tmp.path, length))
	index <- which(tmp.l == min(tmp.l))
	
	tmp.path = tmp.path[index]
	tmp.sum <- unlist(lapply(tmp.path, function(x)return(sum(V(G)[x]$weight))))
	index <- which(tmp.sum == max(tmp.sum))
	
	selected.path = tmp.path[index]
	collect <- unlist(lapply(selected.path, function(x)return(V(G)[x]$name)))
	
	return(collect)
}
#==================================================================================================================#

seedQueryJ_2.4 <-
function (G, seed, search_r = 2, r = 0.1){
	if (!require(igraph)) {
		stop('igraph must be pre-installed!\n')
	}
	net <- G
	d <- search_r
	if (!is.element("name", list.vertex.attributes(net))) {
		stop("Graph must have 'name' attribute")
	}
	if (!is.element("weight", list.vertex.attributes(net))) {
		stop("Graph must have 'weight' attribute")
	}
	subG <- induced.subgraph(net, seed)
	if (!is.connected(subG)) {   		#the seed must be connected
		stop("Input seeds are disjoint")
	}
	in.nodes <- V(subG)$name
	while (TRUE) {
		subx <- V(subG)$name
		for (rad in 1:d) {
			subsum <- sum(V(subG)$weight)/sqrt(length(subx)) #Peilin
			
			tmp.neigh <- unlist(neighborhood(net, order = rad, nodes = V(subG)$name)) 
			pot.nodes <- V(net)[tmp.neigh]$name
			pot.nodes <- setdiff(pot.nodes, in.nodes)
			if (length(pot.nodes) == 0) break
			sub.weg <- V(net)[pot.nodes]$weight
			best.nodes <- pot.nodes[which(sub.weg == max(sub.weg))]
			
			subsum.u <- (sum(V(subG)$weight) + V(net)[best.nodes[1]]$weight)/sqrt(length(subx)+1)
			
			if (subsum.u > subsum * (1 + r)) {
				tmp <- unlist(lapply(best.nodes, function(x) node2treePath(net,V(subG)$name, x)))
				in.nodes <- c(tmp, V(subG)$name)
				subG <- induced.subgraph(net, in.nodes)
				break
			}
		}
		if (length(subx) == vcount(subG)) break
	}
	return(subG)
}
#==================================================================================================================#

globalQueryJ_2.4 <-
function (G, search_r = 2, r = 0.2, min.size=5){
	if (!require(igraph)) {
		stop('igraph must be pre-installed!\n')
	}	
	sublist = list()
	for (node in V(G)$name) {
		ng <- seedQueryJ_2.4(G, node, search_r, r)
		if (vcount(ng) >= min.size) sublist[[node]] <- ng       #minimum size filtering
	}
	return(sublist)
}
#==================================================================================================================#

dms_2.4 <-
function(network, gene2weight, d=2, r=0.1)
{
###########################################################################################################
######################################## Starting dense module searching #########################################
####################################################################################################################
#==================================================================================================================#
library(igraph)

if(min(gene2weight[,2])<=0 | max(gene2weight[,2])>=1)
{
	stop("P values out of range 0<p<1");
}
cat("genes used: ", length(gene2weight[,1]), "\n", sep="")

rawG <- graph.data.frame(network,directed=F)
g.weight <- sapply(as.numeric(gene2weight[,2]),function(x) qnorm(1-x))
intG <- integGM(rawG,as.character(gene2weight[,1]),g.weight)
GWPI <- simplify(intG)

cat("start searching at ", format(Sys.time(), "%H:%M, %b %d %Y"), " ...\n", sep="")
dm.result <- globalQueryJ_2.4(GWPI,search_r=d,r)

#================================================ extract genesets ================================================#
cat("extracting modules...\n", sep="")
genesets <- list()
for(k in 1:length(dm.result))
{
	node = names(dm.result[k])
	g = dm.result[[k]]
	genesets[[node]] <- V(g)$name
}
seed.genes <- names(genesets)

#================================== clean genesets by removing identical records ==================================#
cat("removing identical modules...\n", sep="")
identical.idx <- list()
for(k in 1:length(genesets))
{
	tmp.idx <- c(k);
	for(kt in 1:length(genesets))
	{
		if(kt==k)next()
		genesk <- genesets[[k]]
		genest <- genesets[[kt]]
		if(length(genesk)!=length(genest))next()
		overlap = intersect(genesk, genest)
		if(length(overlap)==length(genest))
		{
			tmp.idx <- c(tmp.idx, kt)
		}
	}
	if(length(tmp.idx)>1)
	{
		tmp.idx <- sort(tmp.idx)
		identical.idx[[seed.genes[k]]] <- tmp.idx
		#cat(k, ".", sep="")
	}
}

toremove.idx <- c()
for(k in 1:length(identical.idx))
{
	tmp.idx <- identical.idx[[k]]
	toremove.idx <- c(toremove.idx, tmp.idx[-1])
}
toremove.idx <- unique(toremove.idx)
genesets.clear <- genesets[-toremove.idx]

#================================================= random network =================================================#
cat("permutation on random network...\n", sep="")
genesets.length <- c()
for(k in 1:length(genesets.clear))
{
	genes <- genesets.clear[[k]]
	genesets.length <- c(genesets.length, length(genes))
}
genesets.length <- unique(genesets.length)

genes.idx <- seq(1, length(V(GWPI)$name))
graph.g.weight = data.frame(GWPIene=V(GWPI)$name, gain.weight=V(GWPI)$weight)
genesets.length.null.dis <- list()
length.max = max(genesets.length)+5
for(k in 5:length.max)
{
	l.zperm <- c()
	for(j in 1:100000)
	{
		idx.pseudo=sample(genes.idx, size=k)
		l.zperm <- c(l.zperm, sum(graph.g.weight[idx.pseudo, 2])/sqrt(length(idx.pseudo)));
	}
	genesets.length.null.dis[[as.character(k)]] = l.zperm
	cat(k, ".", sep="");
}

genesets.length.null.stat <- list()
for(k in 5:length.max)
{
	l.zperm <- genesets.length.null.dis[[as.character(k)]]
	k.mean <- mean(l.zperm)
	k.sd <- sd(l.zperm)
	genesets.length.null.stat[[as.character(k)]] = c(k.mean, k.sd)
}

################################################### Normalization ##################################################
#==================================================================================================================#
zim <- data.frame(gene=names(genesets.clear), Zm=-9, Zn=-9, zcount=-9)
for(k in 1:length(genesets.clear))
{
	genes <- genesets.clear[[k]]
	match(genes, graph.g.weight[,1]) -> idx
	idx <- idx[!is.na(idx)]
	zim[k,2] <- sum(graph.g.weight[idx, 2])/sqrt(length(idx)) 
	
	tmp <- genesets.length.null.stat[[as.character(length(idx))]]
	zim[k, 3]=(zim[k,2]-tmp[1])/tmp[2]
	zim[k, 4]=sum(genesets.length.null.dis[[as.character(length(idx))]]>=zim[k,2]) 
}
zom = zim[order(zim[,3], decreasing=T), ]
#==================================================================================================================#

#================================================== save results ==================================================#
res.list <- list()
res.list[["GWPI"]]                        = GWPI
res.list[["graph.g.weight"]] 		          = graph.g.weight
res.list[["genesets.clear"]] 		          = genesets.clear
res.list[["genesets.length.null.dis"]] 	  = genesets.length.null.dis
res.list[["genesets.length.null.stat"]] 	= genesets.length.null.stat
res.list[["zi.matrix"]]                   = zim
res.list[["zi.ordered"]]                  = zom
save(res.list, file="RESULT.list.RData")
cat("finished at ", format(Sys.time(), "%H:%M, %b %d %Y"), " ...\n", sep="")
return(res.list);
########################################### End of dense module searching ###########################################
}

