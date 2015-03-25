
dms <-
function(network, geneweight, expr1, expr2=NULL, d=1, r=0.1, lambda="default")
{
###########################################################################################################
######################################## Starting dense module searching #########################################
####################################################################################################################
#==================================================================================================================#
library(igraph)

if(is.null(expr1) & is.null(expr2))                   ### compatible with previous versions
{
 res_no_edge<-dms_2.4(network,geneweight,d,r=0.1)
 return (res_no_edge)
}

GWPI <- generate_graph(expr1,expr2,network,geneweight)

if(is.null(expr2))                                    ### release memory         
{rm(expr1); rm(geneweight); rm(network); gc()} else
{rm(expr1); rm(expr2); rm(geneweight); rm(network); gc()}        

if(lambda=="default")
{ 
 lambda_default<-T
 cat("Estimating lambda...\n")
 lambda<-estimate_lambda(GWPI)
} else lambda_default<-F

cat("start searching at ", format(Sys.time(), "%H:%M, %b %d %Y"), " ...\n", sep="")

d<-1        ##order of neighbors
dm.result <- globalQueryJ(GWPI,search_r=d,r,lambda)
if(length(dm.result)==0) stop("No dense modules identified!\n")

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
genesets.clear<-genesets[-toremove.idx]
dm.result <- dm.result[-toremove.idx]

#================================================= random network =================================================#
cat("permutation on random network...\n", sep="")

genesets.length <- unique(as.numeric(lapply(dm.result,vcount)))
genesets.length.null.dis <- list()
cat("module size: ")
for(k in 5:max(genesets.length))
{ 
 cat(k, ".", sep="")
 genesets.length.null.dis[[as.character(k)]] <- random_network(size=k,G=GWPI,lambda)
}

genesets.length.null.stat <- list()
for(k in 5:max(genesets.length))
{
	l.zperm <- genesets.length.null.dis[[as.character(k)]]
	k.mean <- mean(l.zperm)
	k.sd <- sd(l.zperm)
	genesets.length.null.stat[[as.character(k)]] = c(k.mean, k.sd)
}

################################################### Normalization ##################################################
#==================================================================================================================#

calculate_score<-function(g)
{
 if( ecount(g)>0 )                                                        
 subsum <- (1-lambda)*sum(V(g)$weight)/sqrt(vcount(g)) + lambda*sum(E(g)$weight)/sqrt(ecount(g)) else
 subsum <- (1-lambda)*sum(V(g)$weight)/sqrt(vcount(g))
 subsum
}

ms <- data.frame(gene=names(genesets.clear), Sm=-9, Sn=-9)
for(k in 1:length(genesets.clear))          
{
	ms[k,2] <- calculate_score(dm.result[[k]])
     	tmp <- genesets.length.null.stat[[as.character(vcount(dm.result[[k]]))]]
	ms[k, 3]=(ms[k,2]-tmp[1])/tmp[2]
	#zim[k, 4]=sum(genesets.length.null.dis[[as.character(vcount(dm.result[[k]]))]]>=zim[k,2]) 
}
ms_ordered = ms[order(ms[,3], decreasing=T), ]
#==================================================================================================================#

#================================================== save results ==================================================#
res.list <- list()
res.list[["GWPI"]]                                = GWPI
res.list[["genesets.clear"]] 		              = genesets.clear
res.list[["genesets.length.null.dis"]] 	        = genesets.length.null.dis
res.list[["genesets.length.null.stat"]] 	        = genesets.length.null.stat
res.list[["module.score.matrix"]]                 = ms
res.list[["ordered.module.score.matrix"]]         = ms_ordered

if(lambda_default)
save(res.list, file=paste("Lambda_",lambda,"_estimated_by_default_result.RData",sep=""))
if(!lambda_default)
save(res.list, file=paste("Lambda_",lambda,"_specified_by_user_result.RData",sep=""))

cat("finished at ", format(Sys.time(), "%H:%M, %b %d %Y"), " ...\n", sep="")
return(res.list)
########################################### End of dense module searching ###########################################
}

