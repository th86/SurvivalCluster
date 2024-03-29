# Fisher's Exact Test Compututation
# 
#
# Author: Tai-Hsien Ou Yang
###############################################################################

load("allattractome.new.rda")





getGeneSymbols = function(innames){
	outnames = sapply(innames, function(x){
	
		if(regexpr("\\?", x) > 0){
			o = strsplit(x, "\\|")[[1]][2]
		}else{
			o = strsplit(x, "\\|")[[1]][1]
		}
		return (o)

	}
	)
}


fisher.test<-function(k,n,K,N){
	#N = population, K = signature set size, n = significant sample size, k = sample enriched in signature set 
	p = K/N		#uniform distribution
	# calculate p-values
	pvals = pbinom(k-1, prob=p, size=n, lower.tail=F)
	return(pvals)
}



enrichmentCompute<-function(  signature.List  ,attractome , probeSize=20530 ){ 

	significant.List=names(table(unlist(signature.List)))
	
	pval=rep(NA,length(attractome$rnaseq) )
	names(pval)= names(attractome$rnaseq) 

	for(i in 1:length(attractome$rnaseq) ){
		k=length(which( getGeneSymbols(attractome$rnaseq[[i]]) %in% significant.List ))
		print(k)
		n=length(significant.List)
		K=length(attractome$rnaseq[[i]])
		N=probeSize
		pval[i]=fisher.test( k, n, K, N)
	}
	

	print(n)
	
	
	return(pval)
}


load("panCan12v47.probenum.rda")

pvalList=matrix(NA,12,10)


#colnames(pvalList) = names(enrichmentCompute(cv.signature.member, allattractome, probeSize=N[1]   )     )
#rownames(pvalList) = names(N)


for(i in 1:12){
	load(  paste("./tcga_result/", names(N)[i], ".cv.signature.rda",sep="")     )
	pvalList[i,]=enrichmentCompute(cv.signature.member, allattractome, probeSize=N[i]   ) 
	
	}
