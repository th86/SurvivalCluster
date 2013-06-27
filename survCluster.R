# FCMS Implementation
# 
#
# Author: Tai-Hsien Ou Yang
###############################################################################




cvList<-function( population  ,fd=10){
	cv.mat <- matrix(0,nrow=population,ncol=fd)
	
	cv=list()
	for(i in 1:fd){
		cv.mat[,i]= sample(1:population ,replace=FALSE)
		
		cv[[i]]=list(train = cv.mat[ 1: floor( population*(fd-1)/fd )  , i ],
		test  = cv.mat[  ceiling(( population*(fd-1)/fd)) : population , i ]
		)
	}



	return(cv) 
}



#rm(list=ls())
require("survival")



.notrun<-function(){
cat("Load data\n")
e.exp<-read.delim("KIRC.han.exp")
rownames(e.exp)=e.exp[,1]
e.exp=e.exp[,2:ncol(e.exp)]

surv.t<-read.table("KIRC.han.time", col.names=FALSE)
surv.e<-read.table("KIRC.han.event", col.names=FALSE)


e.surv=cbind(surv.t,surv.e)
rownames(e.surv)=colnames(e.exp)
rm(list=c("surv.t", "surv.e"))
}


#survCluster<-function( e, surv, sizelim=c(10,100), corlim=0.6, pvlim=0.05, cvfold=10 ){

	sizelim=c(10,100) 
	corlim=0.6
	pvlim=5e-4
	cvfold=10

	e=e.exp
	survData= Surv(e.surv[,1], e.surv[,2])

.notrun<-function(){
	e.coeff = rep(NA, nrow(e))
	e.pval  = rep(NA, nrow(e))
	names(e.coeff)=rownames(e)
	names(e.pval)=rownames(e)

	cat("Compute single gene coefficients\n")
	
	prgbar = txtProgressBar(style = 3)
	for( i in 1:nrow(e) )
	{	
		e.vec<-cbind(t(e[i,]), e.surv[,1], e.surv[,2] ) 
		colnames(e.vec)=c("signature","time","status")
		e.vec<-data.frame( e.vec )
		
		cox.summary=summary(coxph(Surv(e.vec[,"time"], e.vec[,"status"]) ~ signature, e.vec ))	
		
		e.coeff[i] = cox.summary$coef[1]    
		e.pval[i]  = cox.summary$logtest[3]
		
		        if (i%%100 == 0)
            			setTxtProgressBar(prgbar, i/nrow(e))
        
	}
	cat("\nDONE\n"   )
}


	
	cv.signature.member=list()
	
	cv<-cvList( ncol(e) ,cvfold  )

	#Cross validation runs
	for(cv.itr in 1:cvfold){ 

	cat("Cross validtion run" 	,cv.itr ,"\n")
	cat("Split probes by survival...")
	
	e.pos = e[intersect( which(e.coeff>0 ), which(e.pval<pvlim ) ) ,]
	e.neg = e[intersect( which(e.coeff<0 ), which(e.pval<pvlim ) ) ,]

	cat("DONE\n"   )

	#Transpose to build probe distance, not samples
	cat("Build distance matrix...")
	e.pos.dist =  1 - cor(  t(e.pos[, cv[[cv.itr]]$train ] )  )   #dist() is for eucledian
	e.neg.dist =  1 - cor(  t(e.neg[, cv[[cv.itr]]$train ] )  )  
	#e.neg.dist = dist(e.neg,upper = TRUE,diag = TRUE) 
	cat("DONE\n"   )

	cat("Cluster...")
	e.pos.tree <- hclust( as.dist(  e.pos.dist ), method = "complete")
 	e.neg.tree <- hclust( as.dist(  e.neg.dist ), method = "complete")
 
	#Filter by the correlation limit
	e.pos.tree.cut <- cutree( e.pos.tree , h=corlim ) #k = length(which(e.coeff>0))/sizelim[2] 
	e.neg.tree.cut <- cutree( e.neg.tree , h=corlim )
	cat("DONE\n"   )

	#Filter by size 
	e.pos.candidate.clusters=intersect(which( table(e.pos.tree.cut)> sizelim[1]) , which( table(e.pos.tree.cut)< sizelim[2]  ) )
	e.neg.candidate.clusters=intersect(which( table(e.neg.tree.cut)> sizelim[1]) , which( table(e.neg.tree.cut)< sizelim[2]  ) )

	#Calculate signatures
	cat("Compute signatures...")
	e.pos.signature = matrix(NA, length(e.pos.candidate.clusters), length(cv[[cv.itr]]$train) )
	e.pos.signature.names = rep(NA, length(e.pos.candidate.clusters))
	e.signature.member=list()
	for(k in 1:length(e.pos.candidate.clusters) ){
		cluster.member = names(e.pos.tree.cut[which(e.pos.tree.cut==e.pos.candidate.clusters[k])])
		e.pos.signature[k , ] = colSums(e[cluster.member, cv[[cv.itr]]$train ])/length( cluster.member )
		e.pos.signature.names[k] = cluster.member[1]
		e.signature.member[[ cluster.member[1] ]]= cluster.member
	}
		rownames(e.pos.signature)=e.pos.signature.names

	e.neg.signature = matrix(NA, length(e.neg.candidate.clusters), length(cv[[cv.itr]]$train) )
	e.neg.signature.names = rep(NA, length(e.neg.candidate.clusters))
	
	for(k in 1:length(e.neg.candidate.clusters) ){
		cluster.member = names(e.neg.tree.cut[which(e.neg.tree.cut==e.neg.candidate.clusters[k])])
		e.neg.signature[k , ] = colSums(e[cluster.member, cv[[cv.itr]]$train ])/length( cluster.member )
		e.neg.signature.names[k] = cluster.member[1] 
		e.signature.member[[ cluster.member[1] ]]= cluster.member
	}
		rownames(e.neg.signature)=e.neg.signature.names

	e.signature = data.frame( t(rbind(e.pos.signature ,  e.neg.signature) ) )

	#e.signature=data.frame( t(e.pos.signature) )

	upper = terms(survData[ cv[[cv.itr]]$train ]~(.), data = e.signature)
	model.aic = step(coxph(survData[ cv[[cv.itr]]$train ] ~(.), data=e.signature), scope=upper, direction="both", k=2) #, trace=FALSE
	

	cv.signature.member[[cv.itr]]=e.signature.member[attr(model.aic$terms, "term.labels")]
	
	}#end of crossvalidation runs

	print(sort(table(unlist(cv.signature.member)),decreasing=TRUE)   )
	save(cv.signature.member,file="KIRC.han.cv.signature.rda")

#}

#survCluster(e.exp, e.surv)



