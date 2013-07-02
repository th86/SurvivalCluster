

load("./c0.6_10_100/BLCA.cv.signature.rda")

load("panCan12v47.probenum.rda")
datasetList=names(N)

target.member=cv.signature.member[[2]]

candidate.member=list()
for( k in 1:length(target.member) ){

target.set = unlist(target.member[k])

	for( dataset.id in 1:length(datasetList) ){		#Scan across all cancer types

							
		load(  paste("./c0.6_10_100/", datasetList[dataset.id], ".cv.signature.rda",sep="")     )
		
		candidate.count=0
		for( cv.id in 1:length(cv.signature.member)    ){  					#Scan across all cross validation runs
			if(length(cv.signature.member)>0){
			for( cluster.id in 1:length( cv.signature.member[[cv.id]]    )  ){			#Scan across all clusters
	  			intersect.set = intersect( unlist(cv.signature.member[[cv.id]][cluster.id]), target.set)
				if(length(intersect.set)>5){					
				candidate.count=candidate.count+1		
				candidate.member[[candidate.count]]=list(  member=intersect.set, cv=cv.id, cluster=cluster.id     )	
				
				}
			}
			}
		
		}
		if(candidate.count>3)
		cat(datasetList[dataset.id], "cluster #",k,"count= ",  candidate.count,"\n")		
		
	}
	


}
