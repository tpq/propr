#reminder of R objects st, n1,n2 etc.: ##############################################

M=read.table("gene_matrix.txt")

cere=rownames(M)[1:54]
cort=rownames(M)[55:98]
n1=length(cere)
n2=length(cort)
n=n1+n2

st=matrix(0,(dim(M)[2]*(dim(M)[2]-1)/2),8)
colnames(st)=c("gene1","gene2","theta","lrv","lrv1","lrv2","lrm1","lrm2")
i=0
for (j in 2:dim(M)[2]){
    
    
    for (k in 1:(j-1)){
    	i=i+1
	print(i)
	st[i,1]=as.integer(j)
	st[i,2]=as.integer(k)
    	
	st[i,4]=var(log(M[,j]/M[,k]))
	st[i,5]=var(log(M[cere,j]/M[cere,k]))
	st[i,6]=var(log(M[cort,j]/M[cort,k]))
	st[i,7]=mean(log(M[cere,j]/M[cere,k]))
	st[i,8]=mean(log(M[cort,j]/M[cort,k]))

	st[i,3]=((n1-1)*st[i,5]+(n2-1)*st[i,6])/((n-1)*st[i,4])
    }
    
}

####################################################################################

set=sample(1:dim(st)[1],10000)
mat=matrix(0,length(set),3)
for (i in 1:length(set)){
	mat[i,1]=n1*st[set[i],5]/(n*st[set[i],4])
	mat[i,2]=n2*st[set[i],6]/(n*st[set[i],4])
	mat[i,3]=n1*n2*(st[set[i],8]-st[set[i],7])^2/(n^2*st[set[i],4])
}

library(robCompositions)
ternaryDiag(mat,pch=20,col=rgb(0.1,0.1,0.1,0.1),name=c("group1   ","   group2","between-group"),main="LRV proportions")


thetae=c(1:length(set))
for (i in 1:length(set)){
 	    	
    	print(length(set)-i)
	if ((n1-1)*st[set[i],5]>(n2-1)*st[set[i],6]){
		thetae[i]=1-((n1-1)*st[set[i],5])/((n-1)*st[set[i],4])
	} else{
		thetae[i]=1-((n2-1)*st[set[i],6])/((n-1)*st[set[i],4])
	}
}

ternaryDiag(mat[which(thetae>0.7 | thetae<0.3),],pch=20,col=rgb(0.1,0.1,0.1,0.1),name=c("group1   ","   group2","between-group"),main="Weighted LRV proportions\ncutoffs on theta_e from above and below")


