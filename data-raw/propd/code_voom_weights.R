library(limma)
library(SDMTools)

M=read.table("gene_matrix.txt")
M <- M[, 1:100]

cere=rownames(M)[1:54]
cort=rownames(M)[55:98]
n1=length(cere)
n2=length(cort)
n=n1+n2

counts=t(M[c(cere,cort),])
design=matrix(0,dim(M)[1],2)
design[1:length(cere),1]=rep(1,length(cere))
design[(length(cere)+1):dim(M)[1],2]=rep(1,length(cort))
v=voom(counts, design=design, plot=TRUE)
w=v$weights

ce=which(rownames(M)%in%cere)
co=which(rownames(M)%in%cort)

stw=matrix(0,(dim(M)[2]*(dim(M)[2]-1)/2),7)
colnames(stw)=c("thetaw","lrvw","lrvw1","lrvw2","lrmw1","lrmw2","thetaew")
i=0
for (j in 2:dim(M)[2]){
  print(j)

  for (k in 1:(j-1)){
    i=i+1

    W=w[j,]*w[k,]
    n=sum(W)
    s=sum(W^2)
    p=n-s/n

    n1=sum(W[ce])
    s1=sum(W[ce]^2)
    p1=n1-s1/n1

    n2=sum(W[co])
    s2=sum(W[co]^2)
    p2=n2-s2/n2

    stw[i,5]=wt.mean(log(M[ce,j]/M[ce,k]),W[ce])
    stw[i,6]=wt.mean(log(M[co,j]/M[co,k]),W[co])

    stw[i,2]=wt.var(log(M[,j]/M[,k]),W)
    stw[i,3]=wt.var(log(M[ce,j]/M[ce,k]),W[ce])
    stw[i,4]=wt.var(log(M[co,j]/M[co,k]),W[co])

    stw[i,1]=(stw[i,3]*p1+stw[i,4]*p2)/(stw[i,2]*p)

    if (p1*stw[i,3]>p2*stw[i,4]){
      stw[i,7]=1-(p1*stw[i,3])/(p*stw[i,2])
    } else{
      stw[i,7]=1-(p2*stw[i,4])/(p*stw[i,2])
    }
  }
}
