M = read.table("data-raw/gene_matrix.txt")

M = M[, 1:100]
cere=rownames(M)[1:54]
cort=rownames(M)[55:98]
n1=length(cere)
n2=length(cort)

st=matrix(0,(dim(M)[2]*(dim(M)[2]-1)/2),8)
colnames(st)=c("gene1","gene2","theta","lrv","lrv1","lrv2","lrm1","lrm2")
i=0
for (j in 2:dim(M)[2]){
  print(j)

  for (k in 1:(j-1)){
    i=i+1

    st[i,1]=as.integer(j)
    st[i,2]=as.integer(k)

    st[i,4]=var(log(M[,j]/M[,k]))
    st[i,5]=var(log(M[cere,j]/M[cere,k]))
    st[i,6]=var(log(M[cort,j]/M[cort,k]))
    st[i,7]=mean(log(M[cere,j]/M[cere,k]))
    st[i,8]=mean(log(M[cort,j]/M[cort,k]))

    st[i,3]=((n1-1)*st[i,5]+(n2-1)*st[i,6])/((n1+n2-1)*st[i,4])
  }

}
save(st,file="stats_genepairs.R")


ptheta=matrix(0,dim(st)[1],20)
for (i in 1:dim(st)[1]){
  print(i)
  for (k in 1:20){
    perm=sample(1:98)
    st5=var(log(M[perm[1:n1],st[i,1]]/M[perm[1:n1],st[i,2]]))
    st6=var(log(M[perm[(n1+1):(n1+n2)],st[i,1]]/M[perm[(n1+1):(n1+n2)],st[i,2]]))
    ptheta[i,k]=((n1-1)*st5+(n2-1)*st6)/((n1+n2-1)*st[i,4])
  }
}
save(ptheta,file="theta_20rand.R")


FDR=matrix(0,4,2)
colnames(FDR)=c("cutoff","FDR")
FDR[,1]=seq(from=0.6,to=0.9,by=0.1)
for(i in 1:dim(FDR)[1]){
  print (i)
  true=length(which(st[,3]<FDR[i,"cutoff"]))
  rand=sum(apply(ptheta,1,function(x){length(which(x<FDR[i,"cutoff"]))}))/20
  FDR[i,2]=rand/true
}
save(FDR,file="FDR_co06-09.R")
