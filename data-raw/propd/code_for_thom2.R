M=read.table("gene_matrix.txt")

cere=rownames(M)[1:54]
cort=rownames(M)[55:98]
n1=length(cere)
n2=length(cort)


altlrv=function(x,y,a){
  N=length(x)
  s=sum((x^a/mean(x^a)-y^a/mean(y^a))^2)/((N-1)*a^2)
  return(s)
}

lrvMa=function(x,y,a){
  N=length(x)
  s=sum((x/mean(x)-y/mean(y))^2)/(a^2*(N-1))
  return(s)
}


a=0.1
#pre-calculate all x^a
Ma=M^a
ast=matrix(0,(dim(M)[2]*(dim(M)[2]-1)/2),8)

colnames(ast)=c("gene1","gene2","atheta","alrv","alrv1","alrv2","alrm1","alrm2")
i=0
for (j in 2:dim(M)[2]){
  print(j)

  for (k in 1:(j-1)){
    i=i+1

    ast[i,1]=as.integer(j)
    ast[i,2]=as.integer(k)

    ast[i,4]=lrvMa(Ma[,j],Ma[,k],a)
    ast[i,5]=lrvMa(Ma[cere,j],Ma[cere,k],a)
    ast[i,6]=lrvMa(Ma[cort,j],Ma[cort,k],a)
    ast[i,7]=mean(Ma[cere,j]-Ma[cere,k])/a
    ast[i,8]=mean(Ma[cort,j]-Ma[cort,k])/a

    ast[i,3]=((n1-1)*ast[i,5]+(n2-1)*ast[i,6])/((n1+n2-1)*ast[i,4])
  }

}
save(ast,file="altstats_genepairs.R")

#do all permutations beforehand:
numperms=100
perms=matrix(0,numperms,dim(M)[1])
for (k in 1:numperms){
  perms[k,]=sample(1:dim(M)[1])
}

#optionally (very slow) we could include (uncorrected) pvalues (but would need more permutations maybe, and do it for a subset only)
#pvals=matrix(0,dim(ast)[1],1)

FDR=matrix(0,10,4)
colnames(FDR)=c("cutoff","randcounts","truecounts","FDR")
FDR[,1]=seq(from=0.55,to=1,by=0.05)

for (i in 1:dim(ast)[1]){

  print(i)
  ptheta=c(1:numperms)
  for (k in 1:100){

    perm=perms[k,]

    ast5=lrvMa(Ma[perm[1:n1],ast[i,1]],Ma[perm[1:n1],ast[i,2]],a)
    ast6=lrvMa(Ma[perm[(n1+1):(n1+n2)],st[i,1]],Ma[perm[(n1+1):(n1+n2)],st[i,2]],a)
    ptheta[k]=((n1-1)*ast5+(n2-1)*ast6)/((n1+n2-1)*ast[i,4])
    for (j in 1:dim(FDR)[1]){
      if (ptheta[k]<FDR[j,"cutoff"]){
        FDR[j,"randcounts"]=FDR[j,"randcounts"]+1
      }
    }

  }
  #pvals[i]=which(sort(c(ptheta,ast[i,3]))==ast[i,3])/(numperms+1)
}
FDR[,"randcounts"]=FDR[,"randcounts"]/numperms

for(i in 1:dim(FDR)[1]){
  print (i)
  FDR[i,"truecounts"]=length(which(ast[,3]<FDR[i,"cutoff"]))
  FDR[i,"FDR"]=FDR[i,"randcounts"]/FDR[i,"truecounts"]
}
save(FDR,"file=FDR_co055-1.R")

#optionally, qvalues could be estimated by
#qvals[i]=min(FDR[which(FDR[,"cutoff"]>ast[i,3]),"FDR"])
