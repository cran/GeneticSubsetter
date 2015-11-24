SubsetOptimizerMTK <-
function(genos,subset,mat=NULL,power=10,save=NULL,print=TRUE){
  genos<-as.matrix(genos)
  mode(genos)<-'numeric'
  genos[genos==0]<-NA
  subset.genos<-length(subset)
  extra<-setdiff(colnames(genos),subset)
  extra.genos<-length(extra)
  n.genos<-ncol(genos)
  i<-0 #number since last change 
  j<-1 #subset column of interest
  if(is.null(mat)){
    mat<-A.mat(t(genos))
  }
  mat<-(2*(mat-min(mat))/(max(mat)-min(mat)))^power
  for(i in 1:n.genos){
    mat[i,i]<-0
  }
  kin<-sum(mat[subset,subset])/(n.genos^2)
  if(print==TRUE) print(c("Starting Value:",kin))
  while(i<(subset.genos)){
    if(!subset[j]%in%save){
      for(k in 1:extra.genos){
        temp.log<-sum(mat[c(subset[-j],extra[k]),c(subset[-j],extra[k])])/(n.genos^2)
        if(temp.log<kin){
          i<-0
          kin<-temp.log
          old<-subset[j]
          subset<-c(subset[-j],extra[k])
          extra[k]<-old
          if(print==TRUE) print(c("Value:",kin))
        }
      }
    }
    if(j==subset.genos){
      j<-0
      if(print==TRUE) print("Cycle Completed, Returning to Beginning of Subset")
    }
    j<-j+1
    i<-i+1
  }
  if(print==TRUE) print("Done")
  return(subset)
}
