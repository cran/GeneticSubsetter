SubsetterCombinedMTK <-
function(genos,save=NULL,size=100,power=10,permutations=100,print=TRUE){
  genos<-as.matrix(genos)
  mode(genos)<-'numeric'
  genos[genos==0]<-NA
  n.genos<-ncol(genos)
  mat<-A.mat(t(genos))
  mat.adj<-(2*(mat-min(mat))/(max(mat)-min(mat)))^power
  for(i in 1:n.genos){
    mat.adj[i,i]<-0
  }
  genos<-as.matrix(genos)
  temp<-SubsetterMTK(genos,mat=mat,power=power,save=save)
  set<-SubsetOptimizerMTK(genos,subset=temp[(nrow(temp)-size+1):nrow(temp),2],mat=mat,power=power,save=save,print=FALSE)
  score<-sum(mat.adj[set,set])/n.genos^2
  if(print==TRUE) print(score)
  if(permutations>1){
    for(i in 1:permutations-1){
      #randomize<-colnames(mat.adj)[sample(1:n.genos,n.genos)]
      if(is.null(save)){
        subset<-sample(colnames(mat.adj),size)
      }else{
        subset<-c(save,sample(setdiff(colnames(mat.adj),save),size-length(save)))
      }
      test<-SubsetOptimizerMTK(genos,subset=subset,mat=mat,power=power,save=save,print=FALSE)
      score.test<-sum(mat.adj[test,test])/n.genos^2
      if(print==TRUE) print(score.test)
      if(score.test<score){
        set<-test
        score<-score.test
      }
    }
  }
  print(c("Best Value:",score))
  return(set)
}
