SubsetterCombinedPIC <-
function(genos,save=NULL,size=100,permutations=100){
  Opt<-function(genos,random.lists,save){
    better.randoms<-matrix(nrow=nrow(random.lists),ncol=ncol(random.lists))
    for(i in 1:(ncol(random.lists))){
      better.randoms[,i]<-SubsetOptimizerPIC(genos,random.lists[,i],save=save,print=FALSE)
    }
    return(better.randoms)
  }
  genos<-as.matrix(genos)
  genos[is.na(genos)]<-0
  n.genos<-ncol(genos)
  result.list<-rep(0,permutations)
  random.lists<-matrix(0,nrow=size,ncol=permutations)
  if(is.null(save)){
  for(i in 1:permutations){
  random.lists[,i]<-sample(colnames(genos),size)
    }
  }else{
    random.lists[1:length(save),]<-save
    for(i in 1:permutations){
  random.lists[(length(save)+1):size,i]<-sample(colnames(genos),size-length(save))
    }
  }
  better.randoms<-matrix(nrow=nrow(random.lists),ncol=ncol(random.lists))
  for(i in 1:permutations){
    better.randoms[,i]<-SubsetOptimizerPIC(genos,random.lists[,i],save=save,print=FALSE)
  }
  for(i in 1:permutations){
    result.list[i]<-PicCalc(genos[,better.randoms[,i]])
  }
  print(c("Best PIC:",max(result.list)))
  best<-which.max(result.list)
  return(better.randoms[,best])
}
