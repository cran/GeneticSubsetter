SubsetOptimizerPIC <-
function(genos,subset,save=NULL,print=TRUE){
  genos<-as.matrix(genos)
  genos[is.na(genos)]<-0
  subset.genos<-length(subset)
  extra<-setdiff(colnames(genos),subset)
  extra.genos<-length(extra)
  tot.genos<-ncol(genos)
  genos[is.na(genos)]<-0
  m<-nrow(genos)
  i<-0 #number since last change 
  j<-1 #subset column of interest
  pic<-PicCalc(genos[,subset])
  if(print==TRUE) print(c("Starting PIC:",pic))
  join<-array(0,c(2,m,ncol(genos)))
  join[2,,]<-genos
  mat.a<-colSums(join==1)
  mat.b<-colSums(join==-1)
  names<-colnames(genos)
  colnames(mat.a)<-names
  colnames(mat.b)<-names
  mat.a.extras<-mat.a[,extra]
  mat.b.extras<-mat.b[,extra]
  mat.a.short<-mat.a[,subset]
  mat.b.short<-mat.b[,subset]
  while(i<(subset.genos)){
    if(!subset[j]%in%save){
      a<-rowSums(mat.a.short)-mat.a.short[,j]
      b<-rowSums(mat.b.short)-mat.b.short[,j]
      temp.a<-a+mat.a.extras
      temp.b<-b+mat.b.extras
      temp.log<-colMeans(1-(temp.a^2+temp.b^2)/(temp.a+temp.b)^2,na.rm=TRUE)
      if(max(temp.log)>pic){
        i<-0
        pic<-max(temp.log)
        top<-which.max(temp.log)
        old<-subset[j]
        subset[j]<-extra[top]
        extra[top]<-old
        top.gen.a<-mat.a.extras[,top]
        top.gen.b<-mat.b.extras[,top]
        mat.a.extras[,top]<-mat.a.short[,j]
        mat.b.extras[,top]<-mat.b.short[,j]
        mat.a.short[,j]<-top.gen.a
        mat.b.short[,j]<-top.gen.b
      }
    }
    if(j==subset.genos){
      j<-0
      if(print==TRUE) print("Cycle Completed, Returning to Beginning of Subset")
      if(print==TRUE) print(c("PIC:",pic))
    }
    j<-j+1
    i<-i+1
  }
  if(print==TRUE) print("Done")
  return(subset)
}
