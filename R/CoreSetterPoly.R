CoreSetterPoly<-
function(genos,ploidy=2,save=NULL){
  if(length(setdiff(save,rownames(genos)))>0) stop("One or more saved genotypes not present in genos")
  mode(genos)<-"numeric"
  n.genos<-nrow(genos)
  print(paste(n.genos," Genotypes"))
  length.save<-length(save)
  if(length.save>1){
    stop<-length(save)
  }else{
    stop<-2
  }
  ###Subsetting using HET criterion
  result.list<-matrix(0,ncol=3,nrow=n.genos)
  colnames(result.list)<-c("Rank","Individual","HET")
  result.list[,1]<-c(n.genos:1)
  #result.list[1,3]<-HET(genos) #Fix
  m<-ncol(genos)/ploidy
  print(paste(m," Markers"))
  ##Find the highest number of alleles per marker
  alleles.per.marker<-rep(0,m)
  for(i in 1:m){
    alleles.per.marker[i]<-length(unique(na.omit(c(genos[,(i*ploidy-ploidy+1):(i*ploidy)]))))
  }
  max.alleles<-max(alleles.per.marker)
  ##Build an array for efficient processing
  arr<-array(0,dim=c(m,n.genos,max.alleles))
  for(i in 1:m){
    locus.genos<-genos[,(i*ploidy-ploidy+1):(i*ploidy)]
    locus.alleles<-unique(na.omit(c(locus.genos)))
    for(j in 1:n.genos){
      allele.value<-ploidy/length(na.omit(locus.genos[j,]))
      for(k in 1:ploidy){
        arr[i,j,which(locus.alleles == locus.genos[j,k])]<-arr[i,j,which(locus.alleles == locus.genos[j,k])]+allele.value
      }
    }
  }
  names<-rownames(genos)
  #Find value of full set
  temp.mat<-apply(arr,3,rowSums)
  result.list[1,3]<-mean(1-apply(temp.mat^2,1,sum)/apply(temp.mat,1,sum)^2,na.rm=TRUE)


  ##Systematically remove genotypes
  for(i in 1:(n.genos-stop)){
    temp.arr<-array(0,dim=c(m,n.genos-i+1,max.alleles))
    for(j in 1:max.alleles){
      temp.arr[,,j]<-rowSums(arr[,,j])-arr[,,j] #Temp.arr: each cell corresponds to the number of that allele at that locus if that genotype were removed. 
    }
    temp.log<-colMeans(1-apply(temp.arr^2,c(1,2),sum)/apply(temp.arr,c(1,2),sum)^2,na.rm=TRUE) #Finds which genotype, if removed, would result in the lowest HET
    if(length(save)!=0){
      temp.log[which(names %in% save)]<-0
    }
    remove<-which.max(temp.log)
    result.list[i,2]<-names[remove]
    result.list[i+1,3]<-max(temp.log)
    arr<-arr[,-remove,]
    names<-names[-remove]
  }
    ##Clean Up
    if(length.save>1){
      result.list[(n.genos-stop+1):n.genos,1]<-1
      result.list[(n.genos-stop+1):n.genos,2]<-save
    }else{
      result.list[(n.genos-1):n.genos,1]<-1
      result.list[(n.genos-1):n.genos,2]<-c(names)
    }
  return(result.list)
}