###############################################################################################
#function to calculate sharpness of groups based on method developed by Pillar, 1999 Ecology
###############################################################################################

sharpness<- function(x,k,method= "euclidean",iterations= 999, method.group= "complete"){
  if(!is.matrix(x)){
    stop("x must be a matrix object")
  }
  if(k<=1){
    stop("number of groups must be greater than 1")
  }
  #step 2
  dist_mat<-vegan::vegdist(x,method=method)^2
  #step 3
  ref_partition<-cutree(hclust(d = dist_mat,method = method.group), k = k) #reference partition
  
  prob<- numeric(length=iterations) #object to receive the results of bootstrap procedure
  G_obsAll<- numeric(length=iterations) #object to receive all Gobs in boots procedure
  null_G_All<-numeric(length=iterations) #object to receive all G0
  for(j in 1:iterations){ #inicio da iteração do bootstrap
    #step 4 bootstrap sample
    boot_samp<-x[sample(1:nrow(x),nrow(x), 
                                  replace = T),] #bootstrap sample
    rownames(boot_samp)<-paste("comm",(nrow(x)+1):(nrow(x)+nrow(x)),sep="")
    #step 5
    joint_dist<-as.matrix(vegan::vegdist(rbind(x,boot_samp),method = method)^2)
    #step 6
    boot_partition<-cutree(hclust(d = (vegan::vegdist(boot_samp,method = method)^2),method = method.group), k = k) 
    #i=2
    levels_boot<-(k+1):(k+k)
    for (i in 1:length(levels(as.factor(boot_partition)))){ #generalization to substitute the names of groups in ref by other names in boot partition
      boot_partition<-gsub(pattern = i,replacement = levels_boot[i],x = boot_partition)
    }
    boot_partition<-as.numeric(boot_partition) #partition in bootstrap sample 
    names(boot_partition)<- paste("comm",(nrow(x)+1):(nrow(x)+nrow(x)),sep="")
    int_groups<-interaction(levels(as.factor(ref_partition)), #possible interactions among groups of boots and ref sample
                            levels(as.factor(boot_partition)),sep = ":") 
    matrix_Q<-matrix(NA,nrow = length(levels(as.factor(ref_partition))), ncol = length(levels(as.factor(boot_partition))), 
                     dimnames= list(paste("group",levels(as.factor(ref_partition)),sep=""),paste("group",levels(as.factor(boot_partition)),sep="")),byrow=FALSE) #matrix to receive the values of contrasts 
    for(i in 1:length(levels(int_groups))){ #iterate among the groups in ref and boots sample to calculate matrix of contrasts - matrix Q
      line<- c(names(which(ref_partition==substr(levels(int_groups)[i],start = 1,stop = 1))),names(which(boot_partition==substr(levels(int_groups)[i],start = 3,stop = 3))))
      column<-c(names(which(ref_partition==substr(levels(int_groups)[i],start = 1,stop = 1))),names(which(boot_partition==substr(levels(int_groups)[i],start = 3,stop = 3))))
      matrix_Tj<-as.matrix(joint_dist[line,column]) #matrix of contrasts
      matrix_Tj<-ifelse(lower.tri(matrix_Tj,diag = T)==TRUE,matrix_Tj,0)
      rownames(matrix_Tj)<- line
      colnames(matrix_Tj)<- column
      matrix_Wr<- as.matrix(joint_dist[names(which(ref_partition==substr(levels(int_groups)[i],start = 1,stop = 1))),
                                       names(which(ref_partition==substr(levels(int_groups)[i],start = 1,stop = 1)))])
      matrix_Wr<-ifelse(lower.tri(matrix_Wr,diag = T)==TRUE,matrix_Wr,0) #matrix of within group variation in reeference sample
      rownames(matrix_Wr)<- names(which(ref_partition==substr(levels(int_groups)[i],start = 1,stop = 1)))
      colnames(matrix_Wr)<- names(which(ref_partition==substr(levels(int_groups)[i],start = 1,stop = 1)))
      matrix_Wb<- as.matrix(joint_dist[names(which(boot_partition==substr(levels(int_groups)[i],start = 3,stop = 3))),
                                       names(which(boot_partition==substr(levels(int_groups)[i],start = 3,stop = 3)))])
      matrix_Wb<-ifelse(lower.tri(matrix_Wb,diag = T)==TRUE,matrix_Wb,0) #matrix of within group variation in bootstrap sample
      rownames(matrix_Wb)<- names(which(boot_partition==substr(levels(int_groups)[i],start = 3,stop = 3)))
      colnames(matrix_Wb)<- names(which(boot_partition==substr(levels(int_groups)[i],start = 3,stop = 3)))
      Tj<- sum(matrix_Tj)/nrow(matrix_Tj) #calculation of the sum of boots sample and refs sample- Tj quantity in Pillar
      Wr<- sum(matrix_Wr)/nrow(matrix_Wr)
      Wb<- sum(matrix_Wb)/nrow(matrix_Wb)
      matrix_Q[which(paste("group",as.numeric(substr(levels(int_groups)[i],start = 1,stop = 1)),sep="")==rownames(matrix_Q)),
               which(paste("group",as.numeric(substr(levels(int_groups)[i],start = 3,stop = 3)),sep="")==colnames(matrix_Q))]<- (Tj - (Wb + Wr)) #matrix containing all values of contrasts
    }
    
    #obtaining S and T  total
    matrix_Q.min<- matrix_Q[,clue::solve_LSAP(matrix_Q, maximum = FALSE)]
    S<- sum(diag(matrix_Q.min))
    Ttotal<- sum((vegan::vegdist(rbind(x,boot_samp),method=method)^2))/(2*nrow(x))
    #calculating G* - step 8
    G_obs<- 1-(S/Ttotal)
    
    #calculating G0
    null_int_groups<-numeric(length=nrow(matrix_Q.min))
    for(i in 1:ncol(matrix_Q.min)){
      null_int_groups[i]<- paste(substr(rownames(matrix_Q.min)[i],start = 6,stop = 6),
                                 substr(colnames(matrix_Q.min)[i],start = 6,stop = 6),sep=":")  
    }
    null_int_groups<-as.factor(null_int_groups) #interaction groups for null bootstrap 
    null_Tj_list<-vector(mode = "list",length=length(null_int_groups)) #object to receive null Tj for each interaction
    null_contrast<- numeric(length = length(null_int_groups))
    names(null_contrast)<- null_int_groups
    null_Tj_list<-vector(mode = "list",length=length(null_int_groups))
    
    for(i in 1:length(null_int_groups)){ #iterate among the groups in ref and boots sample to calculate matrix of contrasts - matrix Q
      n_boot<-length(which(boot_partition==substr(null_int_groups[i],start = 3,stop = 3))) #number of elements in boot group
      null_boot<-sample(names(which(ref_partition==substr(levels(null_int_groups)[i],start = 1,stop = 1))),size = n_boot,replace = T) #null bootstrap sample in nearest reference group
      ref<- names(which(ref_partition==substr(levels(null_int_groups)[i],start = 1,stop = 1)))
      null_lineCol<- c(null_boot,ref)
      null_Tj_list[[i]]<- null_boot
      null_matrix_Tj<-as.matrix(joint_dist[null_lineCol,null_lineCol]) #matrix of contrasts between boots and nearest sample in ref 
      null_matrix_Tj<-ifelse(lower.tri(null_matrix_Tj,diag = T)==TRUE,null_matrix_Tj,0)
      rownames(null_matrix_Tj)<- null_lineCol
      colnames(null_matrix_Tj)<- null_lineCol
      null_matrix_Wb<- as.matrix(joint_dist[null_boot,null_boot])
      null_matrix_Wb<-ifelse(lower.tri(null_matrix_Wb,diag = T)==TRUE,null_matrix_Wb,0) #matrix of within group variation in reeference sample
      rownames(null_matrix_Wb)<- null_boot
      colnames(null_matrix_Wb)<- null_boot
      null_matrix_Wr<- as.matrix(joint_dist[ref,ref])
      null_matrix_Wr<-ifelse(lower.tri(null_matrix_Wr,diag = T)==TRUE,null_matrix_Wr,0) #matrix of within group variation in bootstrap sample
      rownames(null_matrix_Wr)<- ref
      colnames(null_matrix_Wr)<- ref
      null_Tj<- sum(null_matrix_Tj)/nrow(null_matrix_Tj) #calculation of the sum of boots sample and refs sample- Tj quantity in Pillar
      null_Wr<- sum(null_matrix_Wr)/nrow(null_matrix_Wr)
      null_Wb<- sum(null_matrix_Wb)/nrow(null_matrix_Wb)
      null_contrast[i]<-(null_Tj - (null_Wb + null_Wr)) #matrix containing all values of contrasts
    }
    null_S<- sum(null_contrast)
    null_boot_samp<-x[do.call("c",null_Tj_list),]
    null_Ttotal_matrix<-vegan::vegdist(rbind(x,null_boot_samp),method=method)^2
    null_Ttotal_matrix<-ifelse(lower.tri(null_Ttotal_matrix)==T,null_Ttotal_matrix,0)
    rownames(null_Ttotal_matrix)<- c(rownames(x),do.call("c",null_Tj_list))
    colnames(null_Ttotal_matrix)<- c(rownames(x),do.call("c",null_Tj_list))
    null_Ttotal<- sum(null_Ttotal_matrix)/(nrow(as.matrix(null_Ttotal_matrix)))
    null_G<-(1-(null_S/null_Ttotal)) #null G -  G0
    null_G_All[j]<- null_G
    G_obsAll[j]<-G_obs
    prob[j]<-ifelse(G_obs>=null_G,1,0)
  }
  p.value<- sum(prob)/iterations
  mean_Gobs<- mean(G_obsAll)
  sd_Gobs<- sd(G_obsAll)
  mean_G0<-mean(null_G_All)
  sd_G0<-sd(null_G_All)
  result_sharpness<-cbind(mean_Gobs,sd_Gobs,mean_G0,sd_G0,p.value)
  return(result_sharpness)
}

#example
run.example= FALSE
if(run.example==TRUE){
  matrix_test<- matrix(c(17,14,27,21,16,5,9,8,5,0,5,8,0,0,10),nrow=5, ncol=3,
                       dimnames= list(paste("comm",1:5, sep=""),
                                      paste("var",1:3, sep="")))
  resul_Ex<-sharpness(x = matrix_test,k = 2, method = "euclidean",iterations = 10000)
}
