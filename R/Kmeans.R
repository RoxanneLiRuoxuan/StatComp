#' @title K-means function
#' @description K-means: a kind of clustering algorithm; unsupervised learning.
#' @param data the data that you want to classify (numeric)
#' @param k The number of categories you want to classify (numeric)
#' @return the value of the (i) item of the Log-likelihood of beta
#' @examples
#' \dontrun{
#'  MyKmeans(swiss,3)
#' }
#' @export
MyKmeans<-function(data,k){
  if(is.na(data) || is.na(k)){
    stop("Parameters are invalid!!")
  }
  #define a function to compute the Euclidean distance
  Eur_dist<-function(x,y){
    dist<-sqrt(sum((x-y)^2))
    return (dist)
  }
  
  Nrow<-nrow(data)
  Ncol<-ncol(data)
  Iter=TRUE
  ori_index<-sample.int(Nrow,size=k)
  ori_point<-data[ori_index,]
  former_point<-ori_point
  iter_point<-matrix(0,nrow=k,ncol=Ncol)
  #dis_matrix is used to record the distance between every point and every center point.
  dis_matrix<-matrix(0,nrow=Nrow,ncol=k)
  while(Iter){
    for(i in 1:Nrow){
      for(j in 1:k){
        dis_matrix[i,j]<-Eur_dist(data[i,],former_point[j,])
      }
    }
    #locate_matrix is used to record the class each point belongs to
    locate_matrix<-matrix(0,nrow=Nrow,ncol=k)
    for(i in 1:Nrow){
      m<-which.min(dis_matrix[i,])
      locate_matrix[i,m]<-1
    }
    #update the location of every center point
    for(i in 1:k){
      iter_point[i,]<-apply(data[which(locate_matrix[,i] == 1),],2,"mean")
    }
    #judge whether the location of center point has changed
    center.true<-c()
    for(i in 1:k){
      if(all(former_point[i,] == iter_point[i,])==T){
        center.true[i]<-TRUE
      }
    }
    former_point = iter_point
    Iter<-ifelse(all(center.true) == T,F,T)
    }
   out=list()
   out[["center_points"]]<-iter_point
   out[["distance"]]<-dis_matrix
   out[["cluster result"]]<-rep(1,Nrow)
   for(i in 1:Nrow){
     out[["cluster"]][i]<-which(locate_matrix[i,] == 1)
   }
   #return the result, including the center points, distances between every point and every center point and the class each point belongs to.
   return(out)
}

