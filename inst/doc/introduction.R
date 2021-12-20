## ----eval=FALSE---------------------------------------------------------------
#  NumericMatrix Mygibbs(double x,double y,int N) {
#    NumericMatrix GX(2,N);
#    int a=2,b=2,n=10;
#    GX(0,0)=x;
#    GX(1,0)=y;
#    for(int i=1;i<N;i++){
#      double x0=GX(0,i-1);
#      GX(1,i)=rbeta(1,x0+a,n-x0+b)[0];
#      double y0=GX(1,i);
#      GX(0,i)=rbinom(1,n,y0)[0];
#    }
#    return GX;
#  }
#  
#  ##vs r function
#  GXfunc<-function(ori,N){
#     a=2
#     b=2
#     n=10
#     GbN<-N
#     GX<-matrix(0,ncol=GbN,nrow=2)
#     GX[,1]<-ori
#  for(i in 2:GbN){
#    x<-GX[,i-1][1]
#    GX[,i][2]<-rbeta(1,shape1=x+a,shape2=n-x+b)
#    y<-GX[,i][2]
#    GX[,i][1]<-rbinom(1,size=n,prob=y)
#  }
#    return(GX)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  LCV<-function(X,Y,h,n){
#    res0<-matrix(0,ncol=n,nrow=length(h))
#    res<-numeric(length(h))
#    for(l in 1:length(h)){
#     for(i in 1:n){
#      Kernel <- dnorm(x = X[i], mean = X[-i], sd = h[l])
#      beta<-nlm(f = function(beta) {
#        sum(Kernel * (Y[-i] * exp(beta[1] + beta[2] * (X[-i] - X[i]))
#                      - ( beta[1] + beta[2] * (X[-i] - X[i]))))
#      },p = c(0, 0))$estimate[1]
#      res0[l,i]<-(beta - exp(beta) * Y[i])
#      }
#    }
#    res<-rowSums(res0)
#    opt_h<-h[which.max(res)]
#    plot(h, res, type = "o")
#    abline(v = opt_h, col = 2)
#    return(opt_h)
#  }
#  

## ----eval=FALSE---------------------------------------------------------------
#  likelifunc<-function(x,X,Y,h,n){
#    opt_h<-LCV(X,Y,h,n)
#  
#    fitNlm.beta0 <- sapply(x, function(x) {
#      K <- dnorm(x = x, mean = X, sd =opt_h)
#      nlm(f = function(beta) {
#        sum(K * (Y * exp(beta[1] + beta[2] * (X - x))
#                 - (beta[1] + beta[2] * (X - x))))
#      }, p = c(0, 0))$estimate[0]
#    })
#    fitNlm.beta1 <- sapply(x, function(x) {
#      K <- dnorm(x = x, mean = X, sd =opt_h)
#      nlm(f = function(beta) {
#        sum(K * (Y * exp(beta[1] + beta[2] * (X - x))
#                 - (beta[1] + beta[2] * (X - x))))
#      }, p = c(0, 0))$estimate[1]
#    })
#    y<-exp(fitNlm.beta0+fitNlm.beta1*x)
#    plot(x,y, col = 2, lwd = 2, lty = 2)
#    }

## ----eval=FALSE---------------------------------------------------------------
#  n <- 200
#  truebeta <- c(2,2)
#  lambda <- function(x) exp(truebeta[1] + truebeta[2] * x)
#  #generate sample from the exponential generalized linear model
#  X <- runif(n = n, -3, 3)
#  Y <- rexp(n = n, rate = lambda(X))
#  h <- seq(0.1,1, by = .1)
#  LCV(X,Y,h,n)
#  likelifunc(x,X,Y,h,n)

## ----eval=FALSE---------------------------------------------------------------
#  MyKmeans<-function(data,k){
#    if(is.na(data) || is.na(k)){
#      stop("Parameters are invalid!!")
#    }
#    #define a function to compute the Euclidean distance
#    Eur_dist<-function(x,y){
#      dist<-sqrt(sum((x-y)^2))
#      return (dist)
#    }
#  
#    Nrow<-nrow(data)
#    Ncol<-ncol(data)
#    Iter=TRUE
#    ori_index<-sample.int(Nrow,size=k)
#    ori_point<-data[ori_index,]
#    former_point<-ori_point
#    iter_point<-matrix(0,nrow=k,ncol=Ncol)
#    #dis_matrix is used to record the distance between every point and every center point.
#    dis_matrix<-matrix(0,nrow=Nrow,ncol=k)
#    while(Iter){
#      for(i in 1:Nrow){
#        for(j in 1:k){
#          dis_matrix[i,j]<-Eur_dist(data[i,],former_point[j,])
#        }
#      }
#      #locate_matrix is used to record the class each point belongs to
#      locate_matrix<-matrix(0,nrow=Nrow,ncol=k)
#      for(i in 1:Nrow){
#        m<-which.min(dis_matrix[i,])
#        locate_matrix[i,m]<-1
#      }
#      #update the location of every center point
#      for(i in 1:k){
#        iter_point[i,]<-apply(data[which(locate_matrix[,i] == 1),],2,"mean")
#      }
#      #judge whether the location of center point has changed
#      center.true<-c()
#      for(i in 1:k){
#        if(all(former_point[i,] == iter_point[i,])==T){
#          center.true[i]<-TRUE
#        }
#      }
#      former_point = iter_point
#      Iter<-ifelse(all(center.true) == T,F,T)
#      }
#     out=list()
#     out[["center_points"]]<-iter_point
#     out[["distance"]]<-dis_matrix
#     out[["cluster result"]]<-rep(1,Nrow)
#     for(i in 1:Nrow){
#       out[["cluster"]][i]<-which(locate_matrix[i,] == 1)
#     }
#     #return the result, including the center points, distances between every point and every center point and the class each point belongs to.
#     return(out)
#  }

