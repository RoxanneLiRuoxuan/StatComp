## -----------------------------------------------------------------------------
Var=c(1,9,49,100)
for(i in Var){
  n<-1000
  u<-runif(n)
  x1<-sqrt(2*i*log(1/(1-u)))
  hist(x1,prob=TRUE,main=paste("var",i))
  y<-seq(0,40,0.1)
  lines(y,y*exp(-y^2/(2*i))/i)
}


## -----------------------------------------------------------------------------
p1=c(0.25,0.4,0.5,0.75)
for(j in p1){
m<-1000
r1<-rnorm(m)
r2<-rnorm(m,3,1)
o<-sample(c(0,1),m,replace=TRUE,prob=c(1-j,j))
z<-o*r1+(1-o)*r2
hist(z,prob=TRUE,ylim=c(0,1),main=paste("p1",j))
y0<-seq(-3,8,0.1)
lines(y0,(j*exp(-y0^2/2)+(1-j)*exp(-(y0-3)^2/2))/sqrt(2*j))
}

## -----------------------------------------------------------------------------
lamda<-c(1,3,5)
for(l in lamda){
q<-1000
t<-10
r<-2
beta<-5
w1<-rpois(q,l*t)
w2<-rgamma(q,r*w1,beta)
hist(w2,main=paste("lamda",l,"compound possion distribution"))
wmean<-mean(w2)
tmean<-l*t*r/beta
wvar<-var(w2)
tvar<-l*t*r*(r+1)/beta^2
cat("lamda",l,"Distributed sample mean:",wmean,"\n")
cat("lamda",l,"True sample mean:",tmean,"\n")
cat("lamda",l,"Distributed sample mean:",wvar,"\n")
cat("lamda",l,"True sample mean:",tvar,"\n")
}

## -----------------------------------------------------------------------------
set.seed(1234)

beta33<-function(x,alpha=3,beta=3){
  n<-10000
  t<-runif(n,min=0,max=x)
  x0<-30*x*(t^2-2*t^3+t^4)
  return(mean(x0))
}
x1<-seq(0.1,0.9,0.1)
y<-numeric(9)
j<-1
for(i in x1){
  y[j]=beta33(i)
  j<-j+1
}

y0<-numeric(9)
k<-1
for(i in x1){
  y0[k]<-pbeta(i,3,3)
  k<-k+1
}


table1=cbind(x1,y,y0)
knitr::kable(table1,align='c',col.names =c('x','simulation','true'))

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(1234)
#  Rayleign<-function(x,sigma,N=10000,anti=TRUE){
#    u<-runif(N/2,min=0,max=x)
#    if (!anti) v <- runif(N/2,min=0,max=x) else v <- x - u
#    u<-c(u,v)
#    y<-x*u*exp(-u^2/2/(sigma^2))/(sigma^2)
#    return(mean(y))
#  }
#  
#  x<-c(1.5,3,10,20)
#  sigma<-2
#  y2<-numeric(length(x))
#  y3<-numeric(length(x))
#  sd1<-numeric(length(x))
#  sd2<-numeric(length(x))
#  mean1<-numeric(length(x))
#  mean2<-numeric(length(x))
#  for(k in 1:length(x)){
#    for(i in 1:1000){
#      y2[i]<-Rayleign(x[k],sigma,anti=TRUE)
#      y3[i]<-Rayleign(x[k],sigma,anti=FALSE)
#    }
#   mean1[k]<-mean(y2)
#   mean2[k]<-mean(y3)
#   sd1[k]<-sd(y2)
#   sd2[k]<-sd(y3)
#  }
#  
#  percent=sd2/sd1
#  
#  table2<-cbind(x,sigma,mean1,sd1,mean2,sd2,percent)
#  knitr::kable(table2,align='c',col.names =c('x','sigma',"Mean value of simulated value of dual variable","Standard deviation of simulated value of dual variable","Mean of simulated values of independent variables","Standard deviation of simulated value of independent variable","Percentage reduction in standard deviation"))
#  

## -----------------------------------------------------------------------------
set.seed(123)
m<-1000
n<-20
u<-sigma<-y<-numeric(m)
for(i in 1:m){
   x<-rchisq(n,df=2)
   u[i]<-mean(x)
   sigma[i]<-sd(x)
   y[i]<-(2<=u[i]-qt(0.025,n-1)*sigma[i]/sqrt(n)&2>=u[i]-qt(0.975,n-1)*sigma[i]/sqrt(n))
}
cp1<-round(mean(y),3)

## -----------------------------------------------------------------------------
set.seed(123)
m<-1000
n<-20
UCL<-y<-numeric(m)
for(i in 1:m){
  x<-rchisq(n,df=2)
  UCL[i]<-(n-1)*var(x)/qchisq(0.05,df=n-1)
  y[i]<-(4<=UCL[i])
}
cp2<-round(mean(y),3)

## -----------------------------------------------------------------------------
set.seed(123)
lambda<-1
m<-1000
n<-c(10,20,50,100,1000)
p.val1<-p.val2<-p.val3<-y<-t<-numeric(m)
tie1<-tie2<-tie3<-numeric(length(n))
for(k in 1:length(n)){
 for(i in 1:m){
   x<-rchisq(n[k],df=lambda)
   t[i]<-sqrt(n[k])*(mean(x)-1)/sd(x)
   p.val1[i]<-2*(1-pt(abs(t[i]),n[k]-1))
  }
  tie1[k]<-mean(p.val1<0.05)
 for(i in 1:m){
   x<-runif(n[k],0,2)
   t[i]<-sqrt(n[k])*(mean(x)-1)/sd(x)
   p.val2[i]<-2*(1-pt(abs(t[i]),n[k]-1))
 }
  tie2[k]<-mean(p.val2<0.05)
 for(i in 1:m){
   x<-rexp(n[k],lambda)
   t[i]<-sqrt(n[k])*(mean(x)-1)/sd(x)
   p.val3[i]<-2*(1-pt(abs(t[i]),n[k]-1))
  }
  tie3[k]<-mean(p.val3<0.05)
}
table<-cbind(n,tie1,tie2,tie3)
knitr::kable(table,col.names=c('sample size','chi_q','uniform','exp'),align='c')

## ----eval=FALSE---------------------------------------------------------------
#  library(MASS)
#  set.seed(1)
#  mu<-c(0,0,0,0)
#  sigma<-matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),nrow=4,ncol=4)
#  m<-1000
#  n<-c(10,20,30,50,100,500)
#  
#  
#  T1ER<-function(data,alpha){
#    nr<-nrow(data)
#    nc<-ncol(data)
#    cdata<-matrix(c(1:nr*nc),nrow=nr,ncol=nc)
#    for(i in 1:nc){
#      cdata[,i]<-data[,i]-rowMeans(data)
#    }
#    sigma0<-cdata%*%t(cdata)/nc
#    a<-t(cdata)%*%solve(sigma0)%*%cdata
#    bstacs<-sum(rowSums(a^{3}))/nc/nc
#    y<-as.integer((nc*bstacs/6)>(qchisq(1-alpha,nr*(nr+1)*(nr+2)/6)))
#    return(y)
#  }
#  
#  tier1<-numeric(6)
#  
#  for(i in 1:length(n)){
#    yy<-numeric(m)
#    for(j in 1:m){
#      x<-mvrnorm(n[i],mu,sigma)
#      yy[j]<-T1ER(data=t(x),alpha=0.05)
#    }
#    tier1[i]<-round(mean(yy),3)
#  }
#  TypeI<-as.vector(tier1)
#  n<-c('10','20','30','50','100','500')
#  table<-cbind(n,TypeI)
#  knitr::kable(t(table),align='c')

## ----eval=FALSE---------------------------------------------------------------
#  library(MASS)
#  set.seed(1)
#  
#  n <- 30
#  m <- 2000
#  
#  epsi <- c(seq(0,0.15,.01),seq(0.15,1,.05))
#  mu<-c(0,0,0,0)
#  sigma1<-matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),nrow=4,ncol=4)
#  sigma2<-matrix(c(100,0,0,0,0,100,0,0,0,0,100,0,0,0,0,100),nrow=4,ncol=4)
#  ou<-numeric(length(epsi))
#  
#  xm<-matrix(c(1:4*n),nrow=4,ncol=n)
#  power<-numeric(m)
#  power_ep<-numeric(length(epsi))
#  
#  for(i in 1:length(epsi)){
#    w<-epsi[i]
#    for(j in 1:m){
#  
#      sigma_d<-sample(c(1,10),size=n,replace=TRUE,prob=c(1-w,w))
#      for(k in 1:n){
#        sigma_dd<-sigma_d[k]
#        if(sigma_dd==1){sigma<-sigma1}
#        else{sigma<-sigma2}
#        xm[,k]<-mvrnorm(1,mu,sigma)
#      }
#  
#      power[j]<-T1ER(data=xm,alpha <- 0.1)
#     }
#    power_ep[i]<-mean(power)
#  }
#  
#  plot(epsi, power_ep, type = "b",
#  xlab = bquote(epsi), ylim = c(0,1),col
#  ='blue',main='epsilon-power')
#  abline(h = .1, lty = 3)
#  se <- sqrt(power_ep * (1-power_ep) / m)
#  lines(epsi, power_ep+se, lty = 4,col='red')
#  lines(epsi, power_ep-se, lty = 4,col='red')

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(1)
#  data(scor, package = "bootstrap")
#  
#  data.cov<-function(data){
#    Nrow<-nrow(data)
#    Ncol<-ncol(data)
#  
#    meansdata<-matrix(NA,nrow=Nrow,ncol=Ncol)
#  
#    for(n in 1:Ncol){
#      meansdata[,n]<-data[,n]-rowMeans(data)
#    }
#    datacov<-meansdata%*%t(meansdata)/Ncol
#    return(datacov)
#  }
#  
#  B<-2000
#  
#  theta_star<-numeric(B)
#  
#  data.Sigma0<-cov(scor)
#  scor.values0<-eigen(data.Sigma0)
#  theta_hat<-scor.values0$values[1]/sum(scor.values0$values)
#  
#  
#  for(i in 1:B){
#    data.ind<-sample(1:88,size=88,replace=TRUE)
#    data.boot<-scor[data.ind,]
#    data.Sigma<-data.cov(t(data.boot))
#    scor.values<-eigen(data.Sigma)
#    theta_star[i]<-scor.values$values[1]/sum(scor.values$values)
#  }
#  theta.bias<-mean(theta_star)-theta_hat
#  theta.sd<-sd(theta_star)
#  table<-cbind(theta.bias,theta.sd)
#  knitr::kable(table,align="c")

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(1)
#  n<-88
#  data.jack<-matrix(NA,ncol=5,nrow=n-1)
#  theta_jack<-numeric(88)
#  
#  scor.Sigma_j<-matrix(NA,ncol=5,nrow=5)
#  
#  for(j in 1:n){
#    scor.jack<-scor[-j,]
#    scor.Sigma_j<-data.cov(t(scor.jack))
#    scor.values_j<-eigen(scor.Sigma_j)
#    theta_jack[j]<-scor.values_j$values[1]/sum(scor.values_j$values)
#  }
#  
#  bias_jack<-(n-1)*(mean(theta_jack)-theta_hat)
#  sd_jack<-sqrt((n-1)*mean((theta_jack-theta_hat)^2))
#  table<-cbind(bias_jack,sd_jack)
#  knitr::kable(table,align="c")

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(1)
#  library(boot)
#  
#  boot.theta<-function(x,ind){
#    data.Sigma0<-cov(x[ind,])
#    patch.values0<-eigen(data.Sigma0)
#    theta_hat<-patch.values0$values[1]/sum(patch.values0$values)
#  }
#  CI_normal<-CI_basic<-CI_percent<-CI_Bca<-matrix(NA,1,2)
#  
#  。
#  scor_theta<-boot(data=scor,statistic=boot.theta, R = 2000)
#  CI<-boot.ci(scor_theta,type=c("norm","basic","perc","bca"))
#  print(CI)

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(1)
#  m<-1000
#  
#  skewfunc<-function(x,ind){
#    skew_x<-mean((x[ind]-mean(x[ind]))^3)
#    skew_y<-mean((x[ind]-mean(x[ind]))^2)^(3/2)
#    skew<-skew_x/skew_y
#    return(skew)
#  }
#  
#  n<-100
#  
#  xnorm.ci.norm<-xnorm.ci.basic<-xnorm.ci.perc<-xchi.ci.norm<-xchi.ci.basic<-xchi.ci.perc<-matrix(NA,m,2)
#  for(i in 1:m){
#    xnorm<-rnorm(n)
#    xchi<-rchisq(n,df=5)
#    xnorm.boot<-boot(data=xnorm,statistic=skewfunc,R=2000)
#    xchi.boot<-boot(data=xchi,statistic=skewfunc,R=2000)
#    norm_CI<-boot.ci(xnorm.boot,type=c("norm","basic","perc"))
#    chi_CI<-boot.ci(xchi.boot,type=c("norm","basic","perc"))
#    xnorm.ci.norm[i,]<-norm_CI$norm[2:3]
#    xnorm.ci.basic[i,]<-norm_CI$basic[4:5]
#    xnorm.ci.perc[i,]<-norm_CI$percent[4:5]
#    xchi.ci.norm[i,]<-chi_CI$norm[2:3]
#    xchi.ci.basic[i,]<-chi_CI$basic[4:5]
#    xchi.ci.perc[i,]<-chi_CI$percent[4:5]
#  }
#  
#  norm.sk<-0
#  chi.sk<-1.265
#  norm.left<-c(mean(xnorm.ci.norm[,1]>=norm.sk),mean(xnorm.ci.basic[,1]>=norm.sk),mean(xnorm.ci.perc[,1]>=norm.sk))
#  norm.right<-c(mean(xnorm.ci.norm[,2]<=norm.sk),mean(xnorm.ci.basic[,2]<=norm.sk),mean(xnorm.ci.perc[,2]<=norm.sk))
#  chi.left<-c(mean(xchi.ci.norm[,1]>=chi.sk),mean(xchi.ci.basic[,1]>=chi.sk),mean(xchi.ci.perc[,1]>=chi.sk))
#  chi.right<-c(mean(xchi.ci.norm[,2]<=chi.sk),mean(xchi.ci.basic[,2]<=chi.sk),mean(xchi.ci.perc[,2]<=chi.sk))
#  norm.cov<-1-(norm.left+norm.right)
#  chi.cov<-1-(chi.left+chi.right)
#  
#  ci.method<-c('ci.norm','ci.basic','ci.percentile')
#  table<-cbind(ci.method,norm.left,norm.right,chi.left,chi.right,norm.cov,chi.cov)
#  knitr::kable(table,align='c')

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(1)
#  s1<-runif(20,min=0,max=10)
#  s2<-rnorm(20)
#  s3<-rt(20,df=1)
#  spearman01<-cor(s1,s2,method='spearman')
#  spearman02<-cor(s3,s2,method='spearman')
#  s<-c(s1,s2)
#  ss<-c(s2,s3)
#  B<-1000
#  spearman1<-numeric(B)
#  spearman2<-numeric(B)
#  
#  for(i in 1:B){
#    u<-sample(1:40,size=20,replace=FALSE)
#    z1<-s[u]
#    z2<-s[-u]
#    spearman1[i]<-cor(z1,z2,method='spearman')
#  }
#  p_value_spearman1<-mean(abs(spearman1)>abs(spearman01))
#  
#  p_value_cor.test.spearman<-cor.test(s1,s2,alternative = c("two.sided"),method="spearman")$p.value
#  p_value_cor.test.pearson<-cor.test(s1,s2,alternative = c("two.sided"),method="pearson")$p.value
#  p_value_cor.test.kendall<-cor.test(s1,s2,alternative = c("two.sided"),method="kendall")$p.value
#  p_value<-cbind(p_value_spearman1,p_value_cor.test.spearman,p_value_cor.test.pearson,p_value_cor.test.kendall)
#  knitr::kable(p_value,align='c',caption='N(0,1) VS U(0,10)')
#  
#  
#  for(i in 1:B){
#    u<-sample(1:40,size=20,replace=FALSE)
#    z1<-ss[u]
#    z2<-ss[-u]
#    spearman2[i]<-cor(z1,z2,method='spearman')
#  }
#  p_value_spearman2<-mean(abs(spearman2)>abs(spearman02))
#  p_value_cor.test.spearman<-cor.test(s3,s2,alternative = c("two.sided"),method="spearman")$p.value
#  p_value_cor.test.pearson<-cor.test(s3,s2,alternative = c("two.sided"),method="pearson")$p.value
#  p_value_cor.test.kendall<-cor.test(s3,s2,alternative = c("two.sided"),method="kendall")$p.value
#  p_value<-cbind(p_value_spearman2,p_value_cor.test.spearman,p_value_cor.test.pearson,p_value_cor.test.kendall)
#  knitr::kable(p_value,align='c',caption='N(0,1) VS t(1)')
#  
#  

## -----------------------------------------------------------------------------
set.seed(1)
library(RANN)
library(boot)


NN3<-function(sample,ind,index1,index2,k0){
  N<-index1+index2
  sample<-data.frame(sample,1)
  s<-sample[ind,]

  kNN<-nn2(data=s,k=k0+1)
  x1<-kNN$nn.idx[1:index1,-1]
  y1<-kNN$nn.idx[(index1+1):N,-1]
  r1<-sum(x1<=index1)
  r2<-sum(y1>index1)
  NNres<-(r1+r2)/(k0*N)
  return(NNres)
}


threetests<-function(oridata,xdata,ydata){
  ind1<-length(xdata)
  ind2<-length(ydata)
 
  t0<-NN3(sample=oridata,index1=ind1,index2=ind2,k0=3)
  NN.boot<- boot(data = oridata, statistic = NN3, R = 10000,sim = "permutation", index1=ind1,index2=ind2,k0=3)
  NN_pvalue<-round(mean(NN.boot$t>=NN.boot$t0),5)

  library(energy)
  boot.obs <- eqdist.etest(x=oridata, sizes=c(ind1,ind2), R=10000)
  Energy_pvalue<- round(boot.obs$p.value,5)
  #Ball检验
  library(Ball)
  Ball_pvalue= round(bd.test(x = xdata, y = ydata, num.permutations=10000)$p.value,5)
  return(c(NN_pvalue,Energy_pvalue,Ball_pvalue))
}

sd<-c(2,10,100)
p_value<-matrix(1,ncol=3,nrow=3)

for(i in 1:length(sd)){
  x<-as.vector(rnorm(20))
  y<-as.vector(rnorm(20,mean=0,sd=sd[i]))
  z<-c(x,y)
  p_value[i,]<-threetests(oridata=z,xdata=x,ydata=y)
}
print(p_value)

sd_value<-c('sd=2','sd=10','sd=100')
table<-cbind(sd_value,p_value)
knitr::kable(table,align='c',col.names=c('sd_value','NN test','Energy test','Ball test'))

## ----eval=FALSE---------------------------------------------------------------
#  
#  set.seed(1)
#  mu<-c(1,2,10)
#  p_value2<-matrix(1,ncol=3,nrow=3)
#  
#  for(i in 1:length(sd)){
#    x<-as.vector(rnorm(20))
#    y<-as.vector(rnorm(20,mean=mu[i],sd=sd[i]))
#    z<-c(x,y)
#    p_value2[i,]<-threetests(oridata=z,xdata=x,ydata=y)
#  }
#  
#  
#  sd_value<-c('mu=1,sd=2','mu=2,sd=5','mu=10,sd=100')
#  table<-cbind(sd_value,p_value2)
#  knitr::kable(table,align='c',col.names=c('sd_value','NN test','Energy test','Ball test'))

## ----eval=FALSE---------------------------------------------------------------
#  
#  set.seed(1)
#  x1<-as.vector(rt(20,df=1))
#  y1<-numeric(20)
#  rind<-sample(1:2,size=20,replace=TRUE,prob=c(0.3,0.7))
#  
#  for(i in 1:20){
#    if(rind[i]==1) y1[i]<-rnorm(1)
#    else y1[i]<-rnorm(1,mean=1,sd=10)
#  }
#  z1<-c(x1,y1)
#  p_value3<-numeric(3)
#  p_value3<-as.vector(threetests(oridata=z1,xdata=x1,ydata=y1))
#  methods<-c('NN test','Energy test','Ball test')
#  
#  
#  table3<-rbind(methods,p_value3)
#  knitr::kable(table3,align='c')

## ----eval=FALSE---------------------------------------------------------------
#  
#  set.seed(1)
#  x2<-as.vector(rt(10,df=1))
#  y2<-as.vector(rt(1000,df=1))
#  z2<-c(x2,y2)
#  p_value4<-numeric(3)
#  p_value4<-as.vector(threetests(oridata=z2,xdata=x2,ydata=y2))
#  
#  
#  methods<-c('NN test','Energy test','Ball test')
#  table4<-rbind(methods,p_value4)
#  knitr::kable(table4,align='c')

## -----------------------------------------------------------------------------
set.seed(1)

target_func<-function(x,theta,eta){
  tar1<-((x-eta)/theta)^2
  tar2<-theta*pi*(1+tar1)
  tar3<-1/tar2
  return(tar3)
}


experiments<-10000
x<-numeric(experiments)

x[1]<-rnorm(1)
reject_times<-0


Cauchyfunc<-function(orix,times){
  x[1]<-orix
  N<-times
  for(i in 2:N){
    u<-runif(1)
 
    y<-rnorm(1,mean=x[i-1],sd=2)
    MHSamde<-target_func(x=y,theta=1,eta=0)
    MHSamda<-target_func(x=x[i-1],theta=1,eta=0)
    if(u<=MHSamde/MHSamda){x[i]<-y}
    else{x[i]<-x[i-1];reject_times<-reject_times+1}
  }
  print(paste('Rejection rate:',reject_times/times))
  return(x)
}
x<-Cauchyfunc(x[1],experiments)

strunc<-1000
x_strunc<-x[strunc+1:experiments]
perc<-seq(0.1,0.9,0.1)
mcmc_perc<-quantile(x_strunc,probs=perc,na.rm=TRUE)
true_perc<-qcauchy(perc,scale=1)
compare_perc<-cbind(mcmc_perc,true_perc)
print(compare_perc)
qqplot(true_perc,mcmc_perc, main="deciles",xlab="Cauchy Deciles", ylab="Sample Deciles")

## -----------------------------------------------------------------------------
set.seed(12345)


GXori<-c(1,0.1)

a<-2
b<-2
n<-10

GXfunc<-function(ori,N){
   GbN<-N
   GX<-matrix(0,ncol=GbN,nrow=2)
   GX[,1]<-ori
for(i in 2:GbN){
  x<-GX[,i-1][1]
  GX[,i][2]<-rbeta(1,shape1=x+a,shape2=n-x+b)
  y<-GX[,i][2]
  GX[,i][1]<-rbinom(1,size=n,prob=y)
}
  return(GX)
}
GX<-GXfunc(GXori,5000)

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(1543231221)
#  
#  GRstat<-function(Cauchydata){
#    n<-ncol(Cauchydata)#the length of every chain is n
#    k<-nrow(Cauchydata)#k chains
#    chainaver<-rowMeans(Cauchydata)
#    allaver<-mean(chainaver)
#    betweenvar<-n*sum((chainaver-allaver)^2)/(k-1)
#    s<-numeric(n)
#    ss<-numeric(k)
#    for(i in 1:k){
#      for(j in 1:n){
#        s[j]<-(Cauchydata[i,j]-chainaver[i])^2
#      }
#      ss[i]<-mean(s)
#    }
#    withinvar<-mean(ss)
#    VAR<-((n-1)*withinvar+betweenvar)/n
#    R_hat<-sqrt(VAR/withinvar)
#    return(R_hat)
#  }
#  
#  x<-numeric(4)
#  x[1]<-rnorm(1,mean=1,sd=1)
#  x[2]<-rnorm(1,mean=2,sd=1)
#  x[3]<-rnorm(1,mean=3,sd=1)
#  x[4]<-rnorm(1,mean=4,sd=1)
#  experi<-10000
#  burn_in<-1000
#  Cdata<-matrix(0,ncol=experi,nrow=length(x))
#  for(i in 1:4){
#    Cdata[i,]<-t(Cauchyfunc(x[i],experi))
#  }
#  Cdata1<-t(apply(Cdata, 1, cumsum))
#     for (i in 1:nrow(Cdata))
#         Cdata1[i,] <- Cdata[i,] / (1:ncol(Cdata))
#  Cdata<-Cdata1[,(burn_in+1):experi]
#  Cauchy_R_hat<-GRstat(Cdata)
#  print(Cauchy_R_hat)
#  r_hat <- rep(0, experi)
#  for (j in (burn_in+1):experi)
#    r_hat[j] <- GRstat(Cdata1[,1:j])
#  plot(r_hat[(burn_in+1):experi], type="l", xlab="", ylab="R")
#  

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(123456)
#  
#  gx<-matrix(0,ncol=4,nrow=2)
#  gx[,1]<-c(1,0.3)
#  gx[,2]<-c(2,0.2)
#  gx[,3]<-c(2,0.1)
#  gx[,4]<-c(1,0.1)
#  
#  GbN<-10000
#  burn_in<-2000
#  
#  GXdata1<-matrix(0,ncol=GbN,nrow=4)
#  GXdata2<-matrix(0,ncol=GbN,nrow=4)
#  for(i in 1:4){
#    GXdata1[i,]<-GXfunc(gx[,i],GbN)[1,]
#    GXdata2[i,]<-GXfunc(gx[,i],GbN)[2,]
#  }
#  
#  GXdata1<-t(apply(GXdata1, 1, cumsum))
#     for (i in 1:nrow(GXdata1))
#         GXdata1[i,] <- GXdata1[i,] / (1:ncol(GXdata1))
#  GXdata2<-t(apply(GXdata2, 1, cumsum))
#     for (i in 1:nrow(GXdata2))
#         GXdata2[i,] <- GXdata2[i,] / (1:ncol(GXdata2))
#  
#  GXdata11<-GXdata1[,(burn_in+1):GbN]
#  GXdata22<-GXdata2[,(burn_in+1):GbN]
#  GXVAR_x<-GRstat(GXdata11)
#  GXVAR_y<-GRstat(GXdata22)
#  print(paste('x:',GXVAR_x))
#  print(paste('y:',GXVAR_y))
#  
#  GX_r_hat <- rep(0, GbN)
#  for (j in (burn_in+1):GbN)
#    GX_r_hat[j] <- GRstat(GXdata1[,1:j])
#  plot(GX_r_hat[(burn_in+1):GbN], type="l", xlab="", ylab="R")
#  abline(h=1.2, lty=2)
#  for (j in (burn_in+1):GbN)
#    GX_r_hat[j] <- GRstat(GXdata2[,1:j])
#  plot(GX_r_hat[(burn_in+1):GbN], type="l", xlab="", ylab="R")
#  abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
set.seed(1)
####(a)
JUAN<-function(k,a,d){
  
}

juan_k<-function(k,a,d){
  if((k<0)|(!k%%1==0)|(!d%%1==0)|(d<1)|(length(a)!=d)|(!is.vector(a)))
    {print("print invalid value")}
  else{
      a_length<-sqrt(a%*%a)
      a1<-(2*k+2)*log(a_length)
      a2<-lgamma((d+1)/2)
      a3<-lgamma(k+3/2)
      b1<-log(factorial(k))
      b2<-k*log(2)
      b3<-log(2*k+1)
      b4<-log(2*k+2)
      b5<-lgamma(k+d/2+1)
      res<-a1+a2+a3-b1-b2-b3-b4-b5
      if(k%%2==0){return(exp(res))}
      else{return(-exp(res))}
 }
}
####(b)
sum_juan_k<-function(k=0,a,d){
  if((k<0)|(!k%%1==0)|(!d%%1==0)|(d<1)|(length(a)!=d)|(!is.vector(a)))
  {print("print invalid value")}
  else{
    x_sum<-0 
    resdu<-100

    while(abs(resdu)>1e-200){
      x<-juan_k(k=k,a=a,d=d)
      x_sum<-x_sum+x
      k<-k+1
      resdu<-x
    }
    return(x_sum)
    }
}
###(c)
a<-as.vector(c(1,2))
d<-length(a)
a_sum<-sum_juan_k(k=0,a,d=2)
print(a_sum)

## -----------------------------------------------------------------------------
set.seed(1)

in_func1<-function(u,m){
  (1+u^2/(m-1))^(-m/2)
}
happy_func<-function(k,a){
 
  l1<-log(2)+lgamma(k/2)
  l2<-0.5*log((pi*(k-1)))+lgamma((k-1)/2)
  c_k_1<-sqrt((a^2*(k-1))/(k-a^2))
  l3<-integrate(in_func1,rel.tol=.Machine$double.eps^0.25
,lower=0,upper=c_k_1,m=k)
  l3<-log(l3$value)
  y1<-exp(l1+l3-l2)
  return(y1)
}

solu_func<-function(k){
  output<-uniroot(function(a){happy_func(k,a)-happy_func(k+1,a)},lower=1,upper=2)
}

k_value<-c(4:25,100,500,1000)
roots<-numeric(length(k_value))
for(j in 1:length(k_value)){
  roots[j]<-solu_func(k_value[j])$root
}
table<-cbind(k_value,roots)
knitr::kable(table,align='c')

## -----------------------------------------------------------------------------
set.seed(1)
lambda2<-0.8
res<-1000
while(res>.Machine$double.eps^0.5){
  lambda1<-lambda2
  lambda2<-(3*(1+lambda1)+3.75)/10
  res<-lambda2-lambda1
}
print(lambda2)

## -----------------------------------------------------------------------------
###EXAMPLE 3
attach(mtcars)
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)
lares<-lapply(formulas,function(x) lm(formula=x,data=mtcars))
##use function rsq and lapply to extract R^2
rsq <- function(mod) summary(mod)$r.squared
r_square1<-lapply(lares,rsq)
Formulas<-c('mpg ~ disp','mpg ~ I(1 / disp)','mpg ~ disp + wt','mpg ~ I(1 / disp) + wt')
#present the result in a table
table1<-cbind(Formulas,r_square1)
knitr::kable(table1,align='c')

## -----------------------------------------------------------------------------
###EXAMPLE 4
set.seed(1)
#use bootstrap to generate data
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})
#use lapply to extract R^2
whp<-lapply(bootstraps,function(x) lm(formula=mpg ~ disp,data=x))
r_square2<-lapply(whp,rsq)
#present the result in a table
table2<-cbind(r_square2)
knitr::kable(table2,align='c')

## -----------------------------------------------------------------------------
set.seed(1)
###a)
#generate a numeric data frame
orinumbers<-data.frame(x=runif(6,min=0,max=10),y=rnorm(6),z=rcauchy(6))
vapply(orinumbers,sd,FUN.VALUE=numeric(1))

## -----------------------------------------------------------------------------
###b)
#generate a mixed data frame
oribnum<-data.frame(
  name=c('Anna','Cindy','Mike','Ross','Joey','Melody'),
  sex=c('F','F','M','M','M','F'),
  height=c(165.0,170.5,188.3,190.1,183.0,174.2),
  weight=c(50.5,56.3,68.0,75.1,65.9,55.4)
  )
#use is.numeric() to identify numerical data
resb<-vapply(oribnum,is.numeric,FUN.VALUE=logical(1))
#use vapply to compute the standard deviation of numerical columns
vapply(oribnum[resb],sd,FUN.VALUE=numeric(1))


## ----eval=FALSE---------------------------------------------------------------
#  set.seed(1)
#  library(parallel)
#  #Calculate the number of computer cores
#  cores <- detectCores(logical = FALSE)
#  
#  #use parSapply in parallel to generate mcsapply
#  mcsapply<-function(x,f){
#    coresnum<- makeCluster(3)
#    median_res<-parSapply(coresnum,x,f)
#    stopCluster(coresnum)
#    return(median_res)
#  }
#  
#  ttestdata<-replicate(10000,t.test(rnorm(20,1,11), runif(30)),simplify=FALSE)
#  system.time(mcsapply(ttestdata,function(x) unlist(x)))
#  system.time(sapply(ttestdata,function(x) unlist(x)))
#  

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(1)
#  library(Rcpp)
#  cppFunction(
#  'NumericMatrix Mygibbs(double x,double y,int N) {
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
#  }')
#  GX_cpp<-Mygibbs(1,0.1,1000)
#  
#  GXfunc<-function(ori,N){
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
#  a<-2
#  b<-2
#  n<-10
#  GXori<-c(1,0.1)
#  GX_r<-GXfunc(GXori,1000)
#  
#  qqplot(GX_cpp[1,],GX_r[1,])
#  qqplot(GX_cpp[2,],GX_r[2,])
#  

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(1)
#  library(microbenchmark)
#  ts <- microbenchmark(Gibbs1=Mygibbs(1,0.1,1000),Gibbs2=GXfunc(GXori,1000))
#  summary(ts)[,c(1,3,5,6)]

