### Set the parameters
desv=sqrt(3)
n=30
b0= 216.694
b1=2.947

M=100
T<-matrix(0,M,18)
 colnames(T)<-c("b0.N","b1.N","Var.N","Li.b0.N","Ls.b0.N","Li.b1.N","Ls.b1.N",
                "Li.Var.N","Ls.Var.N","b0.T","b1.T","Var.T","Li.b0.T","Ls.b0.T",
                "Li.b1.T","Ls.b1.T","Li.Var.T","Ls.Var.T")

#Estudio de simulación
  set.seed(300) 
   x<-rnorm(n,25,4)
   u<-b0+b1*x
   Ys<-matrix(0,n,M)
   library(extraDistr) 
   for (j in 1:M) {
      Ys[,j]<-rslash(n,u,sigma = desv)}
   
par(mfrow=c(3,3))

for (i in 1:9) {
  b<-as.numeric(boxplot(Ys[,i],plot = F)$out)
  a<- which(Ys[,i]%in%b) 
 plot(x,Ys[,i],ylab = paste("y",i,sep = ""))
 points(x[a],Ys[a,i],col="red",coltext = "brown")  
}

   
library(MASS)
   ################################ Bayesian part
##  Normal regression (classic)
L=10000
B0.n<-matrix(0,L,100)
B1.n<-matrix(0,L,100)
var.n<-matrix(0,L,100)
B0.n[1,]<-100
B1.n[1,]<-4
var.n[1,]<-4
X<-as.matrix(cbind.data.frame(1,x))

for (i in 1:M) {
  for(h in 2:L){
    mu=B0.n[h-1,i]+B1.n[h-1,i]*x
    Mv<-solve(t(X)%*%X)*var.n[h-1,i]
    mu.beta=solve(t(X)%*%X)%*%t(X)%*%Ys[,i]
    b<-mvrnorm(1,mu.beta,Sigma = Mv)
    B0.n[h,i]<-b[1]
    B1.n[h,i]<-b[2]
    
    mu=B0.n[h,i]+B1.n[h,i]*x
    a.sigma=0.5*length(x)
    b.sigma=0.5*sum((Ys[,i]-mu)**2)
    var.n[h,i]=1/rgamma(1,a.sigma,b.sigma)
  }}
## Gráficos
par(mfrow=c(3,3))
for (i in 1:9) {
print("Gráfico de cadenas B0")  

ts.plot(B0.n[,i],col="red",ylab="Bo errores normales")  
abline(b0,0)
}

par(mfrow=c(3,3))
print("Gráfico de cadenas B1") 
for (i in 1:9) {
  
  ts.plot(B1.n[,i],col="red",ylab="B1 errores normales")  
  abline(b1,0)
  }
par(mfrow=c(3,3))
print("Gráfico de cadenas var caso normal")  
for (i in 1:9) {

  ts.plot(var.n[,i],col="red",ylab="Var errores normales")  
  abline(3,0)
  }

par(mfrow=c(3,3))
for (i in 1:9) {
  print("Histograma de cadenas B0")  
  
  hist(B0.n[-c(1:1000),i],col="red",xlab="Distr. posteriori B0-caso normal",main = "")  
 abline(v=b0)
  }

par(mfrow=c(3,3))
for (i in 1:9) {
  print("Histograma de cadenas B0")  
  
  hist(B1.n[-c(1:1000),i],col="red",xlab="Distr. posteriori B1-caso normal",main = "")  
  abline(v=b1)
}
par(mfrow=c(3,3))
for (i in 1:9) {
  print("Histograma de cadenas B0")  
  
  hist(var.n[-c(1:1000),i],col="red",xlab="Distr. posteriori Var- caso normal",main = "")  
  abline(v=3)
}

### t-student regression
   L=10000
   B0<-matrix(0,L,100)
   B1<-matrix(0,L,100)
   var<-matrix(0,L,100)
   B0[1,]<-16
   B1[1,]<-2
   var[1,]<-4
   v<-3
   X<-as.matrix(cbind.data.frame(1,x))
  
  for (i in 1:M) {
   for(h in 2:L){
     mu=B0[h-1,i]+B1[h-1,i]*x
     w1<-rgamma(n,(v+1)/2,(v+(Ys[,i]-mu)^2/var[h-1,i])/2)
     w<-diag(w1)
     
     Mv<-solve(t(X)%*%w%*%X)*var[h-1,i]
     mu.beta=solve(t(X)%*%w%*%X)%*%t(X)%*%w%*%Ys[,i]
     b<-mvrnorm(1,mu.beta,Sigma = Mv)
     B0[h,i]<-b[1]
     B1[h,i]<-b[2]
     
     mu=B0[h,i]+B1[h,i]*x
     a.sigma=0.5*length(x)
     b.sigma=0.5*sum((Ys[,i]-mu)**2*w1)
     var[h,i]=1/rgamma(1,a.sigma,b.sigma)
   }}
## graphics   
   par(mfrow=c(3,3))
   print("Gráfico de cadenas B0") 
   for (i in 1:9) {
     ts.plot(B0[-c(1:1000),i],col="red",ylab="Bo errores t-student")  
   abline(b0,0)
     }
   
   par(mfrow=c(3,3))
   print("Gráfico de cadenas B1") 
   for (i in 1:9) {
     
     ts.plot(B1[-c(1:1000),i],col="red",ylab="B1 errores t-student")  
     abline(b1,0)
     }
   par(mfrow=c(3,3))
   print("Gráfico de cadenas var caso normal")  
   for (i in 1:9) {
     
     ts.plot(var[-c(1:1000),i],col="red",ylab="Var errores t-student")  
     abline(3,0)
     }
   #############
   par(mfrow=c(3,3))
   for (i in 1:9) {
     print("Histograma de cadenas B0")  
     
     hist(B0[-c(1:1000),i],col="red",xlab="Distr. posteriori B0-caso t-student",main = "")  
     abline(v=b0)
   }
   
   par(mfrow=c(3,3))
   for (i in 1:9) {
     print("Histograma de cadenas B1")  
     
     hist(B1[-c(1:1000),i],col="red",xlab="Distr. posteriori B1-caso t-student",main = "")  
     abline(v=b1)
   }
   par(mfrow=c(3,3))
   for (i in 1:9) {
     print("Histograma de cadenas B0")  
     
     hist(var[-c(1:1000),i],col="red",xlab="Distr. posteriori Var- caso t-student",main = "")  
     abline(v=3)
   }

###### diagnostic
   for (j in 1:M) {
     T[j,1:9]<-c(mean(B0.n[,j]),
                 mean(B1.n[,j]),
                 mean(var.n[-c(1:1000),j]),
                 quantile(B0.n[,j],0.025),
                 quantile(B0.n[,j],0.975),
                 quantile(B1.n[,j],0.025),
                 quantile(B1.n[,j],0.975),
                 quantile(var.n[,j],0.025),
                 quantile(var.n[,j],0.975))
     
     T[j,10:18]<-c(mean(B0[,j]),
                 mean(B1[,j]),
                 mean(var[-c(1:1000),j]),
                 quantile(B0[,j],0.025),
                 quantile(B0[,j],0.975),
                 quantile(B1[,j],0.025),
                 quantile(B1[,j],0.975),
                 quantile(var[,j],0.025),
                 quantile(var[,j],0.975))
}
T<-as.data.frame(T)
#estimadores
mean(T$b0.N)
mean((T$b0.N-b0)^2)
mean(T$b0.T)
mean((T$b0.T-b0)^2)

mean(T$b1.N)
mean((T$b1.N-b1)^2)
mean(T$b1.T)
mean((T$b1.T-b1)^2)

mean(T$Var.N)
mean((T$Var.N-3)^2)
mean(T$Var.T)
mean((T$Var.T-3)^2)


mean(T$Li.b0.N)
mean(T$Ls.b0.N)

mean(T$Li.b1.N)
mean(T$Ls.b1.N)


## cobertura
mean((b0>T$Li.b0.N)&(b0<T$Ls.b0.N))
mean((b1>T$Li.b1.N)&(b1<T$Ls.b1.N))
mean((3>T$Li.Var.N)&(3<T$Li.Var.N))

mean((b0>T$Li.b0.T)&(b0<T$Ls.b0.T))
mean((b1>T$Li.b1.T)&(b1<T$Ls.b1.T))
mean((3>T$Li.Var.T)&(3<T$Li.Var.T))

# AMPLITUD
mean(T$Ls.b0.N-T$Li.b0.N)
mean(T$Ls.b1.N-T$Li.b1.N)
mean(T$Ls.Var.N-T$Li.Var.N)

mean(T$Ls.b0.T-T$Li.b0.T)
mean(T$Ls.b1.T-T$Li.b1.T)
mean(T$Ls.Var.T-T$Li.Var.T)


