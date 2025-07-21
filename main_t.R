library(foreach)
library(doParallel)
library(parallel)
library(distr, exclude = "rbind")
library(distrEllipse, exclude = "rbind")
library(pROC)  #AUC

set.seed(123)
numCores <- detectCores()
cl <- makeCluster(numCores-1)
registerDoParallel(cl)

AUC_norm <- function(omega,estimator){
  omega[which(omega!=0)]=1
  estimator[which(estimator!=0)]=1
  true_labels <- as.vector(omega)
  predicted_probs <- as.vector(estimator)
  
  # ROC
  roc_obj <- tryCatch({
    roc(true_labels, predicted_probs)
  }, error = function(err) {
    cat("Error in roc calculation:", conditionMessage(err), "\n")
    NULL
  })
  
  # 
  if (!is.null(roc_obj)) {
    # AUC
    auc_value <- auc(roc_obj)
    #print("AUC:", auc_value)
    cat("AUC:", auc_value, "\n")
    
  } else {
    #
    auc_value <- 0
    #print("AUC:", auc_value)
    cat("AUC:", auc_value, "\n")
  }
  
  
  
  return(auc_value)
  
}

F_norm <- function(omega,estimator){
  p <- ncol(omega)
  x <- abs(omega-estimator)
  loss <- sqrt(sum(x^2))/p
  return(loss)
}

te.trans.clime.new <- function(dat.all,const,A.size,epsilon_0=0.01){
  library(mvtnorm)
  library(fastclime)
  library(Matrix)
  library(pcaPP)
  
  ### rank-based scaled clime
  MyRSclime <- function(X, lambda=0.1){
    p <- ncol(X)
    obj <- rep(-1, 2*p)
    S.hat <- cor.fk(X)
    C <- matrix(lambda, nrow=p, ncol=p)
    R <- S.hat + C
    W <- S.hat - C
    mat <- base::rbind(cbind(W,-R), 
                       cbind(-R,W))
    feasible <- TRUE
    Theta.hat <- matrix(ncol=p,nrow=p)
    for(j in 1:p){
      ej <- rep(0, 2*p)
      ej[j] <- 1
      ej[p+j] <- -1
      rhs <- ej
      out.txt <- capture.output(betapm <- fastlp(obj=obj, mat=mat, rhs=rhs))
      if(!grepl("optimal", out.txt) ){
        feasible=FALSE
        break
      }
      Theta.hat[,j] <- betapm[1:p] - betapm[(p+1):(2*p)]
    }
    if(!feasible){
      cat('Theta.hat not found','\n')
    }
    list(Theta.hat=Theta.hat, conv=feasible)
  }
  
  ###algorithm based on fastclime package
  ### min \|Theta\|_1
  ###subject to \|cov(X)%*%Theta-Bmat\|_max <=lambda (if X is the raw sample)
  ###subject to \|X%*%Theta-Bmat\|_max <=lambda (if X is the sample covariance matrix)
  
  Myfastclime.s<-function(X,Bmat,lambda=0.1, scale=T, n){
    p<-ncol(X)
    obj=rep(-1,2*p)
    obj_bar=rep(0,2*p)
    rhs_bar<-rep(1, 2*p)
    if(isSymmetric(X, tol=10^(-4))){
      Sig.hat<-X
    }else{
      Sig.hat<-cov(X)
    }
    
    Sig.hat0<-Sig.hat
    feasible=T
    Theta.hat<-NULL
    Sig.diag<-diag(Sig.hat)
    Sig.hat<-cov2cor(Sig.hat)
    mat=base::rbind(cbind(Sig.hat,-Sig.hat),
                    cbind(-Sig.hat,Sig.hat))
    for(j in 1:p){
      rhs <- c(Bmat[,j],-Bmat[,j])
      out.txt<-capture.output(  fastlp.re<-fastlp(obj=obj, mat=mat, rhs=rhs+rhs_bar*lambda))
      if(!grepl("optimal", out.txt) ){
        feasible=F
        break
      }
      Theta.hat<-cbind(Theta.hat,(fastlp.re[1:p]-fastlp.re[-(1:p)])/sqrt(Sig.diag[j])/sqrt(Sig.diag))
      if(scale & sum(Theta.hat[,j]==0)==p){
        feasible=F
        break
      }else if(scale){
        Theta.hat[,j]<-  as.numeric(Bmat[j,j]/ (Sig.hat0[j,]%*%Theta.hat[,j]))*Theta.hat[,j]

      }
    }
    if(!feasible){
      cat('Theta.hat not found','\n')
      Theta.hat<-solve(cov(X)+diag(lambda,p))%*%Bmat
    }
    list(Theta.hat=Theta.hat, conv=feasible)
  }
  
  ####the main Trans-CLIME algorithm###
  ##X: primary data; X.A: {X^{(k)}, k in A}; lambda:lambda.Theta; 
  ### agg: perform LS aggregation or not; X.til: samples for aggregation; Theta.cl: CLIME estimator
  Trans.CLIME<-function(X,X.A, const, agg=T, X.til=NULL,Theta.cl){
    if(agg &is.null(X.til)){
      cat('no aggregation samples provided.','\n')
    }
    n0=nrow(X)
    nA<-nrow(X.A)
    p<-ncol(X)
    
    lam.delta<-2*sqrt(log(p)/n0) 
    omega.l1<-mean(apply(Theta.cl,2, function(x) sum(abs(x))))
    Delta.re <- Myfastclime.s(X=diag(1,p), Bmat=diag(1,p)-t(Theta.cl)%*%cor.fk(X.A), 
                              lambda=omega.l1*sqrt(log(p)/n0) , scale=F)
    if(Delta.re$conv){Delta.hat<-Delta.re$Theta.hat}else{ Delta.hat<-diag(0,p) }
    Theta.re <- Myfastclime.s(X=cor.fk(X.A), Bmat=diag(1,p)-t(Delta.hat),
                              lambda=2*const*sqrt(log(p)/nA))
    Theta.hat<-Theta.re$Theta.hat
    
    if(agg){
      Omega.hat<-Agg(Theta.init=cbind(Theta.cl, Theta.hat), X.til=X.til)
    }else{
      Omega.hat<-Theta.hat
    }
    Omega.hat
  }
  
  ####LS aggregation function with 
  ### Theta.init=(Omega.clime, Theta.hat) and X.til: some primary samples
  Agg<-function(Theta.init, X.til){
    p<-ncol(X.til)
    n.til<-nrow(X.til)
    v.mat<-sapply(1:p, function(j){
      W.j<-cov(X.til%*%cbind(Theta.init[,j], Theta.init[,p+j]))
      v0=rep(0,2)
      v0[which.min(c(W.j[1,1]-2*Theta.init[j,j], W.j[2,2]-2*Theta.init[j,p+j]))]<-1
      v0
    })
    Theta.hat<-sapply(1:p, function(j) cbind(Theta.init[,j], Theta.init[,p+j])%*% v.mat[,j])
    
    Theta.hat
  }
  
  
  epsilon_0 <- epsilon_0
  
  Theta0<-dat.all$Theta0 ##true Omega
  X.all<-dat.all$X
  
  X0 <- X.all[1:n.vec[1],]
  n00 <- n.vec[1]
  
  ################## split
  # 
  indices <- sample(1:n00)
  
  # 
  half_size <- n00 / 2
  T0_indices <- indices[1:half_size]
  V0_indices <- indices[(half_size + 1):n00]
  
  # 
  T00 <- X0[T0_indices, ]
  V00 <- X0[V0_indices, ]
  
  ####################
  
  
  
  Gamma.re0 <- MyRSclime(X=T00,lambda=const*2*sqrt(log(p)/half_size))
  Gamma.hat0 <- Gamma.re0$Theta.hat
  Gamma.hat0 <- (Gamma.hat0+t(Gamma.hat0))/2
  
  
  if(isSymmetric(V00, tol=10^(-4))){
    Sig.hat.new<-V00
  }else{
    Sig.hat.new<-cov(V00)
  }
  Sig.hat00<-cov2cor(Sig.hat.new)
  
  Loss_Omega <- function(Omega_loss){
    #  trace(Omega %*% S0)
    trace_product <- sum(diag(Omega_loss %*% Sig.hat00))
    
    # log(det(Omega))
    log_det_Omega <- log(det(Omega_loss))
    
    # 
    result <- trace_product - log_det_Omega
    
    return(result)
  }
  
  Loss_Gamma0 <- Loss_Omega(Gamma.hat0)
  
  
  n0=round(n.vec[1]*4/5) #split sample for aggregation
  
  
  Theta.re01 <- MyRSclime(X=X.all[1:n0,], lambda=2*const*sqrt(log(p)/n0))
  Theta.init1<-Theta.re01$Theta.hat
  
  ind2<-(n.vec[1]-n0+1): n.vec[1]
  Theta.re02<-MyRSclime(X=X.all[ind2,], lambda=2*const*sqrt(log(p)/n0))
  Theta.init2<-Theta.re02$Theta.hat
  
  k_set <- NULL
  
  for(k in 1:K){
    
    
    # trans-te-clime
    
    
    k.seq=(sum(n.vec[1:k])+1):sum(n.vec[1:(k+1)])
    
    Omega.tl1 <- Trans.CLIME(X=X.all[1:n0,], X.A=X.all[k.seq,], const=const,
                             X.til= X.all[(n0+1):n.vec[1],], Theta.cl=Theta.init1)
    Omega.tl2 <- Trans.CLIME(X=X.all[ind2,], X.A=X.all[k.seq,], const=const,
                             X.til= X.all[1:(n.vec[1]-n0),], Theta.cl=Theta.init2)
    Omega.tl<-(Omega.tl1+Omega.tl2)/2
    Omega.tl.s <- (Omega.tl+t(Omega.tl))/2
    
    Loss_Omegak <- Loss_Omega(Omega.tl.s)
    
    if(Loss_Omegak<=(1+epsilon_0)*Loss_Gamma0){
      k_set <- c(k_set,k)
    }else{
      k_set <- k_set
    }
    
  }
  
  if(is.null(k_set)){
    k_set <- 1:4
  }else{
    k_set <- k_set
  }
  
  A0 <- 1:A.size
  
  
  # 
  all_set <- 1:K
  
  # Complement of A0
  A0_complement <- setdiff(all_set, A0)
  
  # Type I Error
  false_positives <- setdiff(k_set, A0)
  type1_error <- length(false_positives) / length(A0_complement)
  
  # Power
  true_positives <- intersect(k_set, A0)
  power <- length(true_positives) / length(A0)
  
  return(list(type1_error=type1_error, power=power))
  
  
  
}

ThetaGen <- function(p, type='Toep'){
  if(type=='Toep'){
    Theta<-toeplitz(0.5^(1:p)*2)
    Theta[which(abs(Theta)<=0.05, arr.ind=T)]<- 0
  }else if(type=='Bdiag'){
    Theta<-kronecker(diag(p/5), toeplitz(c(1.5,1.2,0.9,0.6,0.3)))
  }else if(type=='trand'){
    B<-matrix(rbinom(p*p,1,0.01)/2,ncol=p,nrow=p)
    diag(B)<-0
    B[lower.tri(B)]<-0
    B<-B+t(B)
    delta.list<-seq(1,5,by=0.1)
    cn.list<-NULL
    for(delta in delta.list){
      cn.list<-c(cn.list, kappa(B+delta*diag(p)))
    }
    delta<-delta.list[which.max(cn.list>p)]
    Theta<-B+delta*diag(p)
  }else if(type=='graph'){
    Theta<-matrix(nrow=p,ncol=p)
    Y1 <- runif(p)
    Y2 <- runif(p)
    for(i in 2:p){
      for(j in 1:(i-1)){
        pr<-exp(-((Y1[i]-Y1[j])^2+(Y2[i]-Y2[j])^2)/0.25)/sqrt(2*pi)
        Theta[i,j] <- rbinom(1,1,pr)*0.145
      }
    }
    diag(Theta)<-1
  }
  Sig.diag <- diag(solve(Theta))
  Theta <- diag(sqrt(Sig.diag)) %*% Theta %*% diag(sqrt(Sig.diag))
  return(Theta)
}

DataGen<-function(K,A.size,h, n.vec, distr='norm', add.para=NULL, 
                  add.func1=NULL, add.func2=NULL,add.func3=NULL,
                  p=100,Theta=NULL,seed=NULL){
  if(!is.null(seed)){
    set.seed(seed)
  }
  distroptions("WarningArith" = FALSE)
  if(distr == 'norm'){
    rdist <- rmvnorm
  }
  if(distr == 't'){
    rdist <- function(n, mu, sigma){
      rmvt(n, sigma, 3.5, mu)
    }
  }
  if(distr == 'ell'){
    rdist <- function(n, mu, sigma){
      #A <- t(chol(sigma, pivot=TRUE))
      edistr <- EllipticalDistribution(radDistr = add.func1, loc = mu, 
                                       scale = sigma)
      t(edistr@r(n))
    }
  }
  if(distr == 'te'){
    rdist <- function(n, mu, sigma){
      #A <- t(chol(sigma, pivot=TRUE))
      edistr <- EllipticalDistribution(radDistr = add.func3, loc = mu, 
                                       scale = sigma)
      x <- t(edistr@r(n))
      t(apply(x, 1, add.func2))
    }
  }
  if(distr == 'npnorm'){
    rdist <- function(n, mu, sigma){
      x <- rmvnorm(n, mu, sigma)
      t(apply(x,1, add.func2))
    }
  }
  if(distr == 'npn-exp'){
    rdist <- function(n, mu, sigma){
      x <- rmvnorm(n, mu, sigma)
      t(apply(x, 1, function(t) (exp(t)-exp(0.5))/((exp(1)-1)*exp(1))))
    }
  }
  if(distr == 't-1'){
    rdist <- function(n, mu, sigma){
      edistr <- EllipticalDistribution(radDistr = add.func3, loc = mu, 
                                       scale = sigma)
      x <- t(edistr@r(n))
    }
  }
  
  
  Sig<- solve(Theta)
  X<-rdist(n.vec[1],rep(0,p),Sig)
  Theta.out<-diag(1,p)
  Omega.vec<-0
  for(k in 1 : K){
    if(k<=A.size){
      Delta.k<-matrix(rbinom(p^2,size=1,prob=0.1)*runif(p^2,-h/p,h/p),ncol=p)
      #  cat(max(colSums(abs(Delta.k))),'\n')
      Sig.k<-(diag(1,p)+Delta.k)%*%Sig 
      Sig.k<-(Sig.k+t(Sig.k))/2
      if(min(eigen(Sig.k)$values)<0.05){
        Sig.k<-Sig.k+diag(0.1-min(eigen(Sig.k)$values),p)
      }
      X<- base::rbind(X, rdist(n.vec[k+1],rep(0,p),Sig.k))
    }else{
      #Delta.out<-diag(0.5,p)+matrix(rbinom(p^2,size=1,prob=0.4)*runif(p^2,-0.2,0.2),ncol=p)
      Delta.out<-diag(1,p)+matrix(rbinom(p^2,size=1,prob=0.1)*0.4,ncol=p)
      Sig.k<-(diag(1,p)+Delta.out)%*%Sig
      Sig.k<-(Sig.k+t(Sig.k))/2
      if(min(eigen(Sig.k)$values)<0.05){
        Sig.k<-Sig.k+diag(0.1-min(eigen(Sig.k)$values),p)
      }
      X<- base::rbind(X, rdist(n.vec[k+1],rep(0,p),Sig.k))
    }
    Omega.vec<-c(Omega.vec, max(colSums(abs(diag(1,p)-Sig.k%*%Theta))))
    
  }
  #cat(Omega.vec,'\n')
  list(X=X, Theta0=Theta)
  
}

h.1 <- function(x){
  x.r <- numeric(length(x))
  in.m <- 1:(length(x)/5)
  x.r[5*in.m-4] <- x[5*in.m-4]
  x.r[5*in.m-3] <- sign(x[5*in.m-3])*sqrt(abs(x[5*in.m-3]))/sqrt(2/pi)
  x.r[5*in.m-2] <- x[5*in.m-2]^3/sqrt(15)
  x.r[5*in.m-1] <- (pnorm(x[5*in.m-1])-0.5)/(1/12)
  x.r[5*in.m] <- (exp(x[5*in.m])-exp(0.5))/((exp(1)-1)*exp(1))
  return(x.r)
}

##########################################
Theta0.type.list=c("Toep", "Bdiag") # or Bdiag or Rand or Toep or trand pr graph
K=5
n.vec<-c(100,rep(150,K))
A.size.list=c(4)   #c(2,3,4)
h.list <- c(10,20,30)

epsilon_0.list <- c(0.001,0.005,0.01,0.02,0.05,0.1,0.2)


# repeat parameters
models <- c("t")
#models <- c("te")
add.func2 <- h.1
p.list <- c(50,100,150)
rep.times <- 100
#rep.times <- 2
numCores <- detectCores()
const.list <- c(seq(0.001,0.009,by=0.001),
                seq(0.01,0.09,by=0.01), 
                seq(0.1, 2, by=0.05))

##########################################
parafunc <- function(s,h){
  dat.all<-DataGen(K=5, A.size=A.size, p=p, h=h, n.vec=n.vec,
                   distr=model ,add.func1=add.func1,
                   add.func2=h.1,add.func3 = add.func3,Theta=Theta0,seed=s)
  t1 <- tryCatch(te.trans.clime.new(dat.all,const,A.size,epsilon_0),error=function(err) return(NULL))
  
  
  return(c(t1$type1_error,rep(0,5),
           t1$power,rep(0,5)))
  
  
}



total.res <- list()
for(model in models){
  for(h in h.list){
    for(A.size in A.size.list){
      for(Theta0.type in Theta0.type.list){
        for(p in p.list){
          for(epsilon_0 in epsilon_0.list){
          add.func1 <- sqrt(Chisq(p)/Chisq(1))
          add.func3 <- Fd(p,1)
          Theta0 <- ThetaGen(p,type=Theta0.type)
          all.res <- NULL
          for(const in const.list){
            
            clusterExport(cl, varlist=c("add.func2","h.1","add.func3","add.func1",
                                        "te.trans.clime.new","Theta0",
                                        "DataGen","A.size","p","n.vec","const"))
            res<-foreach(s = 1:rep.times, .combine = base::rbind,
                         .packages = c("mvtnorm","fastclime","Matrix","pROC",
                                       "pcaPP","distr","distrEllipse")) %dopar% parafunc(s,h)
            
            res.f <- c(apply(res, 2, mean), apply(res, 2, sd))
            all.res <- base::rbind(all.res, res.f)
            print(c(model, A.size, Theta0.type, p,epsilon_0, const))
          }
          total.res[[paste0(model,".h",h,".",A.size,".",Theta0.type,".",p,".",
                            epsilon_0)]] <- all.res
          }
        }
      }
    }
  }
}



best.res<-matrix(ncol=24,nrow=length(total.res))
for(i in 1:length(total.res)){
  value1 <- total.res[[i]][,1:6]
  bestvalue <- numeric(24)
  for(j in 1:6){
    bv.in <- which.max(value1[,j])
    rel.in <- c(j, 6+j, 12+j, 18+j)
    bestvalue[rel.in] <- total.res[[i]][bv.in, rel.in]
  }
  best.res[i,]<-bestvalue
}

write.csv(best.res, paste(rep.times,models, "best.res.csv", sep = "_"), 
          row.names = names(total.res))

stopCluster(cl)

set.seed(123)








