


Trans.CLIME.rd <- function(X.all,const){
  library(mvtnorm)
  library(fastclime)
  library(Matrix)
  
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
        # Theta.hat[,j]<-as.numeric(Theta.hat[j,j]/ (t(Theta.hat[,j])%*%Sig.hat0%*%Theta.hat[,j]))*Theta.hat[,j]
        
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
    sigA.hat<-mean(apply(X.A, 2, sd))
    sig0.hat<-mean(apply(X, 2, sd))
    
    lam.delta<-2*sig0.hat*sqrt(log(p)/n0) 
    omega.l1<-mean(apply(Theta.cl,2, function(x) sum(abs(x))))
    Delta.re <- Myfastclime.s(X=diag(1,p), Bmat=diag(1,p)-t(Theta.cl)%*%cov(X.A), lambda=omega.l1*sqrt(log(p)/n0) , scale=F)
    if(Delta.re$conv){Delta.hat<-Delta.re$Theta.hat}else{ Delta.hat<-diag(0,p) }
    Theta.re <- Myfastclime.s(X=cov(X.A), Bmat=diag(1,p)-t(Delta.hat),
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
  
  
  ###original clime####
  Theta.re0<-Myfastclime.s(X=X.all[1:n.vec[1],], Bmat=diag(1,p), 
                           lambda=const*2*sqrt(log(p)/n.vec[1]))
  Theta.hat0<-Theta.re0$Theta.hat
  Theta.hat0.s <- (Theta.hat0+t(Theta.hat0))/2
  
  #Trans-CLIME
  n0=round(n.vec[1]*4/5) #split sample for aggregation
  Theta.re0<-Myfastclime.s(X=X.all[1:n0,], Bmat=diag(1,p), lambda=const*2*sqrt(log(p)/n0))
  Theta.init<-Theta.re0$Theta.hat
  Omega.tl1 <- Trans.CLIME(X=X.all[1:n0,], X.A=X.all[-(1:n.vec[1]),], const=const,
                           X.til= X.all[(n0+1):n.vec[1],], Theta.cl=Theta.init)
  ind2<-(n.vec[1]-n0+1): n.vec[1]
  Theta.re0<-Myfastclime.s(X=X.all[ind2,], Bmat=diag(1,p), lambda=const*2*sqrt(log(p)/n0))
  Theta.init<-Theta.re0$Theta.hat
  Omega.tl2 <- Trans.CLIME(X=X.all[ind2,], X.A=X.all[-(1:n.vec[1]),], const=const,
                           X.til= X.all[1:(n.vec[1]-n0),], Theta.cl=Theta.init)
  Omega.tl<-(Omega.tl1+Omega.tl2)/2
  Omega.tl.s <- (Omega.tl+t(Omega.tl))/2
  
  return(list(Omega.ori=Theta.hat0.s, Omega=Omega.tl.s))
}



te.trans.clime.rd <- function(X.all,n.vec,const){
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
        # Theta.hat[,j]<-as.numeric(Theta.hat[j,j]/ (t(Theta.hat[,j])%*%Sig.hat0%*%Theta.hat[,j]))*Theta.hat[,j]
        
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
    #sigA.hat<-mean(apply(X.A, 2, sd))
    #sig0.hat<-mean(apply(X, 2, sd))
    
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
  
  
  Theta.re0 <- MyRSclime(X=X.all[1:n.vec[1],],lambda=const*2*sqrt(log(p)/n.vec[1]))
  Theta.hat0<-Theta.re0$Theta.hat
  Theta.hat0.s <- (Theta.hat0+t(Theta.hat0))/2
  
  if(length(n.vec)>1){
    # trans-te-clime
    
    #Theta.re0<-Myfastclime.s(X=X.all[1:n.vec[1],], Bmat=diag(1,p), lambda=2*sqrt(log(p)/n0), isrank=TRUE)
    n0=round(n.vec[1]*4/5)
    Theta.re01 <- MyRSclime(X=X.all[1:n0,], lambda=2*const*sqrt(log(p)/n0))
    Theta.init1<-Theta.re01$Theta.hat
    ind2<-(n.vec[1]-n0+1): n.vec[1]
    Theta.re02<-MyRSclime(X=X.all[ind2,], lambda=2*const*sqrt(log(p)/n0))
    Theta.init2<-Theta.re02$Theta.hat
    Omega.tl1 <- Trans.CLIME(X=X.all[1:n0,], X.A=X.all[-(1:n.vec[1]),], const=const,
                             X.til= X.all[(n0+1):n.vec[1],], Theta.cl=Theta.init1)
    Omega.tl2 <- Trans.CLIME(X=X.all[ind2,], X.A=X.all[-(1:n.vec[1]),], const=const,
                             X.til= X.all[1:(n.vec[1]-n0),], Theta.cl=Theta.init2)
    Omega.tl<-(Omega.tl1+Omega.tl2)/2
    Omega.tl.s <- (Omega.tl+t(Omega.tl))/2
  }else{
    Omega.tl.s <- Theta.hat0.s
  }
  
  
  return(list(Omega.ori=Theta.hat0.s, Omega=Omega.tl.s))
  
}



te.trans.clime.rd.new <- function(X.all,n.vec,const){
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
        # Theta.hat[,j]<-as.numeric(Theta.hat[j,j]/ (t(Theta.hat[,j])%*%Sig.hat0%*%Theta.hat[,j]))*Theta.hat[,j]
        
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
    #sigA.hat<-mean(apply(X.A, 2, sd))
    #sig0.hat<-mean(apply(X, 2, sd))
    
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
  
  
  Theta.re0 <- MyRSclime(X=X.all[1:n.vec[1],],lambda=const*2*sqrt(log(p)/n.vec[1]))
  Theta.hat0<-Theta.re0$Theta.hat
  Theta.hat0.s <- (Theta.hat0+t(Theta.hat0))/2
  
  
  
  if(length(n.vec)>1){
    # trans-te-clime
    
    epsilon_0 <- 0.01
    
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
      # trace(Omega %*% S0)
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
    
    K <- length(n.vec)-1
    
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
      k_set <- 1:K
    }else{
      k_set <- k_set
    }
    
    k.seq2 <- NULL
    
    for(k in k_set){
      k.seq=(sum(n.vec[1:k])+1):sum(n.vec[1:(k+1)])
      k.seq2 <- c(k.seq2,k.seq)
    }
    
    Omega.tl1 <- Trans.CLIME(X=X.all[1:n0,], X.A=X.all[k.seq2,], const=const,
                             X.til= X.all[(n0+1):n.vec[1],], Theta.cl=Theta.init1)
    Omega.tl2 <- Trans.CLIME(X=X.all[ind2,], X.A=X.all[k.seq2,], const=const,
                             X.til= X.all[1:(n.vec[1]-n0),], Theta.cl=Theta.init2)
    Omega.tl<-(Omega.tl1+Omega.tl2)/2
    Omega.tl.s00 <- (Omega.tl+t(Omega.tl))/2
    
  }else{
    Omega.tl.s00 <- Theta.hat0.s
  }
  
  
  
  return(list(Omega.ori=Theta.hat0.s, Omega=Omega.tl.s00))
  
}













