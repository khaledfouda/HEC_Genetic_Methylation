##########################
#polar decomposition
##########################
PD=function(A){
  X=array(0,c(dim(A),100))
  Y=array(0,c(dim(A),100))
  X[,,1]=A
  k=0
  ratio=1
  alpha=rep(0,100)
  beta=rep(0,100)
  gamma=rep(0,100)
  while(ratio>1e-5){
    k=k+1
    Y[,,k]=solve(X[,,k])
    if (ratio<1e-4) {
      gamma[k]=1
    } else{
      alpha[k]=sqrt(norm(X[,,k],type="1")*norm(X[,,k],type="i"))
      beta[k]=sqrt(norm(Y[,,k],type="1")*norm(Y[,,k],type="i"))
      gamma[k]=sqrt(beta[k]/alpha[k])
    }
    X[,,k+1]=(gamma[k]*X[,,k]+t(Y[,,k])/gamma[k])/2
    ratio=norm(X[,,k+1]-X[,,k],type="1")/norm(X[,,k],type="1")
  }
  U=X[,,k+1]
  H1=t(U)%*%A
  H=(H1+t(H1))/2
  
  return(list(U=U,H=H))
}


##########################
#Fast Randomized SVT
##########################
FRSVT=function(A,tau,k,p,ita){
  n1=dim(A)[1]
  n2=dim(A)[2]
  l=k+p
  
  Omega=matrix(rnorm(n2*l,0,1),nrow=n2)
  Y=A%*%Omega
  
  qrQ=qr(Y)
  Q=qr.Q(qrQ)[,1:qrQ$rank]
  
  
  for (i in 1:ita){
    Q=qr.Q(qr(A%*%t(A)%*%Q))
  }
  
  B=t(Q)%*%A
  
  HC=qr(t(B))
  H=qr.Q(HC)
  C=qr.R(HC)
  
  PDC=PD(C)
  
  W=PDC$U
  P=PDC$H
  
  EDP=eigen(P)
  V=EDP$vectors
  D=EDP$values
  
  SlambdaAu=Q%*%V
  SlambdaAd=pmax(D-tau,0)
  SlambdaAv=H%*%W%*%V
  
  return(list(u=SlambdaAu,d=SlambdaAd,v=SlambdaAv))
}




##########################
#Validation error
##########################
cv_validate_sim=function(gridlength,xgridlength,svdY,alpha,uvalidatemat,
                         Xbeta_validate,umat_validate){
  
  cv_validate=array(0,c(gridlength,xgridlength))
  
  for (i in 1:gridlength){
    
    B_validate=as.vector(SVTE_alpha(svdY$u,svdY$d,svdY$v,svdY$d[i],alpha))[uvalidatemat!=0]
    
    for (j in 1:xgridlength){
      
      cv_validate[i,j]=sum((Xbeta_validate[j,]+B_validate-umat_validate)^2)/sqrt(sum(uvalidatemat!=0))
      
    }
  }
  
  return(cv_validate)
}



##########################
#Scale SVT
##########################
SVTE_alpha=function(svdu,svdd,svdv,lambda,alpha){
  return(svdu%*%(pmax(svdd-lambda*alpha,0)*t(svdv))/(1+2*lambda*(1-alpha)))
}




########################################
#Proposed Method
########################################

MCCIfit=function(Aobs,X,nfolds=5,SVT_type="exact",lambda1_grid=seq(0.9,0.1,length=30),
                 lambda2_gridlength=150,alpha_grid=seq(0.992,1,length=20)){
  n1=dim(Aobs)[1]
  n2=dim(Aobs)[2]
  m=dim(X)[2]
  
  Xtheta=X[,1:3]
  
  PX=X%*%solve(t(X)%*%X)%*%t(X)
  Eye=diag(1,n1)
  PXp=Eye-PX
  
  iomega=as.numeric(Aobs!=0)
  omega=matrix(iomega,n1,n2)
  iobs=which(iomega==1)
  imiss=which(iomega==0)
  
  gammahat=matrix(rep(NA,n1*n2),nrow=n1)
  for (i in 1:n2){
    thetadata=data.frame(cbind(omega[,i],Xtheta))
    thetaglmout=glm(X1~.,family=binomial(logit),data=thetadata)
    gammahat[,i]=predict(thetaglmout,type = "response")
  }
  
  Aobsw=Aobs/gammahat
  
  PXpAobsw=PXp%*%Aobsw
  
  svdPXpAobsw=switch(SVT_type,
                     "exact"=svd(PXpAobsw),
                     "FRSVT"=FRSVT(PXpAobsw,0,lambda2_gridlength,5,2)
  )
  
  folds <- cut(sample(seq(1,length(iobs))),breaks=nfolds,labels=FALSE)
  iomegatrain=matrix(0,n1*n2,nfolds)
  iomegavalidate=matrix(0,n1*n2,nfolds)
  omegatrain=array(NA,c(n1,n2,nfolds))
  omegavalidate=array(NA,c(n1,n2,nfolds))
  
  gammatrainhat=array(NA,c(n1,n2,nfolds))
  Aobswtrain=array(NA,c(n1,n2,nfolds))
  PXpAobswtrain=array(NA,c(n1,n2,nfolds))
  svdPXpAobswtrain=list()
  A_validate=list()
  Xbetavalidate=list()
  
  for (fold_num in 1:nfolds){
    validateIndexes <- which(folds==fold_num,arr.ind=TRUE)
    ivalidate <- iobs[validateIndexes]
    itrain <- iobs[-validateIndexes]
    
    iomegatrain[itrain,fold_num]=1
    omegatrain[,,fold_num]=matrix(iomegatrain[,fold_num],n1,n2)
    
    iomegavalidate[ivalidate,fold_num]=1
    omegavalidate[,,fold_num]=matrix(iomegavalidate[,fold_num],n1,n2)
    
    for (i in 1:n2){
      thetadata=data.frame(cbind(omegatrain[,i,fold_num],Xtheta))
      thetaglmout=glm(X1~.,family=binomial(logit),data=thetadata)
      gammatrainhat[,i,fold_num]=predict(thetaglmout,type = "response")
    }
    
    Aobswtrain[,,fold_num]=omegatrain[,,fold_num]*Aobs/gammatrainhat[,,fold_num]
    PXpAobswtrain[,,fold_num]=PXp%*%Aobswtrain[,,fold_num]
    svdPXpAobswtrain=append(svdPXpAobswtrain,list(switch(SVT_type,"exact"=svd(PXpAobswtrain[,,fold_num]),
                                                         "FRSVT"=FRSVT(PXpAobswtrain[,,fold_num],0,lambda2_gridlength,5,2)))
    )
    
    A_validate=append(A_validate,list(as.vector(Aobs)[iomegavalidate[,fold_num]!=0]))
    
    Xbetavalidate=append(Xbetavalidate,list(matrix(nrow=length(lambda1_grid),
                                                   ncol = sum(iomegavalidate[,fold_num]))))
  }
  
  for (j in 1:length(lambda1_grid)){
    PX=X%*%solve(t(X)%*%X+lambda1_grid[j]*diag(1,m))%*%t(X)
    
    for (fold_num in 1:nfolds){
      Xbetavalidate[[fold_num]][j,]=as.vector(PX%*%Aobswtrain[,,fold_num])[iomegavalidate[,fold_num]!=0]
    }
  }
            
  cv=array(0,c(nfolds,lambda2_gridlength,length(lambda1_grid),length(alpha_grid)))
  for (fold_num in 1:nfolds){
    for (alpha_num in 1:length(alpha_grid)){
    cv[fold_num,,,alpha_num]=cv_validate_sim(lambda2_gridlength,length(lambda1_grid),svdPXpAobswtrain[[fold_num]],
                                   alpha_grid[alpha_num],iomegavalidate[,fold_num],
                                   Xbetavalidate[[fold_num]],A_validate[[fold_num]])
    }
  }
  
  cv_grid=which(apply(cv,2:4,sum)==min(apply(cv,2:4,sum)),arr.ind=T)[1,]
  
  betahat=solve(t(X)%*%X+lambda1_grid[cv_grid[2]]*diag(1,m))%*%t(X)%*%Aobsw
  Xbetahat=X%*%betahat
  Bhat=SVTE_alpha(svdPXpAobsw$u,svdPXpAobsw$d,svdPXpAobsw$v,
                  svdPXpAobsw$d[cv_grid[1]],alpha_grid[cv_grid[3]])
  rankB=sum(pmax(svdPXpAobsw$d-svdPXpAobsw$d[cv_grid[1]]*alpha_grid[cv_grid[3]],0)>0)
  
  Ahat=Xbetahat+Bhat
  
  return(list(Ahat=Ahat,Bhat=Bhat,betahat=betahat,rank=rankB+m))
}
