#' The function for adaptive sum of powered score (aSPU) test 
#' 
#' It gives the simulation-based p-values of the sum of powered score (SPU) with different power and aSPU test.
#' @param betas the vector (length = p) of maximum likelihood estimates; p is the number of paramters to be tested.
#' @param vcov_mat the p by p variance-covariance matrix of MLE.
#' @param pow the power vector for aSPU tests.
#' @param B number of simulations 
#' @return  simulation-based p-values of the sum of powered score (SPU) with different power and aSPU test
#' @details the suggested power term is c(1:8, Inf)
#' @author Zhiyuan (Jason) Xu, Yiwei Zhang and Wei Pan
#' @import MASS lme4 
#' @references 
#' Wei Pan, Junghi Kim, Yiwei Zhang, Xiaotong Shen and Peng Wei (2014) A powerful and adaptive
#' association test for rare variants, Genetics, 197(4), 1081-95
#' 
#' Yiwei Zhang, Zhiyuan Xu, Xiaotong Shen, Wei Pan (2014) Testing for association with multiple
#' traits in generalized estimation equations, with application to neuroimaging data. Neuroimage.
#' 96:309-25
#' @examples 
#' # install.packages(lme4)
#' library(lme4)
#' library(MASS)
#' data(exdat_GLMM)
#' fit1<-glmer(Y~.-ID + (1|ID), ,family="binomial",data=dat,control=glmerControl(optimizer="bobyqa"))
#' betas<-fixef(fit1)[-1]
#' vcov_mat<-vcov(fit1)[-1,-1]
#' approx_aSPU(betas=betas,vcov_mat = vcov_mat,pow=c(1:8,Inf),B=1000) 
#' @export
approx_aSPU<-function(betas,vcov_mat,pow = c(1:8,Inf),B){
  ## SPU.mulT is the function to calcualte observed SPU statistic and do permutation and calculate p-values for SPU test
  ## This is to implement modified version of SPU test (created by Yiwei Zhang, 
  ##                                                    revised by Zhiyuan Xu 12/2014)
  ## U is arranged as (u_11,...u_i1,...u_1k,...u_ik), where u_ij, i is for SNPs and j is for traits 
  ## SPU=[sum_i (sum_j U_ij^gamma1)^gamma2/gamma1]^{1/gamma2}
  
  #### INPUT:
  ## U: the score vector
  ## V: covariance matrix of the score vector
  ## gamma1: the power parameters for SNPs
  ## gamma2: the power parameters for traits
  ## B: permutation number
  ## K: number of traits
  ## weight: whether wSPU should be outputed
  #### SPU.mulT OUTPUT:
  ## p-values for SPU test
  ####### observed SPU 
  vcov_mat = as.matrix(vcov_mat)
  betas = as.vector(betas)
  V = ginv(vcov_mat)
  U = CovU %*% betas
  K = length(betas)
  gamma2 = pow
  gamma1 = 1
  weight = F
  spu=spuval(U=U,V=V, gamma1=gamma1,gamma2=gamma2,K,weight=weight)
  p=rep(NA,length(gamma1)*length(gamma2))
  ###### NULL Distribution
  T0s=matrix(NA,nrow=B,ncol=length(gamma1)*length(gamma2))
  u.null<-mvrnorm(B,rep(0,length(U)),V)
  for (b in 1:B) T0s[b,]=spuval(U=u.null[b,],V=V,gamma1=gamma1,gamma2=gamma2,K=K,weight=weight)
  for (g in 1:ncol(T0s))  p[g]<-sum(abs(spu[g])<abs(T0s[,g]))/B
  
  ##### p-values
  P0s<-PermPvs(T0s,B=B)
  minp0<-apply(P0s,1,min)
  Paspu<-sum(minp0<min(p))/B
  pvs<-data.frame(matrix(c(p,Paspu),nrow=1))
  
  names(pvs) = c(paste0("SPU_",pow),"aSPU")
  return(pvs)
}

## the SPU statistic is matrix, rows for gamma2 and columns for gamma1, where gamma2 is for traits and gamma1 is for SNPs

## PermPvs is the function to calcualte p-values for SPU test by permutation
PermPvs<-function(T0s,B){
  #P0s<-matrix(NA,ncol=ncol(T0s),nrow=B)
  #for (g in 1:ncol(T0s))
  #  for (b in 1:B)
  #      P0s[b,g]=sum(abs(T0s[b,g])<abs(T0s[-b,g]) )/(B-1)
  P0s=apply(T0s,2,function(z)(B-rank(abs(z)))/(B-1)   ) 
  return(P0s)
}


## spuval is the function to calculate SPU statistic values
spuval<-function(U,V,gamma1,gamma2, K,weight=F){
  if (weight==F) U2=matrix(U,nrow=K,byrow=T)  ## rows for traits and columns for snps
  if (weight==T) U2=matrix(U/sqrt(diag(V)),nrow=K,byrow=T) 
  spumat<-matrix(NA,ncol=length(gamma1),nrow=length(gamma2))
  for (g1 in 1:length(gamma1)){
    if (gamma1[g1]<Inf) spu.g1<-apply(U2,1,function(z)sum(z^gamma1[g1])) else spu.g1<-apply(U2,1,function(z)max(abs(z)))
    for (g2 in 1:length(gamma2)){
      if (gamma2[g2]<Inf) spumat[g2,g1]=sum(spu.g1^(gamma2[g2]/gamma1[g1])) else spumat[g2,g1]=max(abs(spu.g1))
    }
  }
  spu0=c(t(spumat))    
  spu0
}
