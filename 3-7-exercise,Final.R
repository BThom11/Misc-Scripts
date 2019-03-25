## 3/7/19 205a In Class exercise, subliminal priming.
##Likelihood is combination of binomial and normal so the following uses metropolis hastings sampling.
##Hypotheses, if priming effect exists, then intercept (psi), will be below .5, if not it will be =.5
##Most other parameters are not of interest.

#setup
set.seed(132)
library(R2jags)
library(MCMCpack)

#Data. Prac dat is 2 column table of accuracy on 100 trials, and measured effect size of priming.
n=100
Data=read.table("C:/Users/rjtho/Downloads/PracDat.txt",header=T,sep=" ")
y=Data$acc*100
t=Data$eff


#Self Roll, Beginning of actual model####

#likelihood functions. myexp is the mean of the data (mu= beta*exp((theta-psi)/tau)), equivilent to following equation
myExp=function(theta,beta,psi,tau) beta*pexp(theta-psi,rate=1/tau) #function for mu

likeTheta=function(theta,y,t,beta,psi,tau,s2){        #function for likelihood of y and t given all, binomial pdf plus normal pdf
  y*log(theta)+(n-y)*log(1-theta)-(t-myExp(theta,beta,psi,tau))^2/(2*s2)
}
likeBeta=function(theta,t,beta,psi,tau,s2){  #function for LL of beta, second term is pdf of prior on beta (N(50,10)), since beta is independent of y
  -sum((t-myExp(theta,beta,psi,tau))^2/(2*s2))-(beta-50)^2/(2*20^2) }

likePsi=function(theta,t,beta,psi,tau,s2){  #function for LL of beta, second term is pdf of prior on psi (N(.5,.1)), log 
    -sum((t-myExp(theta,beta,psi,tau))^2/(2*s2))-(psi-.5)^2/(2*.1^2) }


#initialize. Just creating vectors to fill, and SD for metropolis hastings. Can tune SDs if acceptance rate too high/low
M=10000
psi=1:M
beta=1:M
theta=matrix(nrow=M,ncol=n)
tau=.06
beta[1]=50
s2=20 ^2  
psi[1]=.5  
theta[1,]=pmax(Data[,1],.5)
counter=rep(0,n)
thetaSD=.08
psiSD=.05
psiCount=0
betaSD=5
betaCount=0

#Sampling with metropolis hastings, could vectorize theta, but fast enough to be irrelevant

for (i in 2:M){
  #for theta
  theta[i,]=theta[i-1,]
  cand=rnorm(n,theta[i-1,],thetaSD) #adding noise to get candidate
  for(j in 1:n){
    if(cand[j]>.5 & cand[j]<1){ #Limit, not accepting candidates outside of these bounds.
      cand.ll=likeTheta(cand[j],y[j],t[j],beta[i-1],psi[i-1],tau,s2) #Comparing LL of candidate to current
      cur.ll=likeTheta(theta[i-1,j],y[j],t[j],beta[i-1],psi[i-1],tau,s2)
      prob=min(1,exp(cand.ll-cur.ll))  #prob is maxed at 1, or equal to difference between cand and cur.
      if(rbinom(1,1,prob)){theta[i,j]=cand[j]; counter[j] = counter[j] +1 } #coin flip to accept candidate or not.
    } 
  }
  
  #for beta
  beta[i]=beta[i-1]
  cand=rnorm(1,beta[i],betaSD)
  cand.ll=likeBeta(theta[i,],t,cand,psi[i-1],tau,s2)
  cur.ll=likeBeta(theta[i,],t,beta[i],psi[i-1],tau,s2)
  prob=min(1,exp(cand.ll-cur.ll))    #prob is maxed at 1, or equal to difference between cand and cur.
  if(rbinom(1,1,prob)){
    beta[i]=cand
    betaCount=betaCount+1}
  
  #for psi
  psi[i]=psi[i-1]
  cand=rnorm(1,psi[i],psiSD)
  if(cand>0){
  cand.ll=likePsi(theta[i,],t,beta[i],cand,tau,s2)
  cur.ll=likePsi(theta[i,],t,beta[i],psi[i],tau,s2)
  prob=min(1,exp(cand.ll-cur.ll))    #prob is maxed at 1, or equal to difference between cand and cur.
  if(rbinom(1,1,prob)){psi[i]=cand;psiCount=psiCount+1}
  }
  

}
counter/M        #counters tell how often a candidate was accepted. Usually aim for .3 to .6, too far out of that range is bad
betaCount/M
psiCount/M
thetameans=apply(theta,2,mean)
plot(Data$acc,thetameans)
abline(0,1)

#More plots
plot(psi,typ='l')
plot(beta,typ='l')

#Data plot with line from posterior estimates
thetaplot=seq(.4,.9,by=.01)
forline=myExp(thetaplot,mean(beta),mean(psi),tau)
plot(Data$acc,Data$eff)
lines(thetaplot,forline,type='l')
