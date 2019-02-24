#Basic Sampling process to infer group means assuming same variance. 
#Hand coded the chain, but could also be done in JAGS/STAN

library(MCMCpack)

#data-----
set.seed(234)
X=rnorm(10,98,15)
Y=rnorm(15,102,15)
meanX=mean(X)
meanY=mean(Y)
NY=length(Y)
NX=length(X)

#priors-----
a=100     # sample mean x & y ~ Normal(a,b)
b=15^2
r=1       # sample variance ~ InvGamma(shape=r,scale=s)
s=10^2

M= 1500
muX=1:M     #Creating vectors to fill
muY=1:M
sig=1:M

#start values-----
muX[1]=100   #initial guesses, can make them a lot worse without any real issues
muY[1]=100
sig[1]=10

#actual sampling and updating
for(i in 2:M){
  vY=1/((NY/sig[i-1])+(1/b))
  cY=((NY*meanY/sig[i-1])+(a/b))
  vX=1/((NX/sig[i-1])+(1/b))
  cX=((NX*meanX/sig[i-1])+(a/b))
  muY[i]=rnorm(1,vY*cY,sqrt(vY))
  muX[i]=rnorm(1,vX*cX,sqrt(vX))
  sig.scale=((sum((Y-muY[i])^2)+sum((X-muX[i])^2))/2)+s
  sig[i]=rinvgamma(1,shape=r+(NY/2+NX/2),scale=sig.scale)
}

##plots of posteriors----
keep=50:M    #remove 50 trial burn in

#Basic plotting function
Hist.Dens=function(data,label){
  hist1=hist(data,plot=F,breaks=20)
  dens=density(data)
  hist(data,probability=T,breaks=20,ylim=c(0,max(hist1$density)+.01),main=label)
  lines(dens,col='red')
  abline(v=mean(data),col='blue')
}

#create hists of sample with density as posterior. Order is muX, muY, Sigma, and Diff of muX muY
Hist.Dens(muX[keep],'muX')
Hist.Dens(muY[keep],'muY') 
Hist.Dens(sqrt(sig[keep]),"Sigma")
Hist.Dens(muX[keep]-muY[keep],'Difference')
