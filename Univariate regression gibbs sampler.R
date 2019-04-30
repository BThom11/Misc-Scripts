## First Gibbs sampler (from Bayes lecture) for 205C.
#Univariate Linear Regression, a good intro case of bayesian parameter estimation

library(MCMCpack) #For invgamma sampling
#Data
x=c(1,2,3,4,5)
y=c(1,2,3,4,10)

#parameters needed for full conditionals from data
n=length(x)
mux=mean(x)
muy=mean(y)

#priors, B1 and B0 normal(M,S) sig2 inverse gamma (shape,scale)
B1M= .5
B1S=.5
B0M=0
B0S=2
shape=2
scale=2


#Initialize values: 
M= 1000 #number of loops
B0= 1:M
B0[1]=0   #B0 B1 and sig2 are parameters for model, arbitrary starting values
B1=1:M
B1[1]=2
sig2=1:M
sig2[1]=8

#Gibbs sampler

for(i in 2:M){
  #for B0. Sampled from normal (A2*B,A2), those parameters defined below and use priors from above
  B0A2=((n/sig2[i-1])+1/(B0S))^-1
  B0B=(n/sig2[i-1])*(muy- B1[i-1]*mux) + (B0M/B0S)
  B0[i]=rnorm(1,B0A2*B0B,sqrt(B0A2))
  
  #for B1
  B1A2=((n/sig2[i-1])+1/(B1S))^-1
  B1B= (n/sig2[i-1])*(muy-B0[i]*mux)+ (B1M/B1S)
  B1[i]=rnorm(1,B1A2*B1B,sqrt(B1A2))
  
  #for sig2
  shapeSig2= shape+(n/2)
  scaleSig2= scale+ sum((y-(B0[i]+B1[i]*x))^2)/2
  sig2[i]=rinvgamma(1,shapeSig2,scale=scaleSig2)
}

keep=50  #burn in period and filtering data respectivley
B0=B0[keep:M]
B1=B1[keep:M]
sig2=sig2[keep:M]

mean(B0)   #mean estimates
mean(B1)
mean(sig2)

par(mfrow=c(2,2))  #make plot output 2x2 grid, one plot in each
plot(x,y)   #raw data and regression line
abline(a=mean(B0),b=mean(B1),col='red')
hist(B0)      #hist for B0 with red line as mean
abline(v=mean(B0),col='red')
hist(B1)      #hist for B1 with red line as mean
abline(v=mean(B1),col='red')
hist(sig2[sig2<=200],breaks=seq(0,200,by=10),xlim=c(0,200))  #hist for sig 2 (removing high values), red as mean, blue median
abline(v=mean(sig2),col='red')
abline(v=median(sig2),col='blue')


