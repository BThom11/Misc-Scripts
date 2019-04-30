#Nelder Mead optimizer on a Simonescu function. Optimizer is not generalizable to N dimensions, but easily coult
#be converted to be. Gif at the end created using gganimate.

library(ggplot2)
library(gganimate)
library(gifski)
library(png)
library(dplyr)

#Creating data to fill Simonescu function plot with parameters rT,rS and n. (see wiki on nelder-mead and follow from there)
#x values are found in coord[,1], y values in coord[,2], z or color found in z as x*y*.1 (with constraints)
x= seq(-1.25,1.25,by=.01)
y=seq(-1.25,1.25,by=.01)
rT=1
rS=.2
n=8
count=1
coord=matrix(0,nrow=(251*251),ncol=2)
#Loop here checks if x and y satisfy constraint of simonescu function with above parameters.
#If they do, coord pair is saved, otherwise not saved.
for(i in 1:length(x)){   
  for(j in 1:length(y)){
    if(y[j]==0){
      next
    }
    if(x[i]^2+y[j]^2<= (rT+rS*cos(n*atan(x[i]/y[j])))^2){
      coord[count,]=c(x[i],y[j])
      count=count+1
    }}}
coord=coord[1:count,]
z=coord[,1]*coord[,2]*.1   #Actual value of third variable given by Simonescu function.



#optimization using Nelder-Mead.

toZ=function(x){             #Function calculates f(x,y) of a matrix with constraints from rT,rS, and n
  val=c(0,0,0)
  for (i in 1:nrow(x)){
    if(x[i,2]==0){
      next
    }
    if((x[i,1]^2+x[i,2]^2)<=(rT+rS*cos(n*atan(x[i,1]/x[i,2])))^2){
      val[i]=x[i,1]*x[i,2]*.1
    }
  }
  return(val)
}
toZvec=function(x){   #Does the same as function above just with single vector, could consolidate into 1 function
  val=0
  if(x[2]==0){
    next
  }
  if((x[1]^2+x[2]^2) <= (rT+rS*cos(n*atan(x[1]/x[2])))^2){
    val=x[1]*x[2]*.1
  }
  return(val)
}
p1=c(0,0)   #p1, 2, and 3 are starting values of points, cannot make them too close together or else local minima
p2=c(.5,.2)  #can change them to have optimizer go towards one of two global minimum.
p3=c(.7,-.3)  
alpha=1   #alpha, gamma, rho, and sigma are all adjustment parameters for reflection, expansion, contraction respectivley 
gamma=2
rho=.5
sigma=.5
loops=25   #Number of optimizer steps

points=matrix(c(p1,p2,p3),nrow=3,byrow=T)  #creates matrix of single state
points_full=array(rep(as.vector(points),loops),dim=c(3,2,loops))   #creates array, 3rd index captures state/frame

for(i in 1:loops){  #Run optimizer, see Nelder-mead Wiki for details. No termination function
  val=toZ(points)
  order=sort(val,index.return=T)[[2]]   #Order
  x0=c((points[order[1],1]+points[order[2],1])/2,(points[order[1],2]+points[order[2],2])/2)
  xR = x0+alpha*(x0-points[order[3],])
  fxR= toZvec(xR)
  
  if(val[order[1]]<=fxR & fxR < val[order[2]]){   #Reflection
    points[order[3],]=xR
    points_full[,,i]=points
    next
  }
  if(fxR< val[order[1]]){    #Expansion
    xE=x0+gamma*(xR-x0)
    fxE= toZvec(xE)
    if(fxE<fxR){
      points[order[3],]=xE
      points_full[,,i]=points
      next
    }
    else{points[order[3],]=xR
    points_full[,,i]=points
    next }}
  if(fxR >= val[order[2]]){      #Contraction and else case handles shrink
    xC=x0+ rho*(points[order[3],]-x0)
    fxC=toZvec(xC)
    if(fxC < val[order[3]]){
      points[order[3],]=xC
      points_full[,,i]=points
      next
    }
    else{ points[order[2],]=points[order[1],]+ sigma*(points[order[2],]-points[order[1],])
    points[order[3],]=points[order[1],]+sigma*(points[order[3],]-points[order[1],])
    points_full[,,i]=points
    next}
  }
}

pointsDF=data.frame(xcoord=as.vector(points_full[,1,]),ycoord=as.vector(points_full[,2,]),frame=rep(1:loops,each=3))
 
#Line above consolodates data from array to DF

ggplot(data=NULL,aes(x=coord[,1],y=coord[,2],colour=z))+   #ggplot and animate call.
  geom_point()+
  geom_point(data=pointsDF,aes(x=xcoord,y=ycoord),col='red',size=2.5)+
  labs(title="State: {previous_state}",x=NULL,y=NULL)+   #stuff below here is gganimate:
  transition_states(frame,transition_length=0,state_length=10)+ # reflects transition state, i.e. each loop through optimizer
  shadow_mark(size=.6,alpha=.6)  #leaves behind smaller faded points of past states.

points #final call showing the current values of the three points (should be basically optimal if ran long enough)


  

