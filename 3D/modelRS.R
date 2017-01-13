library(R2Cuba)
library(mvtnorm)

# Parameters:
args<-commandArgs(trailingOnly = TRUE);
TOL=0.01;
ABSTOL=0.01;
minp=1e-40;
width=3.5;
EV=200000;
seed=as.numeric(args[1]);
set.seed(seed);
verb=F
intverb=0;
pin=.5;

s1=1;
s2=runif(1,0.1,s1);
s3=runif(1,0.1,s2);
alpha=runif(1,0,1);
beta=runif(1,0,1-alpha);
gamma=1-alpha-beta;

a=runif(1,0,2*s3);
b=a;
c=a;

det=-1; while(det<0){
	r12 = runif(1,-.7,.7)
	r13 = runif(1,-.7,.7)
	r23 = runif(1,-.7,.7)
	det = 1+2*r12*r13*r23-r12^2-r13^2-r23^2
}

s_sum=sqrt(alpha^2*s1^2+beta^2*s2^2+gamma^2*s3^2+2*alpha*beta*s1*s2*r12+2*alpha*gamma*s1*s3*r13+2*beta*gamma*s2*s3*r23);

if(verb==TRUE){
	cat(paste("** RUN PARAMS **\n",
          "a = ",a,"\n",
          "s1 = ",s1,"\n",
		  "s2 = ",s2,"\n",
		  "s3 = ",s3,"\n",
		  "r12 = ",r12,"\n",
		  "r13 = ",r13,"\n",
		  "r23 = ",r23,"\n",
		  "alpha = ",alpha,"\n",
		  "beta = ",beta,"\n",
		  "gamma = ",gamma,"\n",
		  "det = ",det,"\n",
		  sep="")
	)
}



C=matrix(c(s1^2,      r12*s1*s2,    r13*s1*s3,
           r12*s1*s2, s2^2,         r23*s2*s3,
		   r13*s1*s3, r23*s2*s3,    s3^2),
		 ncol=3);

pCh<-function(x){
	dmvnorm(x,mean=c(a,b,c),sigma=C);
}

pCl<-function(x){
	dmvnorm(x,mean=c(-a,-b,-c),sigma=C);
}

p<-function(x){
	pin*pCh(x)+(1-pin)*pCl(x);
}

p1Ch<-function(x){
	dnorm(x,mean=a,sd=s1);
}

p1Cl<-function(x){
	dnorm(x,mean=-a,sd=s1);
}

p1<-function(x){
	pin*p1Ch(x)+(1-pin)*p1Cl(x);
}

p2Ch<-function(x){
	dnorm(x,mean=b,sd=s2);
}

p2Cl<-function(x){
	dnorm(x,mean=-b,sd=s2);
}

p2<-function(x){
	pin*p2Ch(x)+(1-pin)*p2Cl(x);
}

p3Ch<-function(x){
	dnorm(x,mean=c,sd=s3);
}

p3Cl<-function(x){
	dnorm(x,mean=-c,sd=s3);
}

p3<-function(x){
	pin*p3Ch(x)+(1-pin)*p3Cl(x);
}

psumCh<-function(x){
	dnorm(x,mean=a,sd=sqrt(alpha^2*s1^2+beta^2*s2^2+gamma^2*s3^2+2*alpha*beta*s1*s2*r12+2*alpha*gamma*s1*s3*r13+2*beta*gamma*s2*s3*r23));
}

psumCl<-function(x){
	dnorm(x,mean=-a,sd=sqrt(alpha^2*s1^2+beta^2*s2^2+gamma^2*s3^2+2*alpha*beta*s1*s2*r12+2*alpha*gamma*s1*s3*r13+2*beta*gamma*s2*s3*r23))
}

psum<-function(x){
	pin*psumCh(x)+(1-pin)*psumCl(x);
}

Iden<-function(x){
	out=0;
	if(pCh(x)>minp){
		out=pin*pCh(x)*log(pCh(x)/p(x),2);
	}
	if(pCl(x)>minp){
		out=out+(1-pin)*pCl(x)*log(pCl(x)/p(x),2);
	}
	return(out);
}

I1den<-function(x){
	out=0;
	if(p1Ch(x)>minp){
		out=pin*p1Ch(x)*log(p1Ch(x)/p1(x),2);
	}
	if(p1Cl(x)>minp){
		out=out+(1-pin)*p1Cl(x)*log(p1Cl(x)/p1(x),2);
	}
	return(out);
}

I2den<-function(x){
	out=0;
	if(p2Ch(x)>minp){
		out=pin*p2Ch(x)*log(p2Ch(x)/p2(x),2);
	}
	if(p2Cl(x)>minp){
		out=out+(1-pin)*p2Cl(x)*log(p2Cl(x)/p2(x),2);
	}
	return(out);
}

I3den<-function(x){
	out=0;
	if(p3Ch(x)>minp){
		out=pin*p3Ch(x)*log(p3Ch(x)/p3(x),2);
	}
	if(p3Cl(x)>minp){
		out=out+(1-pin)*p3Cl(x)*log(p3Cl(x)/p3(x),2);
	}
	return(out);
}

ISden<-function(x){
	out=0;
	if(p1Ch(x[1])*p2Ch(x[2])*p3Ch(x[3])>minp){
		out=pin*p1Ch(x[1])*p2Ch(x[2])*p3Ch(x[3])*log(p1Ch(x[1])*p2Ch(x[2])*p3Ch(x[3])/(pin*p1Ch(x[1])*p2Ch(x[2])*p3Ch(x[3])+(1-pin)*p1Cl(x[1])*p2Cl(x[2])*p3Cl(x[3])),2);
	}
	if(p1Cl(x[1])*p2Cl(x[2])*p3Cl(x[3])>minp){
		out=out+(1-pin)*p1Cl(x[1])*p2Cl(x[2])*p3Cl(x[3])*log(p1Cl(x[1])*p2Cl(x[2])*p3Cl(x[3])/(pin*p1Ch(x[1])*p2Ch(x[2])*p3Ch(x[3])+(1-pin)*p1Cl(x[1])*p2Cl(x[2])*p3Cl(x[3])),2);
	}
	return(out);
}

Isumden<-function(x){
	out=0;
	if(psumCh(x)>minp){
		out=pin*psumCh(x)*log(psumCh(x)/psum(x),2);
	}
	if(psumCl(x)>minp){
		out=out+(1-pin)*psumCl(x)*log(psumCl(x)/psum(x),2);
	}
	return(out);
}

calc<-cuhre(3,1,Iden,lower=c(-a-width*s1,-b-width*s2,-c-width*s3),upper=c(a+width*s1,b+width*s2,c+width*s3),rel.tol=TOL,abs.tol=ABSTOL,flags=list(verbose=intverb),max.eval=EV,key=11);

I<-calc$value;
IS<-cuhre(3,1,ISden,lower=c(-a-width*s1,-b-width*s2,-c-width*s3),upper=c(a+width*s1,b+width*s2,c+width*s3),rel.tol=TOL,abs.tol=ABSTOL,flags=list(verbose=0),max.eval=EV,key=11)$value;


I1<-cuhre(1,1,I1den,lower=-a-width*s1,upper=a+width*s1,flags=list(verbose=0))$value;
I2<-cuhre(1,1,I2den,lower=-b-width*s2,upper=b+width*s2,flags=list(verbose=0))$value;
I3<-cuhre(1,1,I3den,lower=-c-width*s3,upper=c+width*s3,flags=list(verbose=0))$value;
Isum<-cuhre(1,1,Isumden,lower=-a-width*s_sum,upper=a+width*s_sum,flags=list(verbose=0))$value;

RED<-I1+I2+I3-I;
NC<-I-IS;
SC<-I1+I2+I3-IS;
cat(paste(seed,s1,s2,s3,a,r12,r13,r23,alpha,beta,gamma,I,I1,I2,I3,RED,SC,NC,Isum,calc$ifail,"\n",sep=" "));

