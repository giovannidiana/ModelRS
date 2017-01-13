lev<-as.numeric(cut(snr,breaks=50))
minred<-by(tab1$RED,lev,min)
maxred<-by(tab1$RED,lev,max)
s<-seq(0,2,length.out=length(minred))

fitmin <- lm( minred ~ poly(s,7));
fitmax <- lm( maxred ~ poly(s,7));

par(mar=c(6,6,2,2));
plot(0,xlim=c(0,2),ylim=c(-1,2),type='n',xlab="SNR",ylab="Redundancy",cex.axis=2,cex.lab=2);
polygon(c(s,rev(s)),c(predict(fitmin,data.frame(s),interval='confidence',level=0.99)[,1],
                      rev(predict(fitmax,data.frame(s),interval='confidence',level=0.99)[,1])),
		col='#006600');
lines(s,rep(0,length(s)))



