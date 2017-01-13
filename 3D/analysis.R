
tab0<-read.table("output_merge.dat");
#tab0<-read.table("output_SCneg.dat");

names(tab0)<-c("seed","s1","s2","s3","a","r12","r13","r23","alpha","beta","gamma","I","I1","I2","I3","RED","SC","NC","Isum","ifail");
tabnan<-tab0[is.nan(apply(tab0[,10:19],1,sum)),]
rbPal <- colorRampPalette(c('blue','red'))
gPal <- colorRampPalette(c('#00FF00','#006600'))
tab0$Col <- rbPal(20)[as.numeric(cut(tab0$RED,breaks = 20))]
tab0$Col <- rep(NA,nrow(tab0));
tab0$Col[tab0$RED>0] <-rgb(1,0,0,1)
tab0$Col[tab0$RED<0] <-rgb(0,0,1,1)
#tab0<-subset(tab0,ifail==0);
tab1<-tab0[tab0$ifail==0 & tab0$NC>-.5 & tab0$I<1 
          & tab0$Isum<tab0$I & tab0$Isum>0 & tab0$I>tab0$I3,];

show<-function(cond){
    #png("table.png",width=12,height=18,units='in',res=300);
    #pdf("table.pdf",width=12,height=18);
    par(mar=c(6,6,2,2));
	tab<-tab1[cond,];
	layout(matrix(c(1:6),ncol=2));
	plot(tab$I,tab$RED,col=tab$Col,cex=.3,pch=16,cex.axis=2,cex.lab=2,xlab="Information",ylab="Redundancy")
	plot(tab$RED,tab$Isum,col=tab$Col,cex=.3,pch=16,cex.axis=2,cex.lab=2,xlab="Redundancy",ylab="Transmitted Information")
	plot(tab$I,tab$Isum,col=tab$Col,cex=.3,pch=16,cex.axis=2,cex.lab=2,xlab="Information",ylab="Transmitted Information")
	plot(tab$SC,tab$NC,col=tab$Col,cex=.3,pch=16,cex.axis=2,cex.lab=2,xlab="Signal Correlation",ylab="Noise Correlation")
	plot(tab$I2,tab$I3,col=tab$Col,cex=.3,pch=16,cex.axis=2,cex.lab=2,xlab=expression(I[2]),ylab=expression(I[3]))
	plot(tab$I,tab$NC,col=tab$Col,cex=.3,pch=16,cex.axis=2,cex.lab=2,xlab="Information",ylab="Noise Correlation")
	#dev.off()
}

as3X2<-function(){
    #png("as3X2.png",width=12,height=6,units='in',res=300);
    pdf("as3X2.pdf",width=12,height=6);
    par(mar=c(6,6,2,2));
	tab_left<-tab1;
	tab_right<-tab1[tab1$I>0.3,];
	layout(matrix(c(1:2),ncol=2));
	plot(tab_left$a,tab_left$s3,col=tab_left$Col,cex=.3,pch=16,cex.axis=2,cex.lab=2,xlab="a",ylab=expression(sigma[3]))
	plot(tab_right$a,tab_right$s3,col=tab_right$Col,cex=.3,pch=16,cex.axis=2,cex.lab=2,xlab="a",ylab=expression(sigma[3]))
	dev.off()
}

aNC<-function(cond){
	tab<-tab1[cond,];
	tab<-tab[tab$ifail==0 & tab$NC>-.5 & tab$I<1 
	          & tab$Isum<tab$I & tab$Isum>0 & tab$I>tab$I3,];
	plot(tab$a,tab$NC,col=tab$Col,cex=.3)

}

REDSNR<-function(cond){
	tab<-tab1[cond,];
    REDpalette <- rbPal(50)[as.numeric(cut(tab$NC,breaks = 50))]
	plot(tab$a/tab$s3,tab$RED,pch=19,col=REDpalette,cex=.3,xlab="SNR",ylab="Redundancy",cex.axis=2,cex.lab=2)
}

as3<-function(cond){
	tab<-tab1[cond,];
	plot(tab$a,tab$s3,col=tab$Col,cex=.3,xlab="Dynamic range",ylab="Noise",cex.axis=2,cex.lab=2)
    lines(seq(0.05,.5,length.out=30),2*(seq(0.05,.5,length.out=30)),lty=2,lwd=3,col="yellow");

}

RedTI<-function(cond){
	tab<-tab1[cond,];
	plot(tab$RED,tab$Isum,col=tab$Col,cex=.3,pch=16,cex.axis=2,cex.lab=2,xlab="Redundancy",ylab="Transmitted Information")
}

png("redsyn.png",bg="transparent",width=6,height=6,units='in',res=300);
#pdf("redsyn.pdf",bg="transparent",width=6,height=6);
par(mar=c(6,6,2,2));
plot(tab1$RED,tab1$Isum,col=tab1$Col,cex=.3,xlab="Redundancy",ylab="Transmitted Information",cex.axis=2,cex.lab=2,xlim=c(-1.8,1.8),pch=16);
dev.off();

png("as3.png",bg="transparent");
par(mar=c(6,6,2,2));
as3(tab1$I>0.4)
dev.off();

inset <- function(){
	hpos<-hist(tab1$Isum[tab1$RED>0],breaks=50);
	hneg<-hist(tab1$Isum[tab1$RED<0],breaks=50);
	hnorm<-hist(tab1$Isum,breaks=50);
	#pdf("inset.pdf");
	par(mar=c(6,6,2,2));
	plot(seq(0,1,length.out=50),hpos$counts/hnorm$counts,type='l',col="red",lwd=3, xlab="Transmitted Information",ylab="",cex.axis=2,cex.lab=2,ylim=c(0,1))
	lines(seq(0,1,length.out=50),hneg$counts/hnorm$counts,col="blue",lwd=3)
	#dev.off();
}

# Fraction of synergy at fixed Info and increasing a/s3
inset2 <- function(max,Imin,B=40,print=F){

	hpos<-hist((tab1$a/tab1$s3)[tab1$RED<0 & tab1$I>Imin],breaks=seq(0,max,length.out=B));
	hnorm<-hist((tab1$a/tab1$s3)[tab1$I>Imin],breaks=seq(0,max,length.out=B));
	if(print==T) pdf("inset2_I0.pdf");
	par(mar=c(6,6,2,2));
	hrat<-hpos$counts/hnorm$counts;
	if(Imin>0)hrat[1]=1;
	plot(seq(0,max,length.out=B-1),hrat,type='l',lwd=3, xlab=expression(a/sigma[3]),ylab="",cex.axis=2,cex.lab=2,ylim=c(0,1))
	cxblue<-c(0,seq(0,max,length.out=B-1),max);
	cyblue<-c(0,hrat,0);
	if(Imin>0)	cyblue[2]=1;

	cxred<-c(0,seq(0,max,length.out=B-1),max);
	cyred<-c(1,hrat,1.0);
	if(Imin>0) cyred[2]=1;
	polygon(cxblue,cyblue,col="blue",lwd=3)
	polygon(cxred,cyred,col="red",lwd=3)
	if(Imin>0) arrows(.5,0,.5,1, lty=2,length=0,lwd=3,col="yellow");
	if(print==T) dev.off();
}

inset3 <- function(max,Imin,B=40){

	hpos<-hist((tab1$a/tab1$s3)[tab1$RED<0 & tab1$I>Imin],breaks=seq(0,max,length.out=B),plot=F);
	hnorm<-hist((tab1$a/tab1$s3)[tab1$I>Imin],breaks=seq(0,max,length.out=B),plot=F);
	par(mar=c(6,6,2,2));
	hrat<-hpos$counts/hnorm$counts;
	if(Imin>0)hrat[1]=1;
	image(seq(0,2,length.out=length(hrat)),0:1,matrix(hrat,ncol=1),
	      col=rgb(hrat,rep(0,length(hrat)),1-hrat),
		  yaxt='n',ylab="",cex.axis=2,xlab="SNR")
    lines(seq(0,2,length.out=length(hrat)),hrat,col="white",cex=3)

}

pdf("inset3.pdf",width=8,height=2)
inset3(2,0.4)
dev.off();


