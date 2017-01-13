require(RColorBrewer);

colpal_196<-brewer.pal(n=7,"Greys");
colpal_402<-brewer.pal(n=7,"Reds");
colpal_404<-brewer.pal(n=7,"Blues");
colpal_435<-brewer.pal(n=7,"Purples");
colpal=c(colpal_196[7], colpal_402[7], colpal_404[7], colpal_435[7]);

snr_adf <- function(gt,t2){
	
	mydata <- read.table("/home/diana/workspace/data/dataset-add-2/expression5.dat", header=FALSE);
	GT=gt
	T2=t2
	data0 <- subset(mydata,V3==GT & V6==T2 & V4==6);
	names(data0)<-c("B","BD","GT","Day","F2","T2","ADF","ASI","NSM");
	w<-table(data0$F2)/nrow(data0)
	means=as.numeric(by(data0$ADF,data0$F2,mean))
	globalmean=mean(data0$ADF);
	vars=as.numeric(by(data0$ADF,data0$F2,var))
	a=sum(w*(means-globalmean)^2)
	s3=sum(w*vars)
	
	return(sqrt(a/s3))
}

b<-c(snr_adf("QL196",20),snr_adf("QL404",20),snr_adf("QL402",20),snr_adf("QL435",20))
pdf("snr.pdf")
par(mar=c(15,6,2,2))
pos<-barplot(b,ylim=c(0,.7),ylab="ADF SNR", cex.axis=2,cex.lab=2,
             col=colpal)
#plot(0,type='n',frame=F,axes=F,ylim=c(0,.1),xlim=range(pos),xlab="",ylab="");
text(x=pos,y=rep(-.06,4),pos=2, labels = c("WT","tph-1(-)","daf-7(-)","tph-1(-);daf-7(-)"),srt = 60, xpd = T, cex=2)

dev.off()

