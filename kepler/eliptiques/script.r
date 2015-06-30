source('elipsefit.r')

tr<-read.table(paste0(nom,".dat"),quote="\"",comment.char = "#");
x=tr$V1;
y=tr$V2;
efit<-fit.ellipse(x,y);
efit
c<-sqrt(efit$major^2-efit$minor^2)
theta<-sqrt(efit$coef[2]/(efit$coef[1]-efit$coef[3]+sign(efit$coef[2])*sqrt(efit$coef[2]^2+(efit$coef[1]-efit$coef[3])^2)))
focusx<-efit$center[1]-c*cos(theta)
focusy<-efit$center[2]-c*sin(theta)
write.table(c(focusx,focusy),file=paste0(nom,"-focus.dat"), row.names = FALSE, col.names = FALSE);
e<-get.ellipse(efit,n=801);
write.table(e,file=paste0(nom,"-fit.dat"));
plot(x,y)
lines(e,col="red")