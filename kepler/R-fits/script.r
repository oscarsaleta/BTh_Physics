source('elipsefit.r')

nom="tr10";
tr<-read.table(paste0(nom,".dat"),quote="\"",comment.char = "#");
x=tr$V1;
y=tr$V2;
efit<-fit.ellipse(x,y);
efit
c<-sqrt(efit$major^2-efit$minor^2)
focus<-efit$center[1]-c
write.table(focus,file=paste0(nom,"-focus.dat"), row.names = FALSE, col.names = FALSE);
e<-get.ellipse(efit,n=801);
write.table(e,file=paste0(nom,"-fit.dat"));
plot(x,y)
lines(e,col="red")

nom="tr050";
tr<-read.table(paste0(nom,".dat"),quote="\"",comment.char = "#");
x=tr$V1;
y=tr$V2;
efit<-fit.ellipse(x,y);
efit
c<-sqrt(efit$major^2-efit$minor^2)
focus<-efit$center[1]-c
write.table(focus,file=paste0(nom,"-focus.dat"), row.names = FALSE, col.names = FALSE);
e<-get.ellipse(efit,n=801);
write.table(e,file=paste0(nom,"-fit.dat"));
plot(x,y)
lines(e,col="red")

nom="tr20";
tr<-read.table(paste0(nom,".dat"),quote="\"",comment.char = "#");
x=tr$V1;
y=tr$V2;
efit<-fit.ellipse(x,y);
efit
c<-sqrt(efit$major^2-efit$minor^2)
focus<-efit$center[1]-c
write.table(focus,file=paste0(nom,"-focus.dat"), row.names = FALSE, col.names = FALSE);
e<-get.ellipse(efit,n=801);
write.table(e,file=paste0(nom,"-fit.dat"));
plot(x,y)
lines(e,col="red")


nom="ab10";
tr<-read.table(paste0(nom,".dat"),quote="\"",comment.char = "#");
x=tr$V1;
y=tr$V2;
efit<-fit.ellipse(x,y);
efit
c<-sqrt(efit$major^2-efit$minor^2)
focus<-efit$center[1]-c
write.table(focus,file=paste0(nom,"-focus.dat"), row.names = FALSE, col.names = FALSE);
e<-get.ellipse(efit,n=801);
write.table(e,file=paste0(nom,"-fit.dat"));
plot(x,y)
lines(e,col="red")