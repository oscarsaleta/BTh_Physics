source('/home/slenderman/git/tdg-fisica/software/R/elipsefit.r')

nom=c('r1mQ')
filedir="/home/slenderman/git/tdg-fisica/octubre/Radis_Bohr/"
nom=paste0(filedir,nom)

for(i in 1:length(nom)) {
  tr<-read.table(paste0(nom[i],".dat"),quote="\"",comment.char = "#");
  x=tr$V1;
  y=tr$V2;
  efit<-fit.ellipse(x,y);
  efit
  c<-sqrt(efit$major^2-efit$minor^2)
  theta<-sqrt(efit$coef[2]/(efit$coef[1]-efit$coef[3]+sign(efit$coef[2])*sqrt(efit$coef[2]^2+(efit$coef[1]-efit$coef[3])^2)))
  focusx<-efit$center[1]-c*cos(theta)
  focusy<-efit$center[2]-c*sin(theta)
  write.table(c(focusx,focusy),file=paste0(nom[i],"-focus.dat"), row.names = FALSE, col.names = FALSE);
  e<-get.ellipse(efit,n=801);
  write.table(e,file=paste0(nom[i],"-fit.dat"));
  plot(x,y)
  lines(e,col="red")
}