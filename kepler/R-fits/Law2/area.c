#include <stdio.h>
#include <math.h>

int main (int argc, char *argv[]) {
    char nomFitxer[30];
    double focus;
    double x[2],y[2],t[2];
    double theta[2],dtheta,r;
    double foo;
    FILE *fitxer;

    //Lectura paràmetres
    if (argc < 2
            || sscanf(argv[1],"%s",nomFitxer)!=1) {
        fprintf(stderr,"%s: fitxer\n", argv[0]);
        return 1;
    }

    fprintf(stderr,"Llegint focus\n");
    fscanf(stdin,"%lf",&focus);
    fprintf(stderr,"Focus llegit\n");

    if ( (fitxer = fopen(nomFitxer,"r")) == NULL ) {
        fprintf(stderr,"Error reading %s\n",nomFitxer);
        return 1;
    }

    fprintf(stdout,"#x y t dA/dt\n");
    
    fscanf(fitxer,"%lf %lf %lf %lf",&x[0],&y[0],&t[0],&foo);
    theta[0] = atan(y[0]/x[0]);
    while (!feof(fitxer)) {
        fscanf(fitxer,"%lf %lf %lf %lf",&x[1],&y[1],&t[1],&foo);
        r = sqrt( (x[0]-focus)*(x[0]-focus)+y[0]*y[0] );
        theta[1] = atan(y[1]/x[1]);
        dtheta = (theta[1]-theta[0])/(t[1]-t[0]);
        fprintf(stdout,"%10.6G %10.6G %10.6G %10.6G\n",x[1],y[1],t[1],0.5*r*r*dtheta);
        x[0]=x[1];
        y[0]=y[1];
        theta[0]=theta[1];
    }

    fclose(fitxer);
    return 0;
}
