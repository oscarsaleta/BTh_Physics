/* gcc -Wall -o area area.c -lm */

#include <stdio.h>
#include <math.h>

int main (int argc, char *argv[]) {
    char nomFitxer[30];
    double focus[2];
    double x[2],y[2],t[2];
    double theta[2],dtheta,r[2];
    double theta_f[2],dtheta_f,r_f;
    double energy,darea;
    FILE *fitxer;

    //Lectura paràmetres
    // Hem de llegir el fitxer on estan guardats els punts de dades
    if (argc < 2
            || sscanf(argv[1],"%s",nomFitxer)!=1) {
        fprintf(stderr,"%s: fitxer < focus > output\n", argv[0]);
        return 1;
    }

    //Lectura del focus (stdin)
    fprintf(stderr,"Llegint focus\n");
    fscanf(stdin,"%lf",&focus[0]);
    fscanf(stdin,"%lf",&focus[1]);
    fprintf(stderr,"Focus llegit: (%g,%g)\n",focus[0],focus[1]);

    //Obrir fitxer
    if ( (fitxer = fopen(nomFitxer,"r")) == NULL ) {
        fprintf(stderr,"Error reading %s\n",nomFitxer);
        return 1;
    }

    fprintf(stdout,"#x y t dA/dt K\n");

    //Càlculs
    fprintf(stderr,"Calculant angles inicials respecte el focus...\n");
    fscanf(fitxer,"%lf %lf %lf %lf",&x[0],&y[0],&t[0],&r[0]); //llegir punt inicial
    theta[0] = atan2(y[0],x[0]); //angle resp (0,0)
    theta_f[0] = atan2(y[0]-focus[1],x[0]-focus[0]); //angle resp focus

    fprintf(stderr,"Processant fitxer...");
    while (!feof(fitxer)) {
        fscanf(fitxer,"%lf %lf %lf %lf",&x[1],&y[1],&t[1],&r[1]);
        r_f = sqrt( (x[0]-focus[0])*(x[0]-focus[0])+(y[0]-focus[1])*(y[0]-focus[1]) );
        theta[1] = atan2(y[1],x[1]);
        dtheta = (theta[1]-theta[0])/(t[1]-t[0]); //angle resp (0,0)
        theta_f[1] = atan2(y[1]-focus[1],x[1]-focus[0]); //angle resp focus
        dtheta_f = (theta_f[1]-theta_f[0])/(t[1]-t[0]);
        darea = 0.5*r_f*r_f*dtheta_f;
        energy = 0.5*r[0]*r[0]*dtheta*dtheta-2./r[0];
        fprintf(stdout,"%10.6G %10.6G %10.6G %10.6G %10.6G %10.6G\n",x[1],y[1],t[1],fabs(darea),energy,r[0]*dtheta);
        x[0]=x[1];
        y[0]=y[1];
        t[0]=t[1];
        r[0]=r[1];
        theta[0]=theta[1];
        theta_f[0]=theta_f[1];
    }
    fprintf(stderr," Done\n");

    fclose(fitxer);
    return 0;
}
