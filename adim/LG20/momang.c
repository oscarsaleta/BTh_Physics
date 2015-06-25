#include <stdio.h>
#include <stdlib.h>

#define M 1

double calc_L (FILE *input, int nlines, double *L) {
    int bsize = 100;
    int i,j;
    char buffer[100],*str;
    double aux[nlines][3];
    double Lmean=0, dx,dy,dt,vx,vy;

    /* Lectura */
    for (i=0; i<nlines; i++) {
        fgets(buffer,bsize,input);
        str=buffer;
        /* aux[0] = x
         * aux[1] = y
         * aux[2] = t
         */
        for (j=0; j<3; j++)
            aux[i][j] = strtod(str,&str);

        if (i>0) {
            dx = aux[i][0]-aux[i-1][0];
            dy = aux[i][1]-aux[i-1][1];
            dt = aux[i][2]-aux[i-1][2];
            vx = dx/dt;
            vy = dy/dt;
            L[i] = aux[i][0]*M*vy-aux[i][1]*M*vx;
            Lmean += L[i];
        }

    }

    return Lmean/((double)nlines);

}
