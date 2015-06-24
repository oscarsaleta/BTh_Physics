#include <stdio.h>
#include <stdlib.h>

#define M 1

double calc_mean (FILE *input, int nlines) {
    int bsize = 100;
    int i,j;
    char buffer[100],*str;
    double aux[2];
    double mean=0;

    /* Lectura */
    for (i=0; i<nlines; i++) {
        fgets(buffer,bsize,input);
        str=buffer;
        for (j=0; j<3; j++)
            aux[j] = strtod(str,&str);
        mean+=aux[1];
    }

    return mean/((double)nlines);

}
