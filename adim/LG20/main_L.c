#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "momang.h"

int main (int argc, char *argv[]) {
    char filename[15];
    int nlines;
    FILE *input;
    double *L,Lmean;

    if (argc < 3
            || sscanf(argv[1],"%s",filename)!=1
            || sscanf(argv[2],"%d",&nlines)!=1) {
        fprintf(stderr,"%s: filename nlines\n",argv[0]);
        return 1;
    }


    L = (double *)malloc((nlines-1)*sizeof(double)); assert(L!=NULL);

    if ( (input = fopen(filename,"r")) == NULL) {
        fprintf(stderr,"%s: can't open file %s\n",argv[0],filename);
    }

    Lmean = calc_L(input,nlines,L);

    /*for (i=0; i<nlines-1; i++) {
        fprintf(stdout,"%10.6G\n",L[i]);
    }*/
    fprintf(stdout,"%10.6G\n",Lmean);

    return 0;
}
