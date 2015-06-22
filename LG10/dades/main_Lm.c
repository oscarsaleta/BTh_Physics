#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "mean.h"

int main (int argc, char *argv[]) {
    char filename[15];
    int nlines;
    FILE *input;
    double mean;

    if (argc < 3
            || sscanf(argv[1],"%s",filename)!=1
            || sscanf(argv[2],"%d",&nlines)!=1) {
        fprintf(stderr,"%s: filename nlines\n",argv[0]);
        return 1;
    }

    if ( (input = fopen(filename,"r")) == NULL) {
        fprintf(stderr,"%s: can't open file %s\n",argv[0],filename);
    }

    mean = calc_mean(input,nlines);

    fprintf(stdout,"%10.6G\n",mean);

    return 0;
}
