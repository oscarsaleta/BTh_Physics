#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#include "numerical_lib.h"
#include "quantum_lib.h"

#define QUANT_h 1
#define QUANT_m 1
#define QUANT_w 1

#define F(i,j) f[(fstruct.xdim*(j)+(i))]

double V(double x, double t);
double Vx(double x, double t);
double Vy(double x, double t);
double U(double x, double y, double t);

complex double stateZero2D(double x, double y, double x0, double y0, double px0, double py0);
complex double stateZero(double x, double x0, double px0);
complex double stateOne(double x, double x0);

int main(int argc, char* argv[]){

    /* Program parameters */
    State fstruct;
    double xi=-8.,xf=8.,yi=-8.,yf=8.;
    complex double *f;
    double r[2]; // posició inicial de la partícula
    double tmax;
    double dt;
    int i;
    /* File for reading */
    char filename[10];
    FILE *input;

    /* READING PARAMETERS */
    if (argc < 3
            || sscanf(argv[1],"%d",&nt)!=1
            || sscanf(argv[5],"%s",filename)!=1) {
        fprintf(stderr,"%s: nt filename\n",argv[0]);
        return 1;
    }

    /* INITIALIZE STRUCTURE AND STATE*/
    fstruct.xdim = fstruct.ydim = 256;
    /* Allocation of arrays and pointers*/
    fstruct.t = (double *)malloc(sizeof(double));
    fstruct.dt = (double *)malloc(sizeof(double));
    fstruct.x = (double *)malloc(fstruct.xdim*sizeof(double)); assert(fstruct.x != NULL);
    fstruct.y = (double *)malloc(fstruct.ydim*sizeof(double)); assert(fstruct.y != NULL);
    f = (complex double *)malloc(fstruct.xdim*fstruct.ydim*sizeof(complex double)); assert(f != NULL);
    /* Initiate state structure */
    *fstruct.t = 0;
    *fstruct.dt = dt;
    fstruct.t_max = tmax;
    vec_create(xi,xf,yi,yf,&fstruct);


    /* Actual commands */
    input = fopen(filename,"r");

    for(i=0;i<nt;i++) {
        fscanf(input,"%lf %lf %lf",&r[0],&r[1],&t);
        fprintf(stdout,"%10.6G %10.6G\n",Qj(r,

    fclose(input);
    free(fstruct.x); free(fstruct.y); free(f);

    return 0;
}

double V(double x, double t){
    return 0.5*x*x;
}

double Vx (double x, double t){
    if (x<0)
        return 0.5*QUANT_m*QUANT_w*QUANT_w*(x+2)*(x+2);
    return 0.5*QUANT_m*QUANT_w*QUANT_w*(x-2)*(x-2);
}

double Vy (double x, double t){
    return 0.5*QUANT_m*QUANT_w*QUANT_w*x*x;
}
double U (double x, double y, double t){
    return 0;
}

complex double stateZero2D(double x, double y, double x0, double y0, double px0, double py0) {
    complex double fx, fy;
    fx = exp(-QUANT_m*QUANT_w*((x-x0)*(x-x0))/(2.*QUANT_h))*cexp(I*px0*x/QUANT_h);
    fy = exp(-QUANT_m*QUANT_w*((y-y0)*(y-y0))/(2.*QUANT_h))*cexp(I*py0*y/QUANT_h);
    return sqrt(QUANT_m*QUANT_w/(M_PI*QUANT_h))*fx*fy;
}

complex double stateZero(double x, double x0, double px0) {
    complex double fx;
    double alpha = QUANT_m*QUANT_w/QUANT_h;
    fx = exp(-alpha*(x-x0)*(x-x0)/2);
    return sqrt(sqrt(alpha/M_PI))*fx;
    //fx = exp(-QUANT_m*QUANT_w*((x-x0)*(x-x0))/(2.*QUANT_h))*cexp(I*px0*x/QUANT_h);
    //return sqrt(sqrt(QUANT_m*QUANT_w/(M_PI*QUANT_h)))*fx;
}

complex double stateOne(double x, double x0) {
    double alpha = QUANT_w*QUANT_m/QUANT_h;
    double y = sqrt(alpha)*(x-x0);
    complex double fx = sqrt(2.)*y*exp(-(y*y)/2.);
    return sqrt(sqrt(alpha/M_PI))*fx;
//    double x0 = sqrt(QUANT_h/QUANT_m/QUANT_w);
//    return 1/sqrt(2.*sqrt(x0*M_PI))*2*x/x0*exp(-x*x/(2.*x0*x0));
}

#undef F

#undef QUANT_h
#undef QUANT_m
#undef QUANT_w
