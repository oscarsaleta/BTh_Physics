#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>

#include "numerical_lib.h"
#include "quantum_lib.h"

#define CENTRE_POU 1.5
#define ALPHA 1
#define IM_MAXIT 5000

#define F(i,j) f[(fstruct.xdim*(j)+(i))]
#define ABSPSISQ(i,j) abspsisq[(fstruct.xdim*(j)+(i))]

double V(double x, double t);
double U(double x, double y, double t);

complex double stateZero2D(double x, double y, double x0, double y0);
complex double stateZeroOne(double x, double y, double fase);
complex double stateZero(double x, double x0);
complex double stateOne(double x, double x0);

int main(int argc, char* argv[]){

    /* Program parameters */
    State fstruct;
    double xi=-8.,xf=8.,yi=-8.,yf=8.;
    complex double *f;
    double dt;
    double tmax;
    double x0=-CENTRE_POU,y0=-0;
    double r[2];
    double fac; //factor de la fase
    /* Dummy variables */
    int i,j;
    complex double *abspsisq; /*Dummy psi storage*/

    /* READING PARAMETERS */
    if (argc < 6
            || sscanf(argv[1],"%lf",&dt)!=1
            || sscanf(argv[2],"%lf",&tmax)!=1
            || sscanf(argv[3],"%lf",&fac)!=1
            || sscanf(argv[4],"%lf",&r[0])!=1
            || sscanf(argv[5],"%lf",&r[1])!=1
            ) {
        fprintf(stderr,"%s: dt tmax fase x0 y0\n",argv[0]);
        return 1;
    }

    /* INITIALIZE STRUCTURE AND STATE*/

    fstruct.xdim = fstruct.ydim = 151;

    /* Allocation of arrays and pointers*/
    fstruct.t = (double *)malloc(sizeof(double));
    fstruct.dt = (double *)malloc(sizeof(double));
    fstruct.x = (double *)malloc(fstruct.xdim*sizeof(double)); assert(fstruct.x != NULL);
    fstruct.y = (double *)malloc(fstruct.ydim*sizeof(double)); assert(fstruct.y != NULL);
    f = (complex double *)malloc(fstruct.xdim*fstruct.ydim*sizeof(complex double)); assert(f != NULL);
    abspsisq = (complex double *)malloc(fstruct.xdim*fstruct.ydim*sizeof(complex double)); assert(abspsisq != NULL);


    /* Initiate state structure */
    *fstruct.t = 0;
    *fstruct.dt = dt;
    fstruct.t_max = tmax;
    vec_create(xi,xf,yi,yf,&fstruct);

    /* Initial State */ 
    for (i=0; i<fstruct.xdim; i++) {
        for (j=0; j<fstruct.ydim; j++) {
            F(i,j) = stateZero2D(fstruct.x[i],fstruct.y[j],x0,y0);
            //F(i,j) = stateZeroOne(fstruct.x[i],fstruct.y[j],fac);
            //F(i,j) = 1/sqrt(2.)*(stateZero(fstruct.x[i],-CENTRE_POU,0)*stateOne(fstruct.y[j],0)
            //        + cexp(2*M_PI*I*fac)*stateZero(fstruct.y[j],0,0)*stateOne(fstruct.x[i],-CENTRE_POU));
            //F(i,j) = (stateZero(fstruct.x[i],x0,0)*stateOne(fstruct.y[j]));
        }
    }
    //print_qwave2D(stdout,f,&fstruct);
    
    fprintf(stderr,"Calculating time evolution...\n");
    //CrNi2Dexact(r,f,V,V,U,&fstruct);
    //CrNi2D_tr(r,f,V,V,U,&fstruct,0);
    CrNi2D_wf(f,V,V,U,&fstruct,stdout);
    fprintf(stderr,"Done.\n");


    //CrNi2Dexact(r,f,V,V,U,&fstruct,0);


    free(fstruct.x); free(fstruct.y); free(f);

    return 0;
}

double V(double x, double t){
    return 0.5*x*x;
}

double U (double x, double y, double t){
    return 0;
}

complex double stateZero2D(double x, double y, double x0, double y0) {
    return 1./sqrt(M_PI)*exp(-0.5*((x-x0)*(x-x0)+(y-y0)*(y-y0)));
}

complex double stateZeroOne(double x, double y, double fase) {
    return stateZero(x,0)*stateOne(y,0)+exp(I*2*M_PI*fase)*stateZero(y,0)*stateOne(x,0);
}

complex double stateZero(double x, double x0) {
    complex double fx = exp(-((x-x0)*(x-x0))*0.5);
    return 1./sqrt(sqrt(M_PI))*fx;
}

complex double stateOne(double x, double x0) {
    complex double fx = sqrt(2.)*(x-x0)*exp(-((x-x0)*(x-x0))*0.5);
    return sqrt(sqrt(1/M_PI))*fx;
}

#undef F

#undef QUANT_h
#undef QUANT_m
#undef QUANT_w
