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

#define QUANT_h 1//1.05457173e-34
#define QUANT_m 1//1.67262178e-27
#define QUANT_w 1
#define CENTRE_POU 1.5

#define F(i,j) f[(fstruct.xdim*(j)+(i))]

double V(double x, double t);
double Vx(double x, double t);
double Vy(double x, double t);
double U(double x, double y, double t);

complex double stateZero2D(double x, double y, double x0, double y0);
complex double stateZero(double x, double x0);
complex double stateOne(double x, double x0);
complex double stateTwo(double x, double x0);

int main(int argc, char* argv[]){

    /* Program parameters */
    State fstruct;
    double xi=-8.,xf=8.,yi=-8.,yf=8.;
    complex double *f,*ft;
    double x0=-3., y0=0.;
    double dt;
    double tmax;
    double r[2];
    double fac; //factor de la fase
    /* Dummy variables */
    int i,j;
    int k;

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
    ft = (complex double *)malloc(fstruct.xdim*fstruct.ydim*sizeof(complex double)); assert(ft != NULL);


    /* Initiate state structure */
    *fstruct.t = 0;
    *fstruct.dt = dt;
    fstruct.t_max = tmax;
    vec_create(xi,xf,yi,yf,&fstruct);

    /* Initial State */ 
    for (i=0; i<fstruct.xdim; i++) {
        for (j=0; j<fstruct.ydim; j++) {
            /* |00> */
            F(i,j) = stateZero2D(fstruct.x[i],fstruct.y[j],-CENTRE_POU,0);

            /* |01>+|10> */
            //F(i,j) = 1/sqrt(2.)*(stateZero(fstruct.x[i],-CENTRE_POU)*stateOne(fstruct.y[j],0)
            //        + cexp(2*M_PI*I*fac)*stateZero(fstruct.y[j],0)*stateOne(fstruct.x[i],-CENTRE_POU));

            /* |01> */
            //F(i,j) = (stateZero(fstruct.x[i],x0)*stateOne(fstruct.y[j]));

            /* LG(2,0) */
            /*F(i,j) = -0.5*stateZero(fstruct.x[i],CENTRE_POU)*stateTwo(fstruct.y[j],0)
                    + I/sqrt(2.)*stateOne(fstruct.x[i],CENTRE_POU)*stateOne(fstruct.y[j],0)
                    + 0.5*stateZero(fstruct.y[j],0)*stateTwo(fstruct.x[i],CENTRE_POU);*/
        }
    }

    //print_qwave2D(stdout,f,&fstruct);


    //CrNi2D_wf(f,Vx,V,U,&fstruct,stdout);
    
    //fprintf(stderr,"CALCULATING TRAJECTORY: phase=%g r=%g\n",fac,r[1]);
    CrNi2D_tr(r,f,Vx,V,U,&fstruct,0);
    //CrNi2Dexact(r,f,V,V,U,&fstruct,0);

    /*for (k=0; k<(int)(tmax/dt); k++) {
        wf_static01(f,ft,k*dt,&fstruct);
        print_qwave2D(stdout,ft,&fstruct);
        *(fstruct.t)+=dt;
    }*/


    free(fstruct.x); free(fstruct.y); free(f);

    return 0;
}

double V(double x, double t){
    return 0.5*x*x;
}

double Vx (double x, double t){
    if (x<0)
        return 0.5*QUANT_m*QUANT_w*QUANT_w*(x+CENTRE_POU)*(x+CENTRE_POU);
    return 0.5*QUANT_m*QUANT_w*QUANT_w*(x-CENTRE_POU)*(x-CENTRE_POU);
}

double Vy (double x, double t){
    return 0.5*QUANT_m*QUANT_w*QUANT_w*x*x;
}
double U (double x, double y, double t){
    return 0;
}

complex double stateZero2D(double x, double y, double x0, double y0) {
    return stateZero(x,x0)*stateZero(y,y0);
    complex double fx, fy;
    fx = exp(-QUANT_m*QUANT_w*((x-x0)*(x-x0))/(2.*QUANT_h));//*cexp(I*px0*x/QUANT_h);
    fy = exp(-QUANT_m*QUANT_w*((y-y0)*(y-y0))/(2.*QUANT_h));//*cexp(I*py0*y/QUANT_h);
    return sqrt(QUANT_m*QUANT_w/(M_PI*QUANT_h))*fx*fy;
}

complex double stateZero(double x, double x0) {
    complex double fx;
    double alpha = QUANT_m*QUANT_w/QUANT_h;
    fx = exp(-alpha*(x-x0)*(x-x0)/2);
    return sqrt(sqrt(alpha/M_PI))*fx;
}

complex double stateOne(double x, double x0) {
    double alpha = QUANT_w*QUANT_m/QUANT_h;
    double y = sqrt(alpha)*(x-x0);
    complex double fx = sqrt(2.)*y*exp(-(y*y)/2.);
    return sqrt(sqrt(alpha/M_PI))*fx;
}

complex double stateTwo(double x, double x0) {
    double alpha = sqrt(QUANT_w*QUANT_m/QUANT_h);
    return sqrt(alpha/(sqrt(M_PI)*8))*(4*alpha*alpha*x*x-2)*exp(-0.5*alpha*alpha*x*x);
}

#undef F

#undef QUANT_h
#undef QUANT_m
#undef QUANT_w
