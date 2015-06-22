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
    int i,j,nt;
    State fstruct;
    double xi=-8.,xf=8.,yi=-8.,yf=8.;
    complex double *f;
    double x0=0,y0=0;
    double r[2]={0.,1.};
    double dt=0.05, tmax=1.00;
    /* File for reading */
    FILE *input,*output;

    /* INITIALIZE STRUCTURE AND STATE*/

    fstruct.xdim = fstruct.ydim = 16;

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

    /* Initial State */ 
    for (i=0; i<fstruct.xdim; i++) {
        for (j=0; j<fstruct.ydim; j++) {
            //F(i,j) = stateZero2D(fstruct.x[i],fstruct.y[j],x0,y0,0,0);
            F(i,j) = 1/sqrt(2.)*(stateZero(fstruct.x[i],x0,0)*stateOne(fstruct.y[j],0)
                    + I*stateZero(fstruct.y[j],y0,0)*stateOne(fstruct.x[i],0));
            //F(i,j) = (stateZero(fstruct.x[i],x0,0)*stateOne(fstruct.y[j]));
        }
    }

    output = fopen("copytest.dat","w");
    CrNi2D_wf(r,f,V,V,U,&fstruct,0,output);
    fclose(output);

    *(fstruct.t) = 0;
    input = fopen("copytest.dat","r");
    nt = ceil(tmax/dt);
    fprintf(stderr,"%d\n",nt);
    for (i=0; i<=nt; i++) {
        //rewind(input);
        read_qwave2D(input,f,&fstruct);
        print_qwave2D(stdout,f,&fstruct);
        *(fstruct.t)+=*(fstruct.dt);
    }

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
