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
#define CENTRE_POU 0
#define IM_MAXIT 1200
#define GRID_DIM 251

#define F1(i,j) f1[(fstruct.xdim*(j)+(i))]
#define F2(i,j) f2[(fstruct.xdim*(j)+(i))]
#define Ft(i,j) ft[(fstruct.xdim*(j)+(i))]

double V(double x, double t);
double Vx(double x, double y, double t);
double Vy(double x, double y, double t);
double U(double x, double y, double t);
int isCenter(double x, double y, double epsilon);

complex double stateZero2D(double x, double y, double x0, double y0);
complex double stateZero(double x, double x0);
complex double stateOne(double x, double x0);
complex double stateTwo(double x, double x0);
complex double donut1(double x, double y);
complex double donut2(double x, double y);

int main(int argc, char* argv[]){

    /* Program parameters */
    State fstruct;
    double xi=-8.,xf=8.,yi=-8.,yf=8.;
    complex double *f1,*f2,*ft;
    double dt;
    double tmax;
    double r[2];
    double a,b;
    /* Dummy variables */
    int i,j;

    /* READING PARAMETERS */
    if (argc < 7
            || sscanf(argv[1],"%lf",&dt)!=1
            || sscanf(argv[2],"%lf",&tmax)!=1
            || sscanf(argv[3],"%lf",&a)!=1
            || sscanf(argv[4],"%lf",&b)!=1
            || sscanf(argv[5],"%lf",&r[0])!=1
            || sscanf(argv[6],"%lf",&r[1])!=1
            ) {
        fprintf(stderr,"%s: dt tmax a b x0 y0\n",argv[0]);
        return 1;
    }

    /* INITIALIZE STRUCTURE AND STATE*/

    fstruct.xdim = fstruct.ydim = GRID_DIM;

    /* Allocation of arrays and pointers*/
    fstruct.t = (double *)malloc(sizeof(double));
    fstruct.dt = (double *)malloc(sizeof(double));
    fstruct.x = (double *)malloc(fstruct.xdim*sizeof(double)); assert(fstruct.x != NULL);
    fstruct.y = (double *)malloc(fstruct.ydim*sizeof(double)); assert(fstruct.y != NULL);
    f1 = (complex double *)malloc(fstruct.xdim*fstruct.ydim*sizeof(complex double)); assert(f1 != NULL);
    f2 = (complex double *)malloc(fstruct.xdim*fstruct.ydim*sizeof(complex double)); assert(f2 != NULL);
    ft = (complex double *)malloc(fstruct.xdim*fstruct.ydim*sizeof(complex double)); assert(ft != NULL);


    /* Initiate state structure */
    *fstruct.t = 0;
    *fstruct.dt = dt;
    fstruct.t_max = tmax;
    vec_create(xi,xf,yi,yf,&fstruct);

    /* Initial State */
    {
        for (i=0; i<fstruct.xdim; i++) {
        #pragma omp parallel num_threads(4)
        {
            #pragma omp for
            for (j=0; j<fstruct.ydim; j++) {
                /* |00> */
                //F(i,j) = stateZero2D(fstruct.x[i],fstruct.y[j],-CENTRE_POU,0);
                F1(i,j) = donut1(fstruct.x[i],fstruct.y[j]);
                F2(i,j) = donut2(fstruct.x[i],fstruct.y[j]);
            }
        }
        }
    }
    /* Imaginary time evolution */
    fprintf(stderr,"CALCULATING INITIAL WAVE FUNCTION...");
    if (a!=0)
        CrNi2D_im_wf(f1,V,V,U,&fstruct,IM_MAXIT);
    if (b!=0)
        CrNi2D_im_wf(f2,V,V,U,&fstruct,IM_MAXIT);
    for (i=0; i<fstruct.xdim; i++) {
        for (j=0; j<fstruct.ydim; j++) {
            Ft(i,j) = a*F1(i,j)+b*F2(i,j);
        }
    }
    //print_qwave2D(stdout,ft,&fstruct);
    fprintf(stderr," Done\n");

    fprintf(stderr,"CALCULATING WAVE FUNCTION EVOLUTION...:");
    CrNi2D_wf(ft,V,V,U,&fstruct,stdout);
    
    //fprintf(stderr,"CALCULATING TRAJECTORY: r=(%g,%g)...",r[0],r[1]);
    //CrNi2D_tr(r,ft,V,V,U,&fstruct,0);
    fprintf(stderr," Done\n");


    free(fstruct.x); free(fstruct.y); free(f1); free(f2); free(ft);

    return 0;
}

double V (double x, double t) {
    return 0;
}

double U (double x, double y, double t){
    if (isCenter(x,y,1e-5))
        return -1e5;
    return -1/(0.5*sqrt(x*x+y*y));
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

complex double donut1 (double x, double y) {
    return 1/sqrt(2.)*( stateZero(x,0)*stateOne(y,0)+I*stateOne(x,0)*stateZero(y,0) );
}
complex double donut2 (double x, double y) {
    return 1/sqrt(2.)*( stateZero(x,0)*stateOne(y,0)-I*stateOne(x,0)*stateZero(y,0) );
}

complex double stateTwo(double x, double x0) {
    double alpha = sqrt(QUANT_w*QUANT_m/QUANT_h);
    return sqrt(alpha/(sqrt(M_PI)*8))*(4*alpha*alpha*x*x-2)*exp(-0.5*alpha*alpha*x*x);
}

int isCenter(double x, double y, double epsilon) {
    return (x*x+y*y < epsilon) ? 1 : 0;
}
            

#undef F

#undef QUANT_h
#undef QUANT_m
#undef QUANT_w
