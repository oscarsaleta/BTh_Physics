#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>

#include "numerical_lib.h"
#include "quantum_lib.h"
#include "functions.h"

#define LIMIT 8.0
#define GRID_DIM 251

#define F1(i,j) f1[(fstruct.xdim*(j)+(i))]
#define F2(i,j) f2[(fstruct.xdim*(j)+(i))]
#define Ft(i,j) ft[(fstruct.xdim*(j)+(i))]

#define X(i) (-LIMIT+(i)*fstruct.dx)
#define Y(j) (-LIMIT+(j)*fstruct.dy)

int main(int argc, char* argv[]){

    /* Program parameters */
    State fstruct;
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
    fstruct = struct_init(LIMIT,GRID_DIM,tmax,dt);
    fprintf(stderr,"xdim=%d ydim=%d lim=%g tmax=%g t=%g dt=%g dx=%g dy=%g\n",
            fstruct.xdim,fstruct.ydim,fstruct.limit,fstruct.t_max,*fstruct.t,
            *fstruct.dt,fstruct.dx,fstruct.dy);

    /* Allocation of arrays and pointers*/
    f1 = (complex double *)malloc(fstruct.xdim*fstruct.ydim*sizeof(complex double)); assert(f1 != NULL);
    f2 = (complex double *)malloc(fstruct.xdim*fstruct.ydim*sizeof(complex double)); assert(f2 != NULL);
    ft = (complex double *)malloc(fstruct.xdim*fstruct.ydim*sizeof(complex double)); assert(ft != NULL);


    /* Initial State */
    for (i=0; i<fstruct.xdim; i++) {
        for (j=0; j<fstruct.ydim; j++) {
            F1(i,j) = donut1(X(i),Y(j));
            F2(i,j) = donut3(X(i),Y(j));
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

    fprintf(stderr,"CALCULATING TRAJECTORY: r=(%g,%g)...",r[0],r[1]);
    CrNi2D_tr(r,ft,V,V,U,&fstruct,0);
    fprintf(stderr," Done\n");


    free(f1); free(f2); free(ft);

    return 0;
}


#undef F

#undef QUANT_h
#undef QUANT_m
#undef QUANT_w

