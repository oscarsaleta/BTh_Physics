/*NUMERICAL LIB*/
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include "numerical_lib.h"

/* TRIDIAGONAL MATRIX METHOD */
/* Given a system with diagonal b, upper diagonal c,
 * and lower diagonal a, with independet term r, it
 * solves the system [ a b c ] * u = r, for u.
 *
 * IMPORTANT NOTE: 
 * since a and c have one less component a[0] and c[n-1]
 * can have any value, because they won't be used.
 * Values for a go from 1 to n-1 and c from 0 to n-2.
 */ 
void tridag(complex double *a, complex double *b, complex double *c, complex double *r, complex double *u, complex double *gam, int n){
    int j;
    complex double bet;


    if (b[0] == 0.0){
        fprintf(stderr,"Error fatal!\n");
        system("kill -2");
    }

    u[0] = r[0]/(bet=b[0]);
    for (j=1;j<n;j++) {
        gam[j]=c[j-1]/bet;
        bet=b[j]-a[j]*gam[j];
        if(bet == 0.0){
            fprintf(stderr,"Error fatal!\n");
            system("kill -2");
        }
        u[j]=(r[j]-a[j]*u[j-1])/bet;
    }
    for (j=(n-2);j>=0;j--)
        u[j] -= gam[j+1]*u[j+1];

}


/* QUADRATIC 1D INTERPOLATION */
/* Interpolates a given function in order to evaluate it at
 * point x0 between two x array points.
 * Variables:
 *  x = spatial vector
 *  fx = function array of lenght 3 (only 3 points of fx are needed for interpolation)
 *  N = length of x
 *  j = index of closest point of x to x0
 *  x0 = point where f has to be evaluated
 *  delta = spatial step
 * */
double interpol(double* x, double* fx, int N, int j, double x0, double delta){

    return fx[0]+(fx[1]-fx[0])/delta*(x0-x[j-1])+
               (fx[2]-2*fx[1]+fx[0])/(2.*delta*delta)*(x0-x[j-1])*(x0-x[j]);

}


/* BILINEAL 2D INTERPOLATION */
/* Interpolates a 2D function using the bilineal method in order to evaluate
 * it at (x0,y0).
 * Variables:
 *  f = vector f of length 4 (only 4 points needed for bilineal interpolation)
 *    f[0] = f(0,0)  (Q_11)
 *    f[1] = f(1,0)  (Q_21)
 *    f[2] = f(0,1)  (Q_12)
 *    f[3] = f(1,1)  (Q_22)
 *  i = closest x index to x0 (floor)
 *  j = closest y index to y0 (floor)
 *  x0, y0 = point where f has to be evaluated
 *  prm = structure containing state parameters
 * */
double interpol2D (double *f, int i, int j, double x0, double y0, State *prm) {
    double *x = prm->x;
    double *y = prm->y;
    double dx = prm->dx; 
    double dy = prm->dy;

    return 1./(dx*dy)*(f[0]*(x[i+1]-x0)*(y[j+1]-y0)+f[1]*(x0-x[i])*(y[j+1]-y0)
            +f[2]*(x[i+1]-x0)*(y0-y[j])+f[3]*(x0-x[i])*(y0-y[j]));
}


#define F(i,j) f[(xdim*(j)+(i))]
/* SIMPSON METHOD FOR 2d INTEGRATION */
/* Integrates numerically a function given f and steps dx dy.
 * Since we need an odd number of points, and x,f have an 
 * even number of points (usually in QM), we neglect,f(i,n) and f(m,j),
 * since they are  0 anyways because of absorbing boundries .
 * Variables:
 * f = complex image of x,y equally spaced in x and y
 * prm = structure containing state parameters
 */
complex double simpson2d (complex double *f, State *prm){
    int xdim = prm->xdim;//(*((struct state_struct*)prm)).xdim;
    int ydim = prm->ydim;//(*((struct state_struct*)prm)).ydim;
    double dx = prm->dx;//(*((struct state_struct*)prm)).dx; 
    double dy = prm->dy;//(*((struct state_struct*)prm)).dy;
    complex double aux, sum = F(0,0)+F(0,ydim-2)+F(xdim-2,0)+F(xdim-2,ydim-2);
    int i, j;

    aux = 0;
    for(j=1;j<=(ydim-2)/2;j++)
        aux += F(0,2*j-1);
    sum += 4.*aux;


    aux = 0;
    for(j = 1;j <= (ydim-2)/2; j++)
        aux += F(xdim-2,2*j-1);
    sum += 4.*aux;

    aux = 0;
    for(j = 1; j <= (ydim-2)/2-1; j++)
        aux += F(0,2*j);
    sum += 2.*aux;

    
    aux = 0;
    for(j = 1; j <= (ydim-2)/2-1; j++)
        aux += F(xdim-2,2*j);
    sum += 2.*aux;
    
    /**/

    aux = 0;
    for(j=1;j<=(xdim-2)/2;j++)
        aux += F(2*j-1,0);
    sum += 4.*aux;


    aux = 0;
    for(j = 1;j <= (xdim-2)/2; j++)
        aux += F(2*j-1,ydim-2);
    sum += 4.*aux;

    aux = 0;
    for(j = 1; j <= (xdim-2)/2-1; j++)
        aux += F(2*j,0);
    sum += 2.*aux;

    
    aux = 0;
    for(j = 1; j <= (xdim-2)/2-1; j++)
        aux += F(2*j,ydim-2);
    sum += 2.*aux;
    
    /**/

    aux = 0;
    for(j = 1; j <= (ydim-2)/2; j++)
        for(i=1; i <= (xdim-2)/2; i++)
            aux += F(2*i-1,2*j-1);
    sum += 16.*aux;


    aux = 0;
    for(j = 1; j <= (ydim-2)/2-1; j++)
        for(i=1; i <= (xdim-2)/2; i++)
            aux += F(2*i-1,2*j);
    sum += 8.*aux;

    aux = 0;
    for(j = 1; j <= (ydim-2)/2; j++)
        for(i=1; i <= (xdim-2)/2-1; i++)
            aux += F(2*i,2*j-1);
    sum += 8.*aux;

    aux = 0;
    for(j = 1; j <= (ydim-2)/2-1; j++)
        for(i=1; i <= (xdim-2)/2-1; i++)
            aux += F(2*i,2*j);
    sum += 4.*aux;

    return sum*dx*dy/9.;
} 
#undef F


/* VEC_CREATE */
/* Creates vectors of x and y for 2d utility.
 * Variables:
 * x0, xf = boundaries of x
 * y0, yf = boundaries of y
 * prm = structure containing state parameters
 */
void vec_create (double x0, double xf, double y0, double yf, State *prm){
    int i;
    int xdim = prm->xdim;
    int ydim = prm->ydim;
    
    prm->dx = (xf-x0)/(xdim-1);
    prm->dy = (yf-y0)/(ydim-1);
    
    for(i=0; i<xdim ; i++)
            prm->x[i] = x0+i*prm->dx;

    for(i=0; i<ydim ; i++)
            prm->y[i] = y0+i*prm->dy;
} 
