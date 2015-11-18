/* QUANTUM_LIB */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include "numerical_lib.h"
#include "quantum_lib.h"

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#define QUANT_h 1//1.05457173e-34
#define QUANT_m 1//1.67262178e-27
#define QUANT_w 1


/* PRINT_QWAVE 2D */
/* Prints a 2D wavefunction into file 
 * with columns x,y,real(f(x,y)),imag(f(x,y)),abs(f(x,y))^2,V
 * Variables:
 *  name = name of output file
 *  f = wave function
 *  prm =  structur containing state parameters
 */
#define F(i,j) f[(xdim*(j)+(i))]
#define X(i) (-(state->limit)+(i)*(state->dx))
#define Y(i) (-(state->limit)+(i)*(state->dy))
void print_qwave2D (FILE *output, complex double *f, State *state){
    int i,j;
    int xdim = state->xdim;
    int ydim = state->ydim;

    fprintf(output,"# t %10.6G\n",*(state->t));

    for(i=0; i<xdim; i++) {
        for (j=0; j<ydim; j++) {
            fprintf(output,"%15.11G %15.11G %15.11G %15.11G %15.11G\n",
                X(i),Y(j),cabs(F(i,j))*cabs(F(i,j)),creal(F(i,j)),
                cimag(F(i,j)));
        }
        fprintf(output,"\n");
    }
    fprintf(output,"\n");
}
#undef F
#undef X
#undef Y


#define F(i,j) f[xdim*(j)+(i)]
void read_qwave2D (FILE *input, complex double *f, State *state) {
    int i,j,k,dummy=0;
    int xdim=state->xdim;
    int ydim=state->ydim;
    double aux[5];
    char buffer[100],*str;

    fgets(buffer,100,input);
    for (i=0; i<xdim; i++) {
        for (j=0; j<ydim; j++) {
            // agafem una línia, després llegim valors del buffer
            fgets(buffer,100,input);
            str=buffer;
            for (k=0; k<5; k++) {
                aux[k]=strtod(str,&str);
            }
            F(i,j)=aux[3]+I*aux[4];
            dummy++;
        }
        // saltar una línia
        fgets(buffer,100,input);
    }
    // saltar una línia
    fgets(buffer,100,input);
    
}
#undef F


/* Absorbing Boundary Conditions */
/* Variables:
 *  N = length of x
 *  perc = percentage of x to absorbe in each side
 *  nu = power of the sinus
 * */
void abs_mask(int N, double perc, int nu, double *mask) {
    int d; /* Points that are absorbed in each side */
    int i;
    
    d = (int) floor(N*perc);

    for (i=0; i<d; i++) {
        mask[i] = pow(sin(M_PI*i/(2.*d)),nu);
        mask[N-1-i] = mask[i];
    }

    for (i=d; i<N-1-d; i++) {
        mask[i] = 1;
    }
}



/* Bohm velocity of a particle in 2d*/
/* Variables:
 *  v = vector (vx,vy)
 *  psi = wavefunction (xdim x ydim)
 *  position = position to find the velocities (px,py)
 *  prm = structure containing state parameters 
 */
#define PSI(i,j) psi[(xdim*(j)+(i))]
#define X(i) (-(prm->limit)+(i)*dx)
#define Y(i) (-(prm->limit)+(i)*dy)
void bohm_vel2D(double *v, complex double *psi,  double *position, State *prm){

    double dx = prm->dx;
    double dy = prm->dy;
//  double *x = prm->x;
//  double *y = prm->y;
    int xdim = prm->xdim;
    int ydim = prm->ydim;
    double vel[4];
    complex double psi_aux[3], dpsi;
    int i,j,k;

    /* Outside the grid speed is zero */
    if ( position[0]<X(0)
            || position[0]>X(xdim-1)
            || position[1]<Y(0)
            || position[1]>Y(ydim-1) ) {
        v[0] = v[1] = 0;
        return;
    }

    i = (int) floor((position[0]-X(0))/dx);
    j = (int) floor((position[1]-Y(0))/dy);

    /* Velocity in X-direction */

    /* Vx(0,0) */

    for (k=0; k<3; k++){
        if( i+k-1 < 0 || i+k-1 > xdim-1 )
            psi_aux[k] = 0;
        else
            psi_aux[k] = PSI(i+k-1,j);
    }
    
    if (cabs(psi_aux[1]) > 1e-30){
        dpsi = (psi_aux[2]-psi_aux[0])/(2.*dx);
        vel[0] = QUANT_h/QUANT_m*cimag(dpsi/psi_aux[1]);
    }else
        vel[0] = 0;

    /* Vx(1,0) */

    for(k=0;k<3;k++){
        if( i+k < 0 || i+k > xdim-1 )
            psi_aux[k] = 0;
        else
            psi_aux[k] = PSI(i+k,j);
    }
    
    if(cabs(psi_aux[1]) > 1e-30){
        dpsi = (psi_aux[2]-psi_aux[0])/(2.*dx);
        vel[1] = QUANT_h/QUANT_m*cimag(dpsi/psi_aux[1]);
    }else
        vel[1] = 0;

    /* Vx(0,1) */
    
    for(k=0;k<3;k++){
        if( i+k-1 < 0 || i+k-1 > xdim-1 || j+1 > ydim-1)
            psi_aux[k] = 0;
        else
            psi_aux[k] = PSI(i+k-1,j+1);
    }
    
    if(cabs(psi_aux[1]) > 1e-30){
        dpsi = (psi_aux[2]-psi_aux[0])/(2.*dx);
        vel[2] = QUANT_h/QUANT_m*cimag(dpsi/psi_aux[1]);
    }else
        vel[2] = 0;

    /* Vx(1,1) */
    
    for(k=0;k<3;k++){
        if( i+k < 0 || i+k > xdim-1 || j+1 > ydim-1)
            psi_aux[k] = 0;
        else
            psi_aux[k] = PSI(i+k,j+1);
    }
    
    if(cabs(psi_aux[1]) > 1e-30){
        dpsi = (psi_aux[2]-psi_aux[0])/(2.*dx);
        vel[3] = QUANT_h/QUANT_m*cimag(dpsi/psi_aux[1]);
    }else
        vel[3] = 0;

    /* Velocity X in (x0,y0) */
    
    v[0] = interpol2D(vel,i,j,position[0],position[1],prm);

    /* Velocity in Y-direction */

    /* Vy(0,0) */

    for(k=0;k<3;k++){
        if( j+k-1 < 0 || j+k-1 > ydim-1 )
            psi_aux[k] = 0;
        else
            psi_aux[k] = PSI(i,j+k-1);
    }
    
    if(cabs(psi_aux[1]) > 1e-30){
        dpsi = (psi_aux[2]-psi_aux[0])/(2.*dy);
        vel[0] = QUANT_h/QUANT_m*cimag(dpsi/psi_aux[1]);
    }else
        vel[0] = 0;

    /* Vy(1,0) */

    for(k=0;k<3;k++){
        if( j+k-1 < 0 || j+k-1 > ydim-1  || i+1 > xdim-1)
            psi_aux[k] = 0;
        else
            psi_aux[k] = PSI(i+1,j+k-1);
    }
    
    if(cabs(psi_aux[1]) > 1e-30){
        dpsi = (psi_aux[2]-psi_aux[0])/(2.*dy);
        vel[1] = QUANT_h/QUANT_m*cimag(dpsi/psi_aux[1]);
    }else
        vel[1] = 0;

    /* Vy(0,1) */
    
    for(k=0;k<3;k++){
        if( j+k < 0 || j+k > ydim-1 )
            psi_aux[k] = 0;
        else
            psi_aux[k] = PSI(i,j+k);
    }
    
    if(cabs(psi_aux[1]) > 1e-30){
        dpsi = (psi_aux[2]-psi_aux[0])/(2.*dy);
        vel[2] = QUANT_h/QUANT_m*cimag(dpsi/psi_aux[1]);
    }else
        vel[2] = 0;

    /* Vy(1,1) */
    
    for(k=0;k<3;k++){
        if( j+k < 0 || j+k > ydim-1 || i+1 > xdim-1)
            psi_aux[k] = 0;
        else
            psi_aux[k] = PSI(i+1,j+k);
    }
    
    if(cabs(psi_aux[1]) > 1e-30){
        dpsi = (psi_aux[2]-psi_aux[0])/(2.*dy);
        vel[3] = QUANT_h/QUANT_m*cimag(dpsi/psi_aux[1]);
    }else
        vel[3] = 0;

    /* Velocity Y in (x0,y0) */
    
    v[1] = interpol2D(vel,i,j,position[0],position[1],prm);


} 
#undef PSI
#undef X
#undef Y


/* Time evolution (Crank-Nicolson exact) prints wavefunctions*/
/* Variables:
 *  pos = 
 *  psi = wavefunction (m x n)
 *  Vx(x,t) = potential energy for dimension x
 *  Vy(y,t) = potential energy for dimension y
 *  U(x,y,t) = interaction potential energy
 *  prm = pointer to struct with program parameters
 *  output = file to write data (can be stdout)
 */
#define PSI(i,j) psi[(xdim*(j)+(i))]
#define CPSI(i,j) current_psi[(xdim*(j)+(i))]
#define ABSPSISQ(i,j) abspsisq[(xdim*(j)+(i))]
#define VX(i,j) Varray_x[(xdim*(j)+(i))]
#define VY(i,j) Varray_y[(xdim*(j)+(i))]
#define BETAX(i,j) beta_x[(xdim*(j)+(i))]
#define BETAY(i,j) beta_y[(xdim*(j)+(i))]
#define X(i) (-(prm->limit)+(i)*dx)
#define Y(i) (-(prm->limit)+(i)*dy)
void CrNi2D_wf (complex double *psi, double (*Vx)(double, double), double (*Vy)(double, double), double (*U)(double, double, double), State *prm, FILE *output) {
    /*Reading the struct*/
//  double *x = prm->x;
//  double *y = prm->y;
    int xdim = prm->xdim;
    int ydim = prm->ydim;
    double dx = prm->dx;
    double dy = prm->dy;
    double dt = *(prm->dt);
    double *t = prm->t;
    double t_max = prm->t_max;
    /*Variable declaration*/
    int i,j,k=0; /*Loop variables*/
    complex double alpha_x,alpha_y; /*Alpha coefficients of Crank-Nicolson method*/
    complex double *beta_x,*beta_y; /*Beta coefficients of Crank-Nicolson method*/
    complex double *A_x,*A_y,*B; /* No C component defined because A=C in tridiag resolution */
    complex double *r,*gamma; /*Coefficients for triangular system resolution*/
    double *Varray_x, *Varray_y; /*Value of x and y components of V*/
    complex double *abspsisq, *psi_aux; /*Dummy psi storage*/
    double normalisation; /*For normalizing psi*/

    int max = xdim < ydim ? ydim : xdim; /*Maximum lenght between x and y*/

    beta_x = (complex double *)malloc(xdim*ydim*sizeof(complex double)); assert(beta_x != NULL);
    beta_y = (complex double *)malloc(xdim*ydim*sizeof(complex double)); assert(beta_y != NULL);
    A_x = (complex double *)malloc(xdim*sizeof(complex double)); assert(A_x != NULL);
    A_y = (complex double *)malloc(ydim*sizeof(complex double)); assert(A_y != NULL);
    B = (complex double *)malloc(max*sizeof(complex double)); assert(B != NULL);
    r = (complex double *)malloc(max*sizeof(complex double)); assert(r != NULL);
    Varray_x = (double *)malloc(xdim*ydim*sizeof(double)); assert(Varray_x != NULL);
    Varray_y = (double *)malloc(xdim*ydim*sizeof(double)); assert(Varray_y != NULL);
    abspsisq = (complex double *)malloc(xdim*ydim*sizeof(complex double)); assert(abspsisq != NULL);
    psi_aux = (complex double *)malloc(max*sizeof(complex double)); assert(psi_aux != NULL);
    gamma = (complex double *)malloc(max*sizeof(complex double)); assert(gamma != NULL);


    /* Initial potential (it may not vary over time) so we calculate it first */
    for (i=0; i<xdim; i++) {
        for (j=0; j<ydim; j++) {
            VX(i,j) = 0.5*(*Vx)(X(i),0) + 0.25*(*U)(X(i),Y(j),0);
            VY(i,j) = (*Vy)(Y(j),0) + 0.5*(*U)(X(i),Y(j),0);
        }
    }

    /* Print wavefunction */
    print_qwave2D(output,psi,prm);
    
    /* Let the main loop begin!
     *
     * This loop solves the Schrödinger equation for the system using
     *  the Crank-Nicolson method, and computes 6 psi at different times.
     * These psi's are used to take a step in the particles trajectories using
     *  an adaptative Runge Kutta 6th order implicit method with adaptative time-step.
     *  If the method fails to reach the desired accuracy, the CN and RK are repeated
     *  with a new dt (inner loop).
     * The main loop runs over the whole time domain of the problem.
     * */
    while (*t < t_max) {
        fprintf(stderr,"wf:: t = %10.6G\n",*t);

        alpha_x = I*QUANT_h*dt/(8*QUANT_m*dx*dx); /*We're applying only 1/2 of H_x in each x loop*/
        alpha_y = I*QUANT_h*dt/(4*QUANT_m*dy*dy);
        for (i=0; i<xdim; i++)
            A_x[i] = -alpha_x;
        
        for (j=0; j<ydim; j++)
            A_y[j] = -alpha_y;


        for (i=0; i<xdim; i++) {
            for (j=0; j<ydim; j++) {
                BETAX(i,j) = I*dt*VX(i,j)/(2.*QUANT_h);
                BETAY(i,j) = I*dt*VY(i,j)/(2.*QUANT_h);
            }
        }

        /* Evolve 1/2 H_x */
        for (j=0; j<ydim; j++) {
            for (i=0; i<xdim; i++) 
                psi_aux[i] = PSI(i,j);
            for (i=1; i<xdim-1; i++)
                r[i] = (1-2*alpha_x-BETAX(i,j))*psi_aux[i]+alpha_x*(psi_aux[i-1]+psi_aux[i+1]);
            r[0] = (1-2*alpha_x-BETAX(0,j))*psi_aux[0]+alpha_x*psi_aux[1];
            r[xdim-1] = (1-2*alpha_x-BETAX(xdim-1,j))*psi_aux[xdim-1]+alpha_x*psi_aux[xdim-2];
            
            for(i=0;i<xdim;i++){
                BETAX(i,j) = I*dt*VX(i,j)/(2.*QUANT_h);
                B[i] = 1+2*alpha_x+BETAX(i,j);
            }
            /* Solve the tridiagonal system */
            tridag(A_x,B,A_x,r,psi_aux,gamma,xdim);
            /* Reassign partially evolved psi */
            for (i=0; i<xdim; i++)
                PSI(i,j) = psi_aux[i];
        }


        /* Normalization*/
        for (i=0; i<xdim; i++)
            for (j=0; j<ydim; j++)
                ABSPSISQ(i,j) = cabs(PSI(i,j))*cabs(PSI(i,j));

        normalisation = sqrt(simpson2d(abspsisq,prm));

        for (i=0; i<xdim; i++)
            for (j=0; j<ydim; j++)
                PSI(i,j) /= normalisation;

        
        /* Evolve H_y */
        for (i=0; i<xdim; i++) {
            for (j=0; j<ydim; j++) 
                psi_aux[j] = PSI(i,j); 
            
            for (j=1; j<ydim-1; j++)
                r[j] = (1-2*alpha_y-BETAY(i,j))*psi_aux[j]+alpha_y*(psi_aux[j-1]+psi_aux[j+1]);
            r[0] = (1-2*alpha_y-BETAY(i,0))*psi_aux[0]+alpha_y*psi_aux[1];
            r[ydim-1] = (1-2*alpha_y-BETAY(i,ydim-1))*psi_aux[ydim-1]+alpha_y*psi_aux[ydim-2];
            
            
            for (j=0; j<ydim; j++){
                BETAY(i,j) = I*dt*VY(i,j)/(2.*QUANT_h);
                B[j] = 1+2*alpha_y+BETAY(i,j);
            }
            /* Solve the tridiagonal system */
            tridag(A_y,B,A_y,r,psi_aux,gamma,ydim);
            /* Reassign partially evolved psi */
            for (j=0; j<ydim; j++)
                PSI(i,j) = psi_aux[j];
        }

        /* Normalization*/
        for (i=0; i<xdim; i++)
            for (j=0; j<ydim; j++)
                ABSPSISQ(i,j) = cabs(PSI(i,j))*cabs(PSI(i,j));

        normalisation = sqrt(simpson2d(abspsisq,prm));

        for (i=0; i<xdim; i++)
            for (j=0; j<ydim; j++)
                PSI(i,j) /= normalisation;
        
                
        /* Evolve 1/2 H_x */
        for (j=0; j<ydim; j++) {
            for (i=0; i<xdim; i++) 
                psi_aux[i] = PSI(i,j);
            for (i=1; i<xdim-1; i++)
                r[i] = (1-2*alpha_x-BETAX(i,j))*psi_aux[i]+alpha_x*(psi_aux[i-1]+psi_aux[i+1]);
            r[0] = (1-2*alpha_x-BETAX(0,j))*psi_aux[0]+alpha_x*psi_aux[1];
            r[xdim-1] = (1-2*alpha_x-BETAX(xdim-1,j))*psi_aux[xdim-1]+alpha_x*psi_aux[xdim-2];
            
            
            for(i=0;i<xdim;i++){
                BETAX(i,j) = I*dt*VX(i,j)/(2.*QUANT_h);
                B[i] = 1+2*alpha_x+BETAX(i,j);
            }
            /* Solve the tridiagonal system */
            tridag(A_x,B,A_x,r,psi_aux,gamma,xdim);
            /* Reassign partially evolved psi */
            for (i=0; i<xdim; i++)
                PSI(i,j) = psi_aux[i];
        }


        /* Normalization*/
        for (i=0; i<xdim; i++)
            for (j=0; j<ydim; j++)
                ABSPSISQ(i,j) = cabs(PSI(i,j))*cabs(PSI(i,j));

        normalisation = sqrt(simpson2d(abspsisq,prm));

        for (i=0; i<xdim; i++)
            for (j=0; j<ydim; j++)
                PSI(i,j) /= normalisation;


        /* Update time */
        *t += (dt);
        k++;
        /* Print wavefunction */
        if (k % 10 == 0)
            print_qwave2D(output,psi,prm);
        
    }

    free(A_x); free(A_y);
    free(B); free(r);
    free(gamma);
    free(beta_x); free(beta_y);
    free(Varray_x); free(Varray_y);
    free(abspsisq); free(psi_aux);
}
#undef PSI
#undef CPSI
#undef ABSPSISQ
#undef VX
#undef VY
#undef BETAX
#undef BETAY
#undef X
#undef Y


/* Imaginary Time evolution (Crank-Nicolson exact) */
/* Variables:
 *  pos = 
 *  psi = wavefunction (m x n)
 *  Vx(x,t) = potential energy for dimension x
 *  Vy(y,t) = potential energy for dimension y
 *  U(x,y,t) = interaction potential energy
 *  prm = pointer to struct with program parameters
 */
#define PSI(i,j) psi[(xdim*(j)+(i))]
#define CPSI(i,j) current_psi[(xdim*(j)+(i))]
#define ABSPSISQ(i,j) abspsisq[(xdim*(j)+(i))]
#define VX(i,j) Varray_x[(xdim*(j)+(i))]
#define VY(i,j) Varray_y[(xdim*(j)+(i))]
#define BETAX(i,j) beta_x[(xdim*(j)+(i))]
#define BETAY(i,j) beta_y[(xdim*(j)+(i))]
#define X(i) (-(prm->limit)+(i)*dx)
#define Y(i) (-(prm->limit)+(i)*dy)
void CrNi2D_im_wf (complex double *psi, double (*Vx)(double, double), double (*Vy)(double, double), double (*U)(double, double, double), State *prm, int maxit) {
    /*Reading the struct*/
//  double *x = prm->x;
//  double *y = prm->y;
    int xdim = prm->xdim;
    int ydim = prm->ydim;
    double dx = prm->dx;
    double dy = prm->dy;
    double dt = *(prm->dt);
    /*Variable declaration*/
    int i,j,k; /*Loop variables*/
    complex double alpha_x,alpha_y; /*Alpha coefficients of Crank-Nicolson method*/
    complex double *beta_x,*beta_y; /*Beta coefficients of Crank-Nicolson method*/
    complex double *A_x,*A_y,*B; /* No C component defined because A=C in tridiag resolution */
    complex double *r,*gamma; /*Coefficients for triangular system resolution*/
    double *Varray_x, *Varray_y; /*Value of x and y components of V*/
    complex double *abspsisq, *psi_aux; /*Dummy psi storage*/
    double normalisation; /*For normalizing psi*/

    int max = xdim < ydim ? ydim : xdim; /*Maximum lenght between x and y*/

    beta_x = (complex double *)malloc(xdim*ydim*sizeof(complex double)); assert(beta_x != NULL);
    beta_y = (complex double *)malloc(xdim*ydim*sizeof(complex double)); assert(beta_y != NULL);
    A_x = (complex double *)malloc(xdim*sizeof(complex double)); assert(A_x != NULL);
    A_y = (complex double *)malloc(ydim*sizeof(complex double)); assert(A_y != NULL);
    B = (complex double *)malloc(max*sizeof(complex double)); assert(B != NULL);
    r = (complex double *)malloc(max*sizeof(complex double)); assert(r != NULL);
    Varray_x = (double *)malloc(xdim*ydim*sizeof(double)); assert(Varray_x != NULL);
    Varray_y = (double *)malloc(xdim*ydim*sizeof(double)); assert(Varray_y != NULL);
    abspsisq = (complex double *)malloc(xdim*ydim*sizeof(complex double)); assert(abspsisq != NULL);
    psi_aux = (complex double *)malloc(max*sizeof(complex double)); assert(psi_aux != NULL);
    gamma = (complex double *)malloc(max*sizeof(complex double)); assert(gamma != NULL);


    /* Initial potential (it may not vary over time) so we calculate it first */
    for (i=0; i<xdim; i++) {
        for (j=0; j<ydim; j++) {
            VX(i,j) = 0.5*(*Vx)(X(i),0) + 0.25*(*U)(X(i),Y(j),0);
            VY(i,j) = (*Vy)(Y(j),0) + 0.5*(*U)(X(i),Y(j),0);
        }
    }

    
    /* Let the main loop begin!
     *
     * This loop solves the Schrödinger equation for the system using
     *  the Crank-Nicolson method, and computes 6 psi at different times.
     * These psi's are used to take a step in the particles trajectories using
     *  an adaptative Runge Kutta 6th order implicit method with adaptative time-step.
     *  If the method fails to reach the desired accuracy, the CN and RK are repeated
     *  with a new dt (inner loop).
     * The main loop runs over the whole time domain of the problem.
     * */
    for (k=0; k<maxit; k++) {

        alpha_x = QUANT_h*dt/(8*QUANT_m*dx*dx); /*We're applying only 1/2 of H_x in each x loop*/
        alpha_y = QUANT_h*dt/(4*QUANT_m*dy*dy);
        for (i=0; i<xdim; i++)
            A_x[i] = -alpha_x;
        for (j=0; j<ydim; j++)
            A_y[j] = -alpha_y;

        for (i=0; i<xdim; i++) {
            for (j=0; j<ydim; j++) {
                BETAX(i,j) = dt*VX(i,j)/(2.*QUANT_h);
                BETAY(i,j) = dt*VY(i,j)/(2.*QUANT_h);
            }
        }
        /* Evolve 1/2 H_x */
        for (j=0; j<ydim; j++) {
            for (i=0; i<xdim; i++) 
                psi_aux[i] = PSI(i,j);
            for (i=1; i<xdim-1; i++)
                r[i] = (1-2*alpha_x-BETAX(i,j))*psi_aux[i]+alpha_x*(psi_aux[i-1]+psi_aux[i+1]);
            r[0] = (1-2*alpha_x-BETAX(0,j))*psi_aux[0]+alpha_x*psi_aux[1];
            r[xdim-1] = (1-2*alpha_x-BETAX(xdim-1,j))*psi_aux[xdim-1]+alpha_x*psi_aux[xdim-2];
            for(i=0;i<xdim;i++){
                BETAX(i,j) = dt*VX(i,j)/(2.*QUANT_h);
                B[i] = 1+2*alpha_x+BETAX(i,j);
            }
            /* Solve the tridiagonal system */
            tridag(A_x,B,A_x,r,psi_aux,gamma,xdim);
            /* Reassign partially evolved psi */
            for (i=0; i<xdim; i++)
                PSI(i,j) = psi_aux[i];
        }


        /* Normalization*/
        for (i=0; i<xdim; i++)
            for (j=0; j<ydim; j++)
                ABSPSISQ(i,j) = cabs(PSI(i,j))*cabs(PSI(i,j));

        normalisation = sqrt(simpson2d(abspsisq,prm));

        for (i=0; i<xdim; i++)
            for (j=0; j<ydim; j++)
                PSI(i,j) /= normalisation;

        
        /* Evolve H_y */
        for (i=0; i<xdim; i++) {
            for (j=0; j<ydim; j++) 
                psi_aux[j] = PSI(i,j); 
            
            for (j=1; j<ydim-1; j++)
                r[j] = (1-2*alpha_y-BETAY(i,j))*psi_aux[j]+alpha_y*(psi_aux[j-1]+psi_aux[j+1]);
            r[0] = (1-2*alpha_y-BETAY(i,0))*psi_aux[0]+alpha_y*psi_aux[1];
            r[ydim-1] = (1-2*alpha_y-BETAY(i,ydim-1))*psi_aux[ydim-1]+alpha_y*psi_aux[ydim-2];
            
            
            for (j=0; j<ydim; j++){
                BETAY(i,j) = dt*VY(i,j)/(2.*QUANT_h);
                B[j] = 1+2*alpha_y+BETAY(i,j);
            }
            /* Solve the tridiagonal system */
            tridag(A_y,B,A_y,r,psi_aux,gamma,ydim);
            /* Reassign partially evolved psi */
            for (j=0; j<ydim; j++)
                PSI(i,j) = psi_aux[j];
        }

        /* Normalization*/
        for (i=0; i<xdim; i++)
            for (j=0; j<ydim; j++)
                ABSPSISQ(i,j) = cabs(PSI(i,j))*cabs(PSI(i,j));

        normalisation = sqrt(simpson2d(abspsisq,prm));

        for (i=0; i<xdim; i++)
            for (j=0; j<ydim; j++)
                PSI(i,j) /= normalisation;
        
                
        /* Evolve 1/2 H_x */
        for (j=0; j<ydim; j++) {
            for (i=0; i<xdim; i++) 
                psi_aux[i] = PSI(i,j);
            for (i=1; i<xdim-1; i++)
                r[i] = (1-2*alpha_x-BETAX(i,j))*psi_aux[i]+alpha_x*(psi_aux[i-1]+psi_aux[i+1]);
            r[0] = (1-2*alpha_x-BETAX(0,j))*psi_aux[0]+alpha_x*psi_aux[1];
            r[xdim-1] = (1-2*alpha_x-BETAX(xdim-1,j))*psi_aux[xdim-1]+alpha_x*psi_aux[xdim-2];
            
            
            for(i=0;i<xdim;i++){
                BETAX(i,j) = dt*VX(i,j)/(2.*QUANT_h);
                B[i] = 1+2*alpha_x+BETAX(i,j);
            }
            /* Solve the tridiagonal system */
            tridag(A_x,B,A_x,r,psi_aux,gamma,xdim);
            /* Reassign partially evolved psi */
            for (i=0; i<xdim; i++)
                PSI(i,j) = psi_aux[i];
        }


        /* Normalization*/
        for (i=0; i<xdim; i++)
            for (j=0; j<ydim; j++)
                ABSPSISQ(i,j) = cabs(PSI(i,j))*cabs(PSI(i,j));

        normalisation = sqrt(simpson2d(abspsisq,prm));

        for (i=0; i<xdim; i++)
            for (j=0; j<ydim; j++)
                PSI(i,j) /= normalisation;


    }

    free(A_x); free(A_y);
    free(B); free(r);
    free(gamma);
    free(beta_x); free(beta_y);
    free(Varray_x); free(Varray_y);
    free(abspsisq); free(psi_aux);
}
#undef PSI
#undef CPSI
#undef ABSPSISQ
#undef VX
#undef VY
#undef BETAX
#undef BETAY
#undef X
#undef Y
//TODO: podem posar un flag a CrNi2D_wf perque posi una constant =I ò =1 (per im.time o normal)



/* Time evolution (Crank-Nicolson exact + Runge-Kutta 4)*/
/* Variables:
 *  psi = wavefunction (m x n)
 *  Vx(x,t) = potential energy for dimension x
 *  Vy(y,t) = potential energy for dimension y
 *  U(x,y,t) = interaction potential energy
 *  flag = 0 if potential doesn't depend on time
 *         1 if potential depends on time
 *  prm = pointer to struct with program parameters
 */
#define PSI(i,j) psi[(xdim*(j)+(i))]
#define CPSI(i,j) current_psi[(xdim*(j)+(i))]
#define ABSPSISQ(i,j) abspsisq[(xdim*(j)+(i))]
#define VX(i,j) Varray_x[(xdim*(j)+(i))]
#define VY(i,j) Varray_y[(xdim*(j)+(i))]
#define BETAX(i,j) beta_x[(xdim*(j)+(i))]
#define BETAY(i,j) beta_y[(xdim*(j)+(i))]
#define X(i) (-(prm->limit)+(i)*dx)
#define Y(i) (-(prm->limit)+(i)*dy)
void CrNi2D_tr (double *pos, complex double *psi, double (*Vx)(double, double), double (*Vy)(double, double), double (*U)(double, double, double), State *prm, int flag) {
    /*Reading the struct*/
//  double *x = prm->x;
//  double *y = prm->y;
    int xdim = prm->xdim;
    int ydim = prm->ydim;
    double dx = prm->dx;
    double dy = prm->dy;
    double dt = *(prm->dt);
    double *t = prm->t;
    double t_max = prm->t_max;

    /*Variable declaration*/
    int i,j,k; /*Loop variables*/
    complex double alpha_x,alpha_y; /*Alpha coefficients of Crank-Nicolson method*/
    complex double *beta_x,*beta_y; /*Beta coefficients of Crank-Nicolson method*/
    complex double *A_x,*A_y,*B; /* No C component defined because A=C in tridiag resolution */
    complex double *r,*gamma; /*Coefficients for triangular system resolution*/
    double *Varray_x, *Varray_y; /*Value of x and y components of V*/
    complex double *abspsisq, *psi_aux; /*Dummy psi storage*/
    double normalisation; /*For normalizing psi*/
    /*RK of 6th order requires 6 steps (ie. 6 psi's) to be known to perform a step in trajectories*/
    complex double *current_psi;
    complex double *psivec[2];
    double Q; /*quantum potential*/

    int max = xdim < ydim ? ydim : xdim; /*Maximum lenght between x and y*/

    /* Fourth order RK method */
    static double a2=0.5,a3=1.0;
    double dtcoef[]={a2,a3};
    /* These variables are used to do RK steps and track accuracy */
    double vel[2],ak2[2],ak3[2],ak4[2];
    double posTemp[2];


    beta_x = (complex double *)malloc(xdim*ydim*sizeof(complex double)); assert(beta_x != NULL);
    beta_y = (complex double *)malloc(xdim*ydim*sizeof(complex double)); assert(beta_y != NULL);
    psivec[0] = (complex double *)malloc(xdim*ydim*sizeof(complex double)); assert(psivec[0] != NULL);
    psivec[1] = (complex double *)malloc(xdim*ydim*sizeof(complex double)); assert(psivec[1] != NULL);
    A_x = (complex double *)malloc(xdim*sizeof(complex double)); assert(A_x != NULL);
    A_y = (complex double *)malloc(ydim*sizeof(complex double)); assert(A_y != NULL);
    B = (complex double *)malloc(max*sizeof(complex double)); assert(B != NULL);
    r = (complex double *)malloc(max*sizeof(complex double)); assert(r != NULL);
    Varray_x = (double *)malloc(xdim*ydim*sizeof(double)); assert(Varray_x != NULL);
    Varray_y = (double *)malloc(xdim*ydim*sizeof(double)); assert(Varray_y != NULL);
    abspsisq = (complex double *)malloc(xdim*ydim*sizeof(complex double)); assert(abspsisq != NULL);
    psi_aux = (complex double *)malloc(max*sizeof(complex double)); assert(psi_aux != NULL);
    gamma = (complex double *)malloc(max*sizeof(complex double)); assert(gamma != NULL);


    /* Initial potential (it may not vary over time) so we calculate it first */
    for (i=0; i<xdim; i++) {
        for (j=0; j<ydim; j++) {
            VX(i,j) = 0.5*(*Vx)(X(i),0) + 0.25*(*U)(X(i),Y(j),0);
            VY(i,j) = (*Vy)(Y(j),0) + 0.5*(*U)(X(i),Y(j),0);
        }
    }


    /* Let the main loop begin!
     *
     * This loop solves the Schrödinger equation for the system using
     *  the Crank-Nicolson method, and computes 2 psi at different times.
     * These psi's are used to take a step in the particles trajectories using
     *  a Runge Kutta method (4th order) with static time-step.
     * The main loop runs over the whole time domain of the problem.
     * */
    while (*t < t_max) {

        for(k=0;k<2;k++) {
            current_psi = psivec[k];

            alpha_x = I*QUANT_h*dt*dtcoef[k]/(8*QUANT_m*dx*dx); /*We're applying only 1/2 of H_x in each x loop*/
            alpha_y = I*QUANT_h*dt*dtcoef[k]/(4*QUANT_m*dy*dy);
            for (i=0; i<xdim; i++)
                A_x[i] = -alpha_x;
            
            for (j=0; j<ydim; j++)
                A_y[j] = -alpha_y;


            if(flag == 1){
                for (i=0; i<xdim; i++) {
                    for (j=0; j<ydim; j++) {
                        VX(i,j) = 0.5*(*Vx)(X(i),*t) + 0.25*(*U)(X(i),Y(j),*t);
                        VY(i,j) = (*Vy)(Y(j),*t) + 0.5*(*U)(X(i),Y(j),*t);
                    }
                }
            }

            for (i=0; i<xdim; i++) {
                for (j=0; j<ydim; j++) {
                    BETAX(i,j) = I*dt*dtcoef[k]*VX(i,j)/(2.*QUANT_h);
                    BETAY(i,j) = I*dt*dtcoef[k]*VY(i,j)/(2.*QUANT_h);
                }
            }

            /* Evolve 1/2 H_x */
            for (j=0; j<ydim; j++) {
                for (i=0; i<xdim; i++) 
                    psi_aux[i] = PSI(i,j);
                for (i=1; i<xdim-1; i++)
                    r[i] = (1-2*alpha_x-BETAX(i,j))*psi_aux[i]+alpha_x*(psi_aux[i-1]+psi_aux[i+1]);
                r[0] = (1-2*alpha_x-BETAX(0,j))*psi_aux[0]+alpha_x*psi_aux[1];
                r[xdim-1] = (1-2*alpha_x-BETAX(xdim-1,j))*psi_aux[xdim-1]+alpha_x*psi_aux[xdim-2];
                
                if(flag == 1)
                    for(i=0;i<xdim;i++)
                        VX(i,j) = 0.5*(*Vx)(X(i),*t+dt*0.5*dtcoef[k]) + 0.25*(*U)(X(i),Y(j),*t+dt*0.5*dtcoef[k]);
                
                for(i=0;i<xdim;i++){
                    BETAX(i,j) = I*dt*dtcoef[k]*VX(i,j)/(2.*QUANT_h);
                    B[i] = 1+2*alpha_x+BETAX(i,j);
                }
                /* Solve the tridiagonal system */
                tridag(A_x,B,A_x,r,psi_aux,gamma,xdim);
                /* Reassign partially evolved psi */
                for (i=0; i<xdim; i++)
                    CPSI(i,j) = psi_aux[i];
            }


            /* Normalization*/
            for (i=0; i<xdim; i++)
                for (j=0; j<ydim; j++)
                    ABSPSISQ(i,j) = cabs(CPSI(i,j))*cabs(CPSI(i,j));

            normalisation = sqrt(simpson2d(abspsisq,prm));

            for (i=0; i<xdim; i++)
                for (j=0; j<ydim; j++)
                    CPSI(i,j) /= normalisation;

            
            /* Evolve H_y */
            for (i=0; i<xdim; i++) {
                for (j=0; j<ydim; j++) 
                    psi_aux[j] = CPSI(i,j); 
                
                for (j=1; j<ydim-1; j++)
                    r[j] = (1-2*alpha_y-BETAY(i,j))*psi_aux[j]+alpha_y*(psi_aux[j-1]+psi_aux[j+1]);
                r[0] = (1-2*alpha_y-BETAY(i,0))*psi_aux[0]+alpha_y*psi_aux[1];
                r[ydim-1] = (1-2*alpha_y-BETAY(i,ydim-1))*psi_aux[ydim-1]+alpha_y*psi_aux[ydim-2];
                
                if(flag == 1)
                    for (j=0; j<ydim; j++)
                        VY(i,j) = (*Vy)(Y(j),*t+dt*dtcoef[k]) + 0.5*(*U)(X(i),Y(j),*t+dt*dtcoef[k]);
                
                for (j=0; j<ydim; j++){
                    BETAY(i,j) = I*dt*dtcoef[k]*VY(i,j)/(2.*QUANT_h);
                    B[j] = 1+2*alpha_y+BETAY(i,j);
                }
                /* Solve the tridiagonal system */
                tridag(A_y,B,A_y,r,psi_aux,gamma,ydim);
                /* Reassign partially evolved psi */
                for (j=0; j<ydim; j++)
                    CPSI(i,j) = psi_aux[j];
            }

            /* Normalization*/
            for (i=0; i<xdim; i++)
                for (j=0; j<ydim; j++)
                    ABSPSISQ(i,j) = cabs(CPSI(i,j))*cabs(CPSI(i,j));

            normalisation = sqrt(simpson2d(abspsisq,prm));

            for (i=0; i<xdim; i++)
                for (j=0; j<ydim; j++)
                    CPSI(i,j) /= normalisation;
            
            
            /* Evolve 1/2 H_x */
            for (j=0; j<ydim; j++) {
                for (i=0; i<xdim; i++) 
                    psi_aux[i] = CPSI(i,j);
                for (i=1; i<xdim-1; i++)
                    r[i] = (1-2*alpha_x-BETAX(i,j))*psi_aux[i]+alpha_x*(psi_aux[i-1]+psi_aux[i+1]);
                r[0] = (1-2*alpha_x-BETAX(0,j))*psi_aux[0]+alpha_x*psi_aux[1];
                r[xdim-1] = (1-2*alpha_x-BETAX(xdim-1,j))*psi_aux[xdim-1]+alpha_x*psi_aux[xdim-2];
                
                if(flag == 1)
                    for(i=0;i<xdim;i++)
                        VX(i,j) = 0.5*(*Vx)(X(i),*t+dt*dtcoef[k]) + 0.25*(*U)(X(i),Y(j),*t+dt*dtcoef[k]);
                
                for(i=0;i<xdim;i++){
                    BETAX(i,j) = I*dt*dtcoef[k]*VX(i,j)/(2.*QUANT_h);
                    B[i] = 1+2*alpha_x+BETAX(i,j);
                }
                /* Solve the tridiagonal system */
                tridag(A_x,B,A_x,r,psi_aux,gamma,xdim);
                /* Reassign partially evolved psi */
                for (i=0; i<xdim; i++)
                    CPSI(i,j) = psi_aux[i];
            }


            /* Normalization*/
            for (i=0; i<xdim; i++)
                for (j=0; j<ydim; j++)
                    ABSPSISQ(i,j) = cabs(CPSI(i,j))*cabs(CPSI(i,j));

            normalisation = sqrt(simpson2d(abspsisq,prm));

            for (i=0; i<xdim; i++)
                for (j=0; j<ydim; j++)
                    CPSI(i,j) /= normalisation;

        }

        /* RUNGE KUTTA STEPS */
        /* Compute first velocity at time t */
        bohm_vel2D(vel,psi,pos,prm); /*First step*/
        for (i=0;i<2;i++)
            posTemp[i]=pos[i]+vel[i]*dt*0.5;
        bohm_vel2D(ak2,psivec[0],posTemp,prm); /* Second step*/
        for (i=0;i<2;i++)
            posTemp[i]=pos[i]+dt*ak2[i]*0.5;
        bohm_vel2D(ak3,psivec[0],posTemp,prm); /*Third step*/
        for (i=0;i<2;i++) 
            posTemp[i]=pos[i]+dt*ak3[i];
        bohm_vel2D(ak4,psivec[1],posTemp,prm); /*Fourth step*/
        for (i=0;i<2;i++) /*Accumulate increments with proper weights.*/
            pos[i]+=dt*(vel[i]+2*ak2[i]+2*ak3[i]+ak4[i])/6.;
           
        /* Update time */
        *t += dt;
//      fprintf(stderr,"exact:: t = %10.6G\n",*t);

        /* Compute Qj */
        Q=Qj(psi,pos,prm);

        /* Update wavefunction */
        current_psi = psivec[1];
        for(i=0;i<xdim;i++)
            for(j=0;j<ydim;j++)
                PSI(i,j) = CPSI(i,j);
        /* Print position */
        fprintf(stdout, "%15.10g   %15.10g   %15.10g   %15.10G   %15.10G\n",pos[0],pos[1],*t,sqrt(pos[0]*pos[0]+pos[1]*pos[1]),Q);
        
    }

    free(A_x); free(A_y);
    free(B); free(r);
    free(psivec[0]); free(psivec[1]);
    free(gamma);
    free(beta_x); free(beta_y);
    free(Varray_x); free(Varray_y);
    free(abspsisq); free(psi_aux);
}
#undef PSI
#undef CPSI
#undef ABSPSISQ
#undef VX
#undef VY
#undef BETAX
#undef BETAY
#undef X
#undef Y



/* Bohm trajectories for a particle in 2D (using rk4 method)*/
/* Variables:
 *  r = vector (length 2) of position (x,y) of next point in the trajectory
 *  phi1 = wavefunction at time t (m x n)
 *  phi2 = wavefunction at time t+dt (m x n)
 *  phi3 = wavefunction at time t+2*dt (m x n)
 *  prm = state parameters
 */
void bohm_traj2D(double *r, complex double* phi1, complex double* phi2, complex double* phi3, State *prm) {
    double dt = *(prm->dt);
    double k1x, k2x, k3x, k4x;
    double k1y, k2y, k3y, k4y;
    double v[2];
    double posTemp[2];

    posTemp[0]=r[0];
    posTemp[1]=r[1];
    bohm_vel2D(v,phi1,posTemp,prm);
    k1x = 2*v[0]*dt;
    k1y = 2*v[1]*dt;

    posTemp[0]=r[0]+k1x/2.;
    posTemp[1]=r[1]+k1y/2.;
    bohm_vel2D(v,phi2,posTemp,prm);
    k2x = 2*v[0]*dt;
    k2y = 2*v[1]*dt;

    posTemp[0]=r[0]+k2x/2.;
    posTemp[1]=r[1]+k2y/2.;
    bohm_vel2D(v,phi2,posTemp,prm);
    k3x = 2*v[0]*dt;
    k3y = 2*v[1]*dt;

    posTemp[0]=r[0]+k3x;
    posTemp[0]=r[1]+k3y;
    bohm_vel2D(v,phi3,posTemp,prm);
    k4x = 2*v[0]*dt;
    k4y = 2*v[1]*dt;

    r[0] += (k1x+2*k2x+2*k3x+k4x)/6.;
    r[1] += (k1y+2*k2y+2*k3y+k4y)/6.; 
}


#define PSI(i,j) psi[(j)*xdim+(i)]
#define X(i) (-(prm->limit)+(i)*dx)
#define Y(i) (-(prm->limit)+(i)*dy)
double Qj (complex double *psi, double *r, State *prm) {
    int i,j;
//  double *x=prm->x;
//  double *y=prm->y;
    double dx=prm->dx;
    double dy=prm->dy;
    int xdim=prm->xdim;
    double Qaux[4];
    double constant=-(QUANT_h*QUANT_h)/(2.*QUANT_m);

    i=(int)floor((r[0]-X(0))/dx);
    j=(int)floor((r[1]-Y(0))/dy);

    /* Podem necessitar Qj a un punt que no és del grid.
     * Calculem Qj als 4 punts més propers i interpolem (2D)
     * r es el punt on volem calcular Qj
     * Usem diferencies finites centrades pel Laplacià
     * */

    /* (i,j) */
    Qaux[0]=constant * 1/PSI(i,j)*((PSI(i+1,j)-2*PSI(i,j)+PSI(i-1,j))/(dx*dx)+
            (PSI(i,j+1)-2*PSI(i,j)+PSI(i,j-1))/(dy*dy));
    /* (i+1,j) */
    Qaux[1]=constant * 1/PSI(i+1,j)*((PSI(i+2,j)-2*PSI(i+1,j)+PSI(i,j))/(dx*dx)+
            (PSI(i+1,j+1)-2*PSI(i+1,j)+PSI(i+1,j-1))/(dy*dy));
    /* (i,j+1) */
    Qaux[2]=constant * 1/PSI(i,j+1)*((PSI(i+1,j+1)-2*PSI(i,j+1)+PSI(i-1,j+1))/(dx*dx)+
            (PSI(i,j+2)-2*PSI(i,j+1)+PSI(i,j))/(dy*dy));
    /* (i+1,j+1) */
    Qaux[3]=constant * 1/PSI(i+1,j+1)*((PSI(i+2,j+1)-2*PSI(i+1,j+1)+PSI(i,j+1))/(dx*dx)+
            (PSI(i+1,j+2)-2*PSI(i+1,j+1)+PSI(i+1,j))/(dy*dy));

    return interpol2D(Qaux,i,j,r[0],r[1],prm);
}
#undef PSI
#undef X
#undef Y
