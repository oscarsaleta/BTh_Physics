#define SAFETY 0.9
#define ERRCON 1.89e-4 /* (5/SAFETY)^(1/PGROW) */
#define PGROW -0.2
#define PSHRNK -0.25
#define PSI(i,j) psi[(xdim*(j)+(i))]
#define CPSI(i,j) current_psi[(xdim*(j)+(i))]
#define ABSPSISQ(i,j) abspsisq[(xdim*(j)+(i))]
#define VX(i,j) Varray_x[(xdim*(j)+(i))]
#define VY(i,j) Varray_y[(xdim*(j)+(i))]
#define BETAX(i,j) beta_x[(xdim*(j)+(i))]
#define BETAY(i,j) beta_y[(xdim*(j)+(i))]
/* Time evolution (Crank-Nicolson exact + adaptative Runge-Kutta)*/
/* Variables:
 *  psi = wavefunction (m x n)
 *  Vx(x,t) = potential energy for dimension x
 *  Vy(y,t) = potential energy for dimension y
 *  U(x,y,t) = interaction potential energy
 *  flag = 0 if potential doesn't depend on time
 *         1 if potential depends on time
 *  prm = pointer to struct with program parameters
 */
void CrNi2D_tr (double *pos, complex double *psi, double (*Vx)(double, double), double (*Vy)(double, double), double (*U)(double, double, double), State *prm, int flag) {
    /*Reading the struct*/
    double *x = prm->x;
    double *y = prm->y;
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
    complex double *psivec[5];

    int max = xdim < ydim ? ydim : xdim; /*Maximum lenght between x and y*/

    /* Fourth order RK method */
    static double a2=0.5,a5=1.0;
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
            VX(i,j) = 0.5*(*Vx)(x[i],0) + 0.25*(*U)(x[i],y[j],0);
            VY(i,j) = (*Vy)(y[j],0) + 0.5*(*U)(x[i],y[j],0);
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
                        VX(i,j) = 0.5*(*Vx)(x[i],*t) + 0.25*(*U)(x[i],y[j],*t);
                        VY(i,j) = (*Vy)(y[j],*t) + 0.5*(*U)(x[i],y[j],*t);
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
                        VX(i,j) = 0.5*(*Vx)(x[i],*t+(*dt)*0.5*dtcoef[k]) + 0.25*(*U)(x[i],y[j],*t+dt*0.5*dtcoef[k]);
                
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
                        VY(i,j) = (*Vy)(y[j],*t+*dt*dtcoef[k]) + 0.5*(*U)(x[i],y[j],*t+*dt*dtcoef[k]);
                
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
                        VX(i,j) = 0.5*(*Vx)(x[i],*t+(*dt)*dtcoef[k]) + 0.25*(*U)(x[i],y[j],*t+dt*dtcoef[k]);
                
                for(i=0;i<xdim;i++){
                    BETAX(i,j) = I*(*dt)*dtcoef[k]*VX(i,j)/(2.*QUANT_h);
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
            posTemp[i]=pos[i]+vel[i]*dt;
        bohm_vel2D(ak2,psivec[0],posTemp,prm); /* Second step*/
        for (i=0;i<2;i++)
            posTemp[i]=pos[i]+dt*ak2[i];
        bohm_vel2D(ak3,psivec[0],posTemp,prm); /*Third step*/
        for (i=0;i<2;i++) 
            posTemp[i]=pos[i]+dt*ak3[i];
        bohm_vel2D(ak4,psivec[1],posTemp,prm); /*Fourth step*/
        for (i=0;i<2;i++) /*Accumulate increments with proper weights.*/
            pos[i]=pos[i]+3*dt*(vel[i]+2*ak2[i]+2*ak3[i]+ak4[i]);
           
        /* Update time */
        *t += dt;
//        fprintf(stderr,"exact:: t = %10.6G\n",*t);
        /* Update wavefunction */
        current_psi = psivec[1];
        for(i=0;i<xdim;i++)
            for(j=0;j<ydim;j++)
                PSI(i,j) = CPSI(i,j);
        /* Print position */
        fprintf(stdout, "%15.10g   %15.10g   %15.10g   %15.10G\n",pos[0],pos[1],*t,sqrt(pos[0]*pos[0]+pos[1]*pos[1]));
        
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

