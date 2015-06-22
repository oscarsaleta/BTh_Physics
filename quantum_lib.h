void print_qwave2D (FILE *output, complex double *f, State *state);

void read_qwave2D (FILE *input, complex double *f, State *state);



/* Absorbing Boundary Conditions */
/* Variables:
 *  N = length of x
 *  perc = percentage of x to absorbe in each side
 *  nu = power of the sinus
 * */
void abs_mask(int N, double perc, int nu, double *mask);


/* Bohm velocity of a particle in 2d*/
/* Variables:
 *  v = vector (vx,vy)
 *  psi = wavefunction (xdim x ydim)
 *  position = position to find the velocities (px,py)
 *  prm = structure containing state parameters 
 */
void bohm_vel2D(double *v, complex double *psi,  double *position, void *prm);


void CrNi2D_im_wf (complex double *psi, double (*Vx)(double, double), double (*Vy)(double, double), double (*U)(double, double, double), State *prm);

/* Wavefunction time evolution (Crank-Nicolson exact)*/
/* Variables:
 *  psi = wavefunction (m x n)
 *  Vx(x,t) = potential energy for dimension x
 *  Vy(y,t) = potential energy for dimension y
 *  U(x,y,t) = interaction potential energy
 *  prm = pointer to struct with program parameters
 *  flag = 0 if potential doesn't depend on time
 *         1 if potential depends on time
 *  output = file to write data (can be stdout)
 */
void CrNi2D_wf (complex double *psi, double (*Vx)(double, double), double (*Vy)(double, double), double (*U)(double, double, double), State *prm, FILE *output);


void wf_static01 (complex double *psia, complex double *psib, double t, State *prm);

    
/* Trajectory time evolution (Crank-Nicolson exact + adaptative Runge-Kutta)*/
/* Variables:
 *  psi = wavefunction (m x n)
 *  Vx(x,t) = potential energy for dimension x
 *  Vy(y,t) = potential energy for dimension y
 *  U(x,y,t) = interaction potential energy
 *  flag = 0 if potential doesn't depend on time
 *         1 if potential depends on time
 *  prm = pointer to struct with program parameters
 */
void CrNi2Dexact (double *pos, complex double *psi, double (*Vx)(double, double), double (*Vy)(double, double), double (*U)(double, double, double), void *prm, int flag);



void CrNi2D_tr (double *pos, complex double *psi, double (*Vx)(double, double), double (*Vy)(double, double), double (*U)(double, double, double), State *prm, int flag);


/* Bohm trajectories for a particle in 2D (using rk45 method)*/
/* Variables:
 *  r = vector (length 2) of position (x,y) of next point in the trajectory
 *  x,y = spatial grids (m x n)
 *  phi1 = wavefunction at time t (m x n)
 *  phi2 = wavefunction at time t+dt (m x n)
 *  phi3 = wavefunction at time t+2*dt (m x n)
 *  n = x length
 *  m = y length
 *  x0 = point in x (not necessarily of the grid) where trajectory has to be computed
 *  y0 = point in y (not necessarily of the grid) where trajectory has to be computed
 *  dt = halve of the time step (for runge kutta)
 */
void bohm_traj2D(double *r, complex double* phi1, complex double* phi2, complex double* phi3, State *prm);



void trajectoryFromFile(double *r, FILE *input, FILE *output, State *prm);
void Qj (complex double *psi, double *r, int nt, double *Q, State *prm); 
