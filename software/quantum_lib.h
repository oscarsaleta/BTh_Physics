void print_qwave2D (FILE *output, complex double *f, State *state);

void read_qwave2D (FILE *input, complex double *f, State *state);

void abs_mask(int N, double perc, int nu, double *mask);

void bohm_vel2D(double *v, complex double *psi,  double *position, State *prm);

void CrNi2D_wf (complex double *psi, double (*Vx)(double, double), double (*Vy)(double, double), double (*U)(double, double, double), State *prm, FILE *output);

void CrNi2D_im_wf (complex double *psi, double (*Vx)(double, double), double (*Vy)(double, double), double (*U)(double, double, double), State *prm, int maxit);

void CrNi2D_tr (double *pos, complex double *psi, double (*Vx)(double, double), double (*Vy)(double, double), double (*U)(double, double, double), State *prm, int flag);

void bohm_traj2D(double *r, complex double* phi1, complex double* phi2, complex double* phi3, State *prm);

void trajectoryFromFile(double *r, FILE *input, FILE *output, State *prm);
double Qj (complex double *psi, double *r, State *prm); 
