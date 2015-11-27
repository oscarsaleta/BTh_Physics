/* Structure for storing the state of the system */
typedef struct state_struct {
    int xdim, ydim;
    double limit;
    double t_max, *t, *dt;
    double dx, dy;
} State;

void tridag(complex double *a, complex double *b, complex double *c, complex double *r, complex double *u, complex double *gam, int n);

double interpol(double* x, double* fx, int N, int j, double x0, double delta);

double interpol2D (double *f,  int i, int j, double x0, double y0, State *prm);

complex double simpson2d (complex double *f, State *prm);

State struct_init (double limit, double dim, double tmax, double dt);
