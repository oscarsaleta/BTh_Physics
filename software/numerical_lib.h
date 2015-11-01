/* Structure for storing the state of the system */
typedef struct state_struct {
    int xdim, ydim;
    double t_max, *t, *dt;
    double *x, *y;
    double dx, dy;
} State;

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
void tridag(complex double *a, complex double *b, complex double *c, complex double *r, complex double *u, complex double *gam, int n);


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
double interpol(double* x, double* fx, int N, int j, double x0, double delta);


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
double interpol2D (double *f,  int i, int j, double x0, double y0, State *prm);


/* SIMPSON METHOD FOR 2d INTEGRATION */
/* Integrates numerically a function given f and steps dx dy.
 * Since we need an odd number of points, and x,f have an 
 * even number of points (usually in QM), we neglect,f(i,n) and f(m,j),
 * since they are  0 anyways because of absorbing boundries .
 * Variables:
 * f = complex image of x,y equally spaced in x and y
 * prm = structure containing state parameters
 */
complex double simpson2d (complex double *f, State *prm);


/* VEC_CREATE */
/* Creates vectors of x and y for 2d utility.
 * Variables:
 * x0, xf = boundries of x
 * y0, yf = boundries of y
 * prm = structure containing state parameters
 */
void vec_create (double x0, double xf, double y0, double yf, State *prm);
