#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define QUANT_h 1//1.05457173e-34
#define QUANT_m 1//1.67262178e-27
#define QUANT_w 1
#define CENTRE_POU 0
#define IM_MAXIT 1200
//#define GRID_DIM 501

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
complex double donut3(double x, double y);

