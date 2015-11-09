#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "functions.h"
#include "numerical_lib.h"
#include "quantum_lib.h"



double V (double x, double t) {
    return 0;
}

double U (double x, double y, double t){
    if (isCenter(x,y,1e-5))
        return -1e5;
    return -2./(sqrt(x*x+y*y));
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
complex double donut3 (double x, double y) {
    return 1/sqrt(2.*M_PI)*(x+I*y)*(x+I*y)*cexp(-0.5*(x*x+y*y));
}

complex double stateTwo(double x, double x0) {
    double alpha = sqrt(QUANT_w*QUANT_m/QUANT_h);
    return sqrt(alpha/(sqrt(M_PI)*8))*(4*alpha*alpha*x*x-2)*exp(-0.5*alpha*alpha*x*x);
}

int isCenter(double x, double y, double epsilon) {
    return (x*x+y*y < epsilon) ? 1 : 0;
}
            

