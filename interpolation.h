#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#define P 4
#include "ndarray.h"
#include <iostream>

void initG(double *state, double (*f_test)(double, double),
           double (*df_test_dx)(double, double),
           double (*df_test_dy)(double, double),
           double (*d2f_test_dxdy)(double, double), int n, int m, double a,
           double b, double c, double d, double error);

double compute(double *state, int nx, int ny, double x, double y);
double eval(NdArray &F, double xmin, double ymin, double dx, double dy,
            double x, double y);

std::ostream &operator<<(std::ostream &o, NdArray &a);

void product(NdArray &m1, NdArray &m2, NdArray &r);
void initializeG(NdArray &G, double *x, double *y, int nx, int ny,
                 double(f)(double, double), double(dfdx)(double, double),
                 double(dfdy)(double, double), double(d2fdxdy)(double, double),
                 double error);
//                 double error, int n_grid_x, int n_grid_y);
void fillArray(double a, double b, double *x, int n);

#endif // INTERPOLATION_H
