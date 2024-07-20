#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "main.h"

#define eps 0.0000000001

double Pf(double x, double y, double a, double b, double c, double d, int n,
          int m, double *mf);
void Coeff(int n, int m, double *xarr, double *yarr, double *coef, double *Fij,
           double *A, double *Res);
double Value(double x, double y, int n, int m, double *coef, double *xarr,
             double *yarr, double *Temp, double *xxi, double *yyi);
int Init(int n_, double a_, double b_, int m_, double c_, double d_);
void Finalize(void);
void Input(int ch);
void Calc(int ch);
double f_xy(double x, double y);
