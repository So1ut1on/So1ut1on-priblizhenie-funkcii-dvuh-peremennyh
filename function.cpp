#include "function.h"

const double C = 1.;

double expStable(double x);

const double EPS = 1e-15;
double f1_test(double x, double y)
{
    return (x * x + y * y) / C;
}

// double df1_dx(double x, double y)
//{
//    y++;
//    return 2 * x / C;
//}

// double df1_dy(double x, double y)
//{
//    x++;
//    return 2 * y / C;
//}

// double d2f1_dxdy(double x, double y)
//{
//    x++;
//    y++;
//    return 0;
//}

double expStable(double x)
{
    if (x < -50) {
        return 0.;
    }
    return std::exp(x);
}

double f1(double x, double y)
{
    return 1;
}

double df1_dx(double x, double y)
{
    return 0;
}

double df1_dy(double x, double y)
{
    return 0;
}

double d2f1_dxdy(double x, double y)
{
    return 0;
}

double f_test(double x, double y)
{
    //    int sig = 2*(x > 0) - 1;
    //    double s2  = sig*sin(x)/(std::abs(x) + EPS);
    double val = y * y + x * x + EPS;
    double r = expStable(-(val));
    return r;
}

double df_test_dx(double x, double y)
{
    //    return (x* cos(x) - sin(x))/(x*x + EPS);

    double val = y * y + x * x + EPS;
    double r = expStable(-(val));
    return -r * 2 * x;
}

double df_test_dy(double x, double y)
{

    double val = y * y + x * x + EPS;
    double r = expStable(-(val));
    return -r * 2 * y;
}

double d2f_test_dxdy(double x, double y)
{
    double val = y * y + x * x + EPS;
    double r = expStable(-(val));
    return -r * 4 * x * y;
}

double f2(double x, double y)
{
    return x;
}

double df2_dx(double x, double y)
{
    return 1;
}

double df2_dy(double x, double y)
{
    return 0;
}

double d2f2_dxdy(double x, double y)
{
    return 0;
}

double f3(double x, double y)
{
    return y;
}

double df3_dx(double x, double y)
{
    return 0;
}

double df3_dy(double x, double y)
{
    return 1;
}

double d2f3_dxdy(double x, double y)
{
    return 0;
}

double f4(double x, double y)
{
    return x + y;
}

double df4_dx(double x, double y)
{
    return 1;
}

double df4_dy(double x, double y)
{
    return 1;
}

double d2f4_dxdy(double x, double y)
{
    return 0;
}

double f5(double x, double y)
{
    return sqrt(f6(x, y));
}

double df5_dx(double x, double y)
{
    return x / (f5(x, y) + EPS);
}

double df5_dy(double x, double y)
{
    return y / (f5(x, y) + EPS);
}

double d2f5_dxdy(double x, double y)
{
    double f5_val = f5(x, y);
    int sign = ((int)(f5_val >= 0) * 2 - 1);
    double eps = sign * EPS;
    return -2 * x * y / (f5_val * f5_val * f5_val + eps);
}

double f6(double x, double y)
{
    return x * x + y * y;
}

double df6_dx(double x, double y)
{
    return 2 * x;
}

double df6_dy(double x, double y)
{
    return 2 * y;
}

double d2f6_dxdy(double x, double y)
{
    return 0;
}

double f7(double x, double y)
{
    double v = x * x - y * y;
    return expStable(v);
}

double df7_dx(double x, double y)
{
    return f7(x, y) * x * 2;
}

double df7_dy(double x, double y)
{
    return -f7(x, y) * y * 2;
}

double d2f7_dxdy(double x, double y)
{
    return -4 * x * y * f7(x, y);
}

double f8(double x, double y)
{
    double v = x * x + y * y;
    return 1 / (25 * v + 1);
}

double df8_dx(double x, double y)
{
    double f_val = f8(x, y);
    return -50 * x * f_val * f_val;
}

double df8_dy(double x, double y)
{
    double f_val = f8(x, y);
    return -50 * y * f_val * f_val;
}

double d2f8_dxdy(double x, double y)
{
    double f_val = f8(x, y);
    return 5000 * x * y * f_val * f_val * f_val;
}
