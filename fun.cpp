#include "fun.h"
double f(double x, double y)
{

    return x * x + y * y;
}

double df_x(double x, double y)
{
    y++;
    return 2 * x;
}

double df_y(double x, double y)
{
    x++;
    return 2 * y;
}
double df_xy(double x, double y)
{
    x++;
    y++;
    return 0;
}
