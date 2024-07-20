#include "interpolation.h"
//#include "function.h"
#include "ndarray.h"
#include <assert.h>
#include <fstream>
#include <iostream>

#define SAVE_ARRAYS true

void product(NdArray &m1, NdArray &m2, NdArray &r)
{
    assert(m1.shape[0] == m2.shape[1]);
    assert(m1.shape[0] == r.shape[0]);
    assert(m2.shape[1] == r.shape[1]);

    for (int i = 0; i < r.shape[0]; i++) {
        for (int j = 0; j < r.shape[1]; j++) {
            double val = 0;
            for (int k = 0; k < m1.shape[0]; k++) {
                double m1ik = m1.at(i, k);
                double m2kj = m2.at(k, j);
                val += m1ik * m2kj;
            }
            r.at(i, j) = val;
        }
    }
}

void initializeG(NdArray &G, double *x, double *y, int nx, int ny,
                 double(f)(double, double), double(dfdx)(double, double),
                 double(dfdy)(double, double), double(d2fdxdy)(double, double),
                 double error)
//                 double error, int n_grid_x, int n_grid_y)
{
    // Создаем массивы (матрицы), необходимые для вычисления коэффициентов
    NdArray Ax(NULL, P, P);
    NdArray Ay(NULL, P, P);
    NdArray Tij(NULL, P, P);
    NdArray Fij(NULL, P, P);

    //    std::ofstream out("/home/vlatse/Temp/G00.txt");
    //    std::ofstream outf("/home/vlatse/Temp/F00.txt"); std::ofstream
    //    outax("/home/vlatse/Temp/Ax00.txt"); std::ofstream
    //    outay("/home/vlatse/Temp/Ay00.txt"); std::ofstream
    //    outt("/home/vlatse/Temp/T00.txt");

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double xi = x[i];
            double yj = y[j];
            double xip1 = x[i + 1];
            double yjp1 = y[j + 1];

            double dx = xip1 - xi;
            double dx2 = dx * dx;
            double dx3 = dx2 * dx;

            double dy = yjp1 - yj;
            double dy2 = dy * dy;
            double dy3 = dy2 * dy;

            Ax.at(0, 0) = 1;
            Ax.at(0, 1) = 0;
            Ax.at(0, 2) = 0;
            Ax.at(0, 3) = 0;

            Ax.at(1, 0) = 0;
            Ax.at(1, 1) = 1;
            Ax.at(1, 2) = 0;
            Ax.at(1, 3) = 0;

            Ax.at(2, 0) = -3 / dx2;
            Ax.at(2, 1) = -2 / dx;
            Ax.at(2, 2) = 3 / dx2;
            Ax.at(2, 3) = -1 / dx;

            Ax.at(3, 0) = 2 / dx3;
            Ax.at(3, 1) = 1 / dx2;
            Ax.at(3, 2) = -2 / dx3;
            Ax.at(3, 3) = 1 / dx2;

            Ay.at(0, 0) = 1;
            Ay.at(1, 0) = 0;
            Ay.at(2, 0) = 0;
            Ay.at(3, 0) = 0;

            Ay.at(0, 1) = 0;
            Ay.at(1, 1) = 1;
            Ay.at(2, 1) = 0;
            Ay.at(3, 1) = 0;

            Ay.at(0, 2) = -3 / dy2;
            Ay.at(1, 2) = -2 / dy;
            Ay.at(2, 2) = 3 / dy2;
            Ay.at(3, 2) = -1 / dy;

            Ay.at(0, 3) = 2 / dy3;
            Ay.at(1, 3) = 1 / dy2;
            Ay.at(2, 3) = -2 / dy3;
            Ay.at(3, 3) = 1 / dy2;

            double fij = f(xi, yj);
            if ((i == nx / 2) & (j == ny / 2)) {
                fij += error;
            }

            double dfdyij = dfdy(xi, yj);
            double fijp1 = f(xi, yjp1);

            if ((i == nx / 2) & (j + 1 == ny / 2)) {
                fijp1 += error;
            }

            double fip1jp1 = f(xip1, yjp1);
            if ((i + 1 == nx / 2) & (j + 1 == ny / 2)) {
                fip1jp1 += error;
            }

            double fip1j = f(xip1, yj);
            if ((i + 1 == nx / 2) & (j == ny / 2)) {
                fip1j += error;
            }

            double dfdyijp1 = dfdy(xi, yjp1);
            Fij.at(0, 0) = fij;
            Fij.at(0, 1) = dfdyij;
            Fij.at(0, 2) = fijp1;
            Fij.at(0, 3) = dfdyijp1;

            Fij.at(1, 0) = dfdx(xi, yj);
            Fij.at(1, 1) = d2fdxdy(xi, yj);
            Fij.at(1, 2) = dfdx(xi, yjp1);
            Fij.at(1, 3) = d2fdxdy(xi, yjp1);

            Fij.at(2, 0) = fip1j;
            Fij.at(2, 1) = dfdy(xip1, yj);
            Fij.at(2, 2) = fip1jp1;
            Fij.at(2, 3) = dfdy(xip1, yjp1);

            Fij.at(3, 0) = dfdx(xip1, yj);
            Fij.at(3, 1) = d2fdxdy(xip1, yj);
            Fij.at(3, 2) = dfdx(xip1, yjp1);
            Fij.at(3, 3) = d2fdxdy(xip1, yjp1);

            NdArray Gij(&G.at(i, j), P, P);
            product(Ax, Fij, Tij);
            product(Tij, Ay, Gij);
            //            if(i == 0 && j == 1){
            //               out << Gij;
            //               outf << Fij;
            //               outax << Ax;
            //               outay << Ay;
            //               outt << Tij;
            //            }
        }
    }
}

void fillArray(double a, double b, double *x, int n)
{
    double dx = (b - a) / (n - 2);
    for (int i = 0; i < n; i++) {
        x[i] = a + i * dx;
    }
}

void initG(double *state, double (*f_test)(double, double),
           double (*df_test_dx)(double, double),
           double (*df_test_dy)(double, double),
           double (*d2f_test_dxdy)(double, double), int n, int m, double a,
           double b, double c, double d, double error)
{
    NdArray G(state, n, m, P, P);
    double *xdata = state + G.size;
    fillArray(a, b, xdata, n + 1);
    double *ydata = xdata + n + 1;
    fillArray(c, d, ydata, m + 1);
    initializeG(G, xdata, ydata, n, m, f_test, df_test_dx, df_test_dy,
                d2f_test_dxdy, error);
}

double compute(double *state, int nx, int ny, double x, double y)
{
    NdArray G(state, nx, ny, P, P);
    double *xdata = state + G.size;
    double *ydata = xdata + nx + 1;
    double xminVal = xdata[0];
    double yminVal = ydata[0];
    double dx = xdata[1] - xdata[0];
    double dy = ydata[1] - ydata[0];
    return eval(G, xminVal, yminVal, dx, dy, x, y);
}

double eval(NdArray &F, double xmin, double ymin, double dx, double dy,
            double x, double y)
{
    int px = (x - xmin) / dx;
    int py = (y - ymin) / dy;
    if (px < 0) {
        px = 0;
    }
    if (py < 0) {
        py = 0;
    }
    if (py >= F.shape[1])
        py = F.shape[1] - 1;
    if (px >= F.shape[0])
        px = F.shape[0] - 1;

    double deltax = x - px * dx - xmin;
    double deltay = y - py * dy - ymin;

    double dxPowK = 1;
    double result = 0;

    for (int k = 0; k < P; k++) {
        double dyPowL = 1;
        for (int l = 0; l < P; l++) {
            double gamma = F.at(px, py, k, l);
            result += gamma * dxPowK * dyPowL;
            dyPowL *= deltay;
        }
        dxPowK *= deltax;
    }
    return result;
}

std::ostream &operator<<(std::ostream &o, NdArray &a)
{
    bool is2d = (a.shape[2] == 1) && (a.shape[3] == 1);
    for (int i = 0; i < a.shape[0]; i++) {
        for (int j = 0; j < a.shape[1]; j++) {
            for (int k = 0; k < a.shape[2]; k++) {
                for (int l = 0; l < a.shape[3]; l++) {
                    o << a.at(i, j, k, l) << ' ';
                }
            }
            if (!is2d)
                o << '\n';
        }
        if (is2d)
            o << '\n';
    }
    return o;
}
