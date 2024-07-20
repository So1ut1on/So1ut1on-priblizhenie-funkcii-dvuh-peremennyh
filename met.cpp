#include "met.h"

static int n;
static double a;
static double b;

static int m;
static double c;
static double d;

static double *xarr = NULL;
static double *yarr = NULL;
static double *f_ij = NULL;
static double *Fij = NULL;
static double *coef = NULL;
static double *A = NULL;
static double *Res = NULL;
static double *Temp = NULL;
static double *xxi = NULL;
static double *yyi = NULL;

//начало общей части

int Init(int n_, double a_, double b_, int m_, double c_, double d_)
{
    n = n_;
    a = a_;
    b = b_;

    m = m_;
    c = c_;
    d = d_;

    xarr = (double *)malloc(n * sizeof(double));
    yarr = (double *)malloc(m * sizeof(double));
    f_ij = (double *)malloc(n * m * sizeof(double));
    Fij = (double *)malloc(16 * sizeof(double));
    A = (double *)malloc(16 * sizeof(double));
    Res = (double *)malloc(16 * sizeof(double));
    coef = (double *)malloc(16 * n * m * sizeof(double));
    Temp = (double *)malloc(16 * n * m * sizeof(double));
    xxi = (double *)malloc(4 * sizeof(double));
    yyi = (double *)malloc(4 * sizeof(double));

    // if (!(xarr && yarr && f_ij && coef)) return 0;

    return 1;
}

void Finalize(void)
{
    if (xarr) {
        free(xarr);
        xarr = NULL;
    }
    if (yarr) {
        free(yarr);
        yarr = NULL;
    }
    // if (f_ij) free(f_ij); f_ij = NULL;
    if (coef) {
        free(coef);
        coef = NULL;
    }
}

void Input(int ch)
{
    int i, j;
    double hx, hy;
    hx = (b - a) / (n - 1);
    // printf("hx: %lf \n", hx);
    hy = (d - c) / (m - 1);

    for (i = 0; i < n; i++) {
        xarr[i] = a + i * hx;
        for (j = 0; j < m; j++) {
            yarr[j] = c + j * hy;
            f_ij[i * m + j] = f(xarr[i], yarr[j], ch);
            // printf("x: %lf, y: %lf\n", xarr[i], yarr[j]);
        }
    }
}

static void calcFij(int i, int j, double *yarr, double *xarr, double *Fij,
                    int ch)
{
    Fij[0] = f(xarr[i], yarr[j], ch);
    Fij[1] = df_y(xarr[i], yarr[j], ch);
    Fij[2] = f(xarr[i], yarr[j + 1], ch);
    Fij[3] = df_y(xarr[i], yarr[j + 1], ch);

    Fij[4] = df_x(xarr[i], yarr[j], ch);
    Fij[5] = df_xy(xarr[i], yarr[j], ch);
    Fij[6] = df_x(xarr[i], yarr[j + 1], ch);
    Fij[7] = df_xy(xarr[i], yarr[j + 1], ch);

    Fij[8] = f(xarr[i + 1], yarr[j], ch);
    Fij[9] = df_y(xarr[i + 1], yarr[j], ch);
    Fij[10] = f(xarr[i + 1], yarr[j + 1], ch);
    Fij[11] = df_y(xarr[i + 1], yarr[j + 1], ch);

    Fij[12] = df_x(xarr[i + 1], yarr[j], ch);
    Fij[13] = df_xy(xarr[i + 1], yarr[j], ch);
    Fij[14] = df_x(xarr[i + 1], yarr[j + 1], ch);
    Fij[15] = df_xy(xarr[i + 1], yarr[j + 1], ch);
    // printf("Fij: %lf, x: %lf, y: %lf, df: %lf\n",Fij[0], xarr[i], yarr[j],
    // Fij[5]);
}

//с трансп тоже(2 р.)
static void transpose(double *matrix)
{
    int t;
    for (int i = 0; i < 4; ++i) {
        for (int j = i; j < 4; ++j) {
            t = matrix[i * 4 + j];
            matrix[i * 4 + j] = matrix[j * 4 + i];
            matrix[j * 4 + i] = t;
        }
    }
}

//(вроде тоже норм)
static void mult(double *A, double *B, double *Res)
{
    int i, j, k;

    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++) {
            Res[i * 4 + j] = 0.;
            for (k = 0; k < 4; k++) {
                Res[i * 4 + j] += A[i * 4 + k] * B[k * 4 + j];
                // if(i == 3 && j == 3)
                // printf("Elem: %lf\n",A[i * 4 + k] * B[k * 4 + j]);
            }
        }
}

//(вроде тоже)
static void calcA(int i, double *arr, double *M)
{
    double h = arr[i + 1] - arr[i];
    // printf("h: %lf \n", h);
    M[0] = 1;
    M[1] = 0;
    M[2] = 0;
    M[3] = 0;

    M[4] = 0;
    M[5] = 1;
    M[6] = 0;
    M[7] = 0;

    M[8] = -3 / (h * h * h);
    M[9] = -2 / (h * h);
    M[10] = 3 / (h * h * h);
    M[11] = -1 / (h * h);

    M[12] = 2 / (h * h * h);
    M[13] = 1 / (h * h);
    M[14] = -2 / (h * h * h);
    M[15] = 1 / (h * h);
    // printf("h: %lf\n",h);
}

void Coeff(int n, int m, double *xarr, double *yarr, double *coef, double *Fij,
           double *A, double *Res, int ch)
{
    int i, j, k;

    for (i = 0; i < n - 1; i++)
        for (j = 0; j < m - 1; j++) {
            calcFij(i, j, yarr, xarr, Fij, ch);

            calcA(i, xarr, A);

            mult(A, Fij, Res);

            calcA(j, yarr, A);

            transpose(A);

            mult(Res, A, Fij);

            for (k = 0; k < 16; k++) {
                coef[16 * (i * m + j) + k] = Fij[k];
                // printf("k: %d,i: %d, j: %d, %lf: \n",k,i ,j, A[k]);
            }
        }
}

double Value(double x, double y, int n, int m, double *coef, double *xarr,
             double *yarr, double *Temp, double *xxi, double *yyi)
{
    int i, j, k, l;
    double f_xy;
    // printf("fxy: %lf: \n",f_xy);
    for (i = 0; i < n - 2; i++)
        if (x >= xarr[i] && xarr[i + 1] >= x)
            break;
    for (j = 0; j < m - 2; j++)
        if (y >= yarr[j] && yarr[j + 1] >= y)
            break;
    // printf("x: %lf, xarr[i]: %lf, xarr[i + 1]: %lf, i: %d \n",x, xarr[i],
    // xarr[i + 1], i);

    // printf("y: %lf, yarr[j]: %lf, yarr[j + 1]: %lf, j: %d \n",y, yarr[j],
    // yarr[j + 1], j);
    for (k = 0; k < 16; k++) {
        Temp[k] = coef[16 * (i * m + j) + k];
        // printf("fxy: %lf \n",Temp[k]);
    }
    xxi[0] = 1;
    xxi[1] = (x - xarr[i]);
    xxi[2] = (x - xarr[i]) * (x - xarr[i]);
    xxi[3] = (x - xarr[i]) * (x - xarr[i]) * (x - xarr[i]);
    yyi[0] = 1;
    yyi[1] = (y - yarr[j]);
    yyi[2] = (y - yarr[j]) * (y - yarr[j]);
    yyi[3] = (y - yarr[j]) * (y - yarr[j]) * (y - yarr[j]);
    // printf("xxi[0]: %lf, xxi[1]: %lf, xxi[2]: %lf, xxi[3]: %lf\n",xxi[0],
    // xxi[1], xxi[2], xxi[3]); printf("yyi[0]: %lf, yyi[1]: %lf, yyi[2]: %lf,
    // yyi[3]: %lf\n",yyi[0], yyi[1], yyi[2], yyi[3]);
    f_xy = 0.;
    for (k = 0; k < 4; k++) {
        for (l = 0; l < 4; l++)
            f_xy += Temp[4 * k + l] * xxi[k] * yyi[l];
    }
    // printf("fxy: %lf, x: %lf,y: %lf \n",f_xy, x, y);
    return f_xy;
}

void Calc(int ch)
{
    Coeff(n, m, xarr, yarr, coef, Fij, A, Res, ch);
}

double f_xy(double x, double y)
{
    // printf("fxy: %lf \n", Value(x, y, n, m, coef, xarr, yarr, Temp, xxi,
    // yyi));
    return Value(x, y, n, m, coef, xarr, yarr, Temp, xxi, yyi);
}
