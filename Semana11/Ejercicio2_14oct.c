/* gauss_2x2.c — Resolver los dos sistemas del pizarrón con Gauss (pivoteo parcial) */

#include <stdio.h>
#include <math.h>

#define N 2
#define EPS 1e-12

static void swap_rows(double A[N][N], double b[N], int r1, int r2) {
    if (r1 == r2) return;
    for (int j = 0; j < N; ++j) {
        double t = A[r1][j];
        A[r1][j] = A[r2][j];
        A[r2][j] = t;
    }
    double tb = b[r1]; b[r1] = b[r2]; b[r2] = tb;
}

int gauss_solve(double A[N][N], double b[N], double x[N]) {
    for (int k = 0; k < N; ++k) {
        int piv = k; double maxabs = fabs(A[k][k]);
        for (int i = k+1; i < N; ++i) {
            double v = fabs(A[i][k]);
            if (v > maxabs) { maxabs = v; piv = i; }
        }
        if (maxabs < EPS) return -1;
        swap_rows(A, b, k, piv);
        for (int i = k+1; i < N; ++i) {
            double m = A[i][k] / A[k][k];
            A[i][k] = 0.0;
            for (int j = k+1; j < N; ++j) A[i][j] -= m * A[k][j];
            b[i] -= m * b[k];
        }
    }
    for (int i = N-1; i >= 0; --i) {
        double s = b[i];
        for (int j = i+1; j < N; ++j) s -= A[i][j] * x[j];
        if (fabs(A[i][i]) < EPS) return -1;
        x[i] = s / A[i][i];
    }
    return 0;
}

static void solve_and_print(double A[N][N], double b[N]) {
    double x[N];
    double Ac[N][N], bc[N];
    for (int i=0;i<N;++i){ for(int j=0;j<N;++j) Ac[i][j]=A[i][j]; bc[i]=b[i]; }
    if (gauss_solve(Ac, bc, x) == 0)
        printf("x1 = %.10f\nx2 = %.10f\n\n", x[0], x[1]);
    else
        printf("Sistema singular\n\n");
}

int main(void) {
    // (1) x1 + 2x2 = 10 ; 1.1 x1 + 2x2 = 10.4   -> x1=4, x2=3
    double A1[N][N] = {{1.0, 2.0},
                       {1.1, 2.0}};
    double b1[N] = {10.0, 10.4};

    // (2) x1 + 2x2 = 10 ; 1.05 x1 + 2x2 = 10.4  -> x1=8, x2=1
    double A2[N][N] = {{1.0, 2.0},
                       {1.05, 2.0}};
    double b2[N] = {10.0, 10.4};

    solve_and_print(A1, b1);
    solve_and_print(A2, b2);
    return 0;
}
