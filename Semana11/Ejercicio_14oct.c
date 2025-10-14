/* gauss.c — Eliminación de Gauss con pivoteo parcial */

#include <stdio.h>
#include <math.h>

#define N 3
#define EPS 1e-12

static void swap_rows(double A[N][N], double b[N], int r1, int r2) {
    if (r1 == r2) return;
    for (int j = 0; j < N; ++j) {
        double tmp = A[r1][j];
        A[r1][j] = A[r2][j];
        A[r2][j] = tmp;
    }
    double tb = b[r1];
    b[r1] = b[r2];
    b[r2] = tb;
}

int gauss_solve(double A[N][N], double b[N], double x[N]) {
    for (int k = 0; k < N; ++k) {
        int piv = k;
        double maxabs = fabs(A[k][k]);
        for (int i = k + 1; i < N; ++i) {
            double v = fabs(A[i][k]);
            if (v > maxabs) { maxabs = v; piv = i; }
        }
        if (maxabs < EPS) return -1;
        swap_rows(A, b, k, piv);

        for (int i = k + 1; i < N; ++i) {
            double m = A[i][k] / A[k][k];
            A[i][k] = 0.0;
            for (int j = k + 1; j < N; ++j) A[i][j] -= m * A[k][j];
            b[i] -= m * b[k];
        }
    }
    for (int i = N - 1; i >= 0; --i) {
        double s = b[i];
        for (int j = i + 1; j < N; ++j) s -= A[i][j] * x[j];
        if (fabs(A[i][i]) < EPS) return -1;
        x[i] = s / A[i][i];
    }
    return 0;
}

int main(void) {
    double A[N][N] = {{4,6,7},{0,2,3},{2,1,6}};
    double b[N] = {-3,8,5};
    double x[N];

    if (gauss_solve(A, b, x) != 0) {
        fprintf(stderr, "Sistema singular.\n");
        return 1;
    }
    for (int i = 0; i < N; ++i) printf("x%d = %.10f\n", i+1, x[i]);
    return 0;
}
