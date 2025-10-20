#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define EPS 1e-12

int leer_sistema(const char *ruta, int *n, double ***A, double **b) {
    FILE *f = fopen(ruta, "r");
    if (!f) { perror("No pude abrir el archivo"); return -1; }

    if (fscanf(f, "%d", n) != 1 || *n <= 0) {
        fprintf(stderr, "No pude leer n.\n");
        fclose(f); return -1;
    }
    int N = *n;

    double **AA = (double**) malloc(N * sizeof *AA);
    double *bb  = (double*)  malloc(N * sizeof *bb);
    if (!AA || !bb) { fprintf(stderr, "Memoria insuficiente.\n"); fclose(f); return -1; }

    for (int i = 0; i < N; ++i) {
        AA[i] = (double*) malloc(N * sizeof *AA[i]);
        if (!AA[i]) { fprintf(stderr, "Memoria insuficiente.\n"); fclose(f); return -1; }
    }

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            if (fscanf(f, "%lf", &AA[i][j]) != 1) { fprintf(stderr, "No pude leer A[%d][%d].\n", i, j); fclose(f); return -1; }

    for (int i = 0; i < N; ++i)
        if (fscanf(f, "%lf", &bb[i]) != 1) { fprintf(stderr, "No pude leer b[%d].\n", i); fclose(f); return -1; }

    fclose(f);
    *A = AA; *b = bb;
    return 0;
}

static void swap_rows(double **A, double *b, int r1, int r2, int N) {
    if (r1 == r2) return;
    for (int j = 0; j < N; ++j) {
        double t = A[r1][j]; A[r1][j] = A[r2][j]; A[r2][j] = t;
    }
    double tb = b[r1]; b[r1] = b[r2]; b[r2] = tb;
}

static double norm_inf(double **A, int N) {
    double n = 0.0;
    for (int i = 0; i < N; ++i) {
        double s = 0.0;
        for (int j = 0; j < N; ++j) s += fabs(A[i][j]);
        if (s > n) n = s;
    }
    return n;
}

static void copia_sistema(double **A, double *b, double ***Ao, double **bo, int N) {
    double **Ac = (double**) malloc(N * sizeof *Ac);
    double *bc = (double*) malloc(N * sizeof *bc);
    for (int i = 0; i < N; ++i) {
        Ac[i] = (double*) malloc(N * sizeof *Ac[i]);
        for (int j = 0; j < N; ++j) Ac[i][j] = A[i][j];
        bc[i] = b[i];
    }
    *Ao = Ac; *bo = bc;
}

static void libera_sistema(double **A, double *b, int N) {
    for (int i = 0; i < N; ++i) free(A[i]);
    free(A); free(b);
}



int gauss_eliminacion(double **A, double *b, double *x, int N, double *rcond_est_out) {
    double A_norm = norm_inf(A, N);
    double min_pivot = INFINITY;

    for (int k = 0; k < N; ++k) {
        // pivoteo parcial
        int piv = k; double maxabs = fabs(A[k][k]);
        for (int i = k+1; i < N; ++i) {
            double v = fabs(A[i][k]);
            if (v > maxabs) { maxabs = v; piv = i; }
        }
        if (piv != k) {
            printf("Gauss: intercambio de filas R%d <-> R%d\n", k+1, piv+1);
            swap_rows(A, b, k, piv, N);
        }
        if (fabs(A[k][k]) < EPS) { fprintf(stderr, "Gauss: pivote ~0 en k=%d\n", k); return -1; }

        // normalización fila pivote
        double pivot = A[k][k];
        min_pivot = fmin(min_pivot, fabs(pivot));
        for (int j = k; j < N; ++j) A[k][j] /= pivot;
        b[k] /= pivot;
        printf("Gauss: normalizo R%d (pivot=%.4e)\n", k+1, pivot);

        // eliminación por debajo
        for (int i = k+1; i < N; ++i) {
            double m = A[i][k]; // pivote ya es 1
            if (fabs(m) > 0) {
                for (int j = k; j < N; ++j) A[i][j] -= m * A[k][j];
                b[i] -= m * b[k];
            }
        }
    }

    // sustitución hacia atrás (diagonal = 1)
    for (int i = N-1; i >= 0; --i) {
        double s = b[i];
        for (int j = i+1; j < N; ++j) s -= A[i][j] * x[j];
        x[i] = s;
    }

    double rcond_est = min_pivot / (A_norm + 1e-16);
    if (rcond_est_out) *rcond_est_out = rcond_est;
    return 0;
}

/* -------------------- Gauss-Jordan -------------------- */
int gauss_jordan(double **A, double *b, double *x, int N, double *rcond_est_out) {
    double A_norm = norm_inf(A, N);
    double min_pivot = INFINITY;

    // pivoteo + normalización
    for (int k = 0; k < N; ++k) {
        int piv = k; double maxabs = fabs(A[k][k]);
        for (int i = k+1; i < N; ++i) {
            double v = fabs(A[i][k]);
            if (v > maxabs) { maxabs = v; piv = i; }
        }
        if (piv != k) {
            printf("G-J: intercambio de filas R%d <-> R%d\n", k+1, piv+1);
            swap_rows(A, b, k, piv, N);
        }
        if (fabs(A[k][k]) < EPS) { fprintf(stderr, "G-J: pivote ~0 en k=%d\n", k); return -1; }

        double pivot = A[k][k];
        min_pivot = fmin(min_pivot, fabs(pivot));
        for (int j = 0; j < N; ++j) A[k][j] /= pivot;
        b[k] /= pivot;
        printf("G-J: normalizo R%d (pivot=%.4e)\n", k+1, pivot);

        
        for (int i = 0; i < N; ++i) if (i != k) {
            double m = A[i][k];
            if (fabs(m) > 0) {
                for (int j = 0; j < N; ++j) A[i][j] -= m * A[k][j];
                b[i] -= m * b[k];
            }
        }
    }

    
    for (int i = 0; i < N; ++i) x[i] = b[i];

    double rcond_est = min_pivot / (A_norm + 1e-16);
    if (rcond_est_out) *rcond_est_out = rcond_est;
    return 0;
}


static void uso(const char *p) {
    fprintf(stderr,
        "Uso: %s <archivo> [--gauss | --jordan | --ambos]\n"
        "   --gauss   : Eliminacion de Gauss (por defecto)\n"
        "   --jordan  : Gauss-Jordan\n"
        "   --ambos   : Ejecuta ambos metodos\n", p);
}

int main(int argc, char **argv) {
    if (argc < 2) { uso(argv[0]); return 1; }

    int modo = 0; // 0=gauss, 1=jordan, 2=ambos
    if (argc >= 3) {
        if      (!strcmp(argv[2], "--gauss"))   modo = 0;
        else if (!strcmp(argv[2], "--jordan"))  modo = 1;
        else if (!strcmp(argv[2], "--ambos"))   modo = 2;
        else { uso(argv[0]); return 1; }
    }

    int n; double **A0, *b0;
    if (leer_sistema(argv[1], &n, &A0, &b0) != 0) return 1;

    // A y b
    double **A1, *b1, **A2, *b2;
    copia_sistema(A0, b0, &A1, &b1, n);
    copia_sistema(A0, b0, &A2, &b2, n);

    double *x = (double*) calloc(n, sizeof *x);
    if (!x) { fprintf(stderr, "Memoria insuficiente.\n"); return 1; }

    int run_gauss  = (modo == 0 || modo == 2);
    int run_jordan = (modo == 1 || modo == 2);

    if (run_gauss) {
        printf("=== Metodo: Gauss (pivoteo parcial + normalizacion) ===\n");
        double rce;
        if (gauss_eliminacion(A1, b1, x, n, &rce) == 0) {
            if (rce < 1e-8) printf("[AVISO] Sistema mal condicionado (rcond≈%.2e)\n", rce);
            else            printf("Sistema bien condicionado (rcond≈%.2e)\n", rce);
            printf("Solucion (Gauss):\n");
            for (int i = 0; i < n; ++i) printf("x[%d] = %.10f\n", i+1, x[i]);
        } else {
            printf("Gauss: sistema singular o numericamente inestable.\n");
        }
        printf("\n");
    }

    if (run_jordan) {
        
        for (int i = 0; i < n; ++i) x[i] = 0.0;
        printf("=== Metodo: Gauss-Jordan (pivoteo parcial + normalizacion) ===\n");
        double rce;
        if (gauss_jordan(A2, b2, x, n, &rce) == 0) {
            if (rce < 1e-8) printf("[AVISO] Sistema mal condicionado (rcond≈%.2e)\n", rce);
            else            printf("Sistema bien condicionado (rcond≈%.2e)\n", rce);
            printf("Solucion (Gauss-Jordan):\n");
            for (int i = 0; i < n; ++i) printf("x[%d] = %.10f\n", i+1, x[i]);
        } else {
            printf("Gauss-Jordan: sistema singular o numericamente inestable.\n");
        }
        printf("\n");
    }

    
    libera_sistema(A0, b0, n);
    libera_sistema(A1, b1, n);
    libera_sistema(A2, b2, n);
    free(x);
    return 0;
}
