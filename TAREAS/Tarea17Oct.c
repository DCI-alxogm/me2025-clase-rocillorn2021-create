#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPS 1e-12

/* Lee n, A (n x n) y b (n) desde archivo.
   Devuelve 0 si ok, -1 si error. */
int leer_sistema(const char *ruta, int *n, double ***A, double **b) {
    FILE *f = fopen(ruta, "r");
    if (!f) { perror("No pude abrir el archivo"); return -1; }

    if (fscanf(f, "%d", n) != 1 || *n <= 0) {
        fprintf(stderr, "No pude leer n.\n");
        fclose(f);
        return -1;
    }

    int N = *n;
    double **AA = (double**) malloc(N * sizeof *AA);
    double *bb = (double*) malloc(N * sizeof *bb);
    if (!AA || !bb) { fprintf(stderr, "Memoria insuficiente.\n"); fclose(f); return -1; }

    for (int i = 0; i < N; ++i) {
        AA[i] = (double*) malloc(N * sizeof *AA[i]);
        if (!AA[i]) { fprintf(stderr, "Memoria insuficiente.\n"); fclose(f); return -1; }
    }

    // Leer A
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (fscanf(f, "%lf", &AA[i][j]) != 1) {
                fprintf(stderr, "No pude leer A[%d][%d].\n", i, j);
                fclose(f); return -1;
            }
        }
    }

    // Leer b
    for (int i = 0; i < N; ++i) {
        if (fscanf(f, "%lf", &bb[i]) != 1) {
            fprintf(stderr, "No pude leer b[%d].\n", i);
            fclose(f); return -1;
        }
    }
    fclose(f);
    *A = AA; *b = bb;
    return 0;
}

/* Intercambia filas i y k (A y b) */
static void swap_rows(double **A, double *b, int i, int k, int N) {
    if (i == k) return;
    for (int j = 0; j < N; ++j) {
        double t = A[i][j]; A[i][j] = A[k][j]; A[k][j] = t;
    }
    double tb = b[i]; b[i] = b[k]; b[k] = tb;
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


int gauss_resolver(double **A, double *b, double *x, int N) {
    
    double A_norm = norm_inf(A, N);
    double min_pivot = INFINITY;
    int poor_scaling = 0;

    for (int k = 0; k < N; ++k) {
        
        int piv = k;
        double maxabs = fabs(A[k][k]);
        for (int i = k + 1; i < N; ++i) {
            double v = fabs(A[i][k]);
            if (v > maxabs) { maxabs = v; piv = i; }
        }

        
        if (piv != k) {
            printf("Intercambio de filas R%d <-> R%d\n", k+1, piv+1);
            swap_rows(A, b, k, piv, N);
        }

        // Chequeo de pivote ~0
        if (fabs(A[k][k]) < EPS) {
            fprintf(stderr, "Pivote ~0 en k=%d. Sistema singular o mal condicionado.\n", k);
            return -1;
        }

        // Normalizar fila pivote (pivot -> 1)
        double pivot = A[k][k];
        min_pivot = fmin(min_pivot, fabs(pivot));
        for (int j = k; j < N; ++j) A[k][j] /= pivot;
        b[k] /= pivot;
        printf("Normalizo fila R%d (pivot = %.4e)\n", k+1, pivot);

        // Eliminar por debajo
        for (int i = k + 1; i < N; ++i) {
            double m = A[i][k];          // ya que el pivote es 1, m = A[i][k]
            if (fabs(m) > 0.0) {
                for (int j = k; j < N; ++j) A[i][j] -= m * A[k][j];
                b[i] -= m * b[k];
            }
        }
    }

    // Sustitución hacia atrás (ya triangular superior con 1s en diagonal)
    for (int i = N - 1; i >= 0; --i) {
        double s = b[i];
        for (int j = i + 1; j < N; ++j) s -= A[i][j] * x[j];
        x[i] = s; // diagonal es 1
    }

    
    double rcond_est = min_pivot / (A_norm + 1e-16);
    if (rcond_est < 1e-8) {
        printf("\n[AVISO] El sistema parece MAL CONDICIONADO (rcond≈%.2e).\n", rcond_est);
    } else {
        printf("\nEl sistema parece BIEN CONDICIONADO (rcond≈%.2e).\n", rcond_est);
    }

    (void)poor_scaling;
return 0;

}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Uso: %s <archivo_entrada>\n", argv[0]);
        fprintf(stderr, "Formato:\n n\n A(nxn)\n b(n)\n");
        return 1;
    }

    int n;
    double **A, *b;
    if (leer_sistema(argv[1], &n, &A, &b) != 0) return 1;

    double *x = (double*) calloc(n, sizeof *x);
    if (!x) { fprintf(stderr, "Memoria insuficiente.\n"); return 1; }

    printf("=== Resolviendo A x = b (n=%d) ===\n", n);
    if (gauss_resolver(A, b, x, n) != 0) {
        fprintf(stderr, "No se pudo resolver el sistema.\n");
        return 1;
    }

    
    printf("\nSolucion:\n");
    for (int i = 0; i < n; ++i)
        printf("x[%d] = %.10f\n", i+1, x[i]);

    
    for (int i = 0; i < n; ++i) free(A[i]);
    free(A); free(b); free(x);
    return 0;
}
