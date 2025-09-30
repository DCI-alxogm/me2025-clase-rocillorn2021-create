/* raices_metodos.c
   Un solo programa que calcula la MISMA raíz de f(x)=exp(-x)-x
   por distintos métodos numéricos, reportando el error aproximado
   relativo en cada iteración hasta cumplir una tolerancia.

   Métodos: Bisección, Falsa Posición, Punto Fijo, Newton-Raphson, Secante.
*/
#include <stdio.h>
#include <math.h>
#include <float.h>

#ifndef M_E
#define M_E 2.71828182845904523536
#endif

/* ============================
   Definición de la función
   ============================ */
double f(double x) {
    /* Función NO polinómica */
    return exp(-x) - x;
}
double df(double x) {
    /* Derivada para Newton */
    return -exp(-x) - 1.0;
}

/* ============================
   Utilidades
   ============================ */
double rel_error(double nuevo, double viejo, int have_old) {
    if (!have_old || nuevo == 0.0) return -1.0; /* -1 indica n/a */
    return fabs((nuevo - viejo) / nuevo);
}

void print_header(const char *h1, const char *h2, const char *h3,
                  const char *h4, const char *h5, const char *h6) {
    printf("%-4s | %-14s | %-14s | %-14s | %-14s | %-12s\n",
           h1, h2, h3, h4, h5, h6);
    printf("-----+----------------+----------------+----------------+----------------+--------------\n");
}

/* ============================
   Métodos
   ============================ */
double biseccion(double a, double b, double tol, int maxit) {
    double fa = f(a), fb = f(b);
    if (fa * fb > 0.0) {
        printf("Bisección: el intervalo inicial no encierra una raíz.\n");
        return NAN;
    }
    printf("\n>> Bisección  [a=%.6f, b=%.6f]\n", a, b);
    print_header("it","a","b","c","f(c)","e_a (rel)");

    double prev = 0.0; int have_prev = 0;
    for (int k=1; k<=maxit; ++k) {
        double c = 0.5*(a+b);
        double fc = f(c);
        double ea = rel_error(c, prev, have_prev);

        printf("%-4d | %-14.8f | %-14.8f | %-14.8f | %-14.8e | %-12s\n",
               k, a, b, c, fc, (ea<0? "-" : (ea<1e-99? "0" : ""))), printf("");

        if (ea >= 0 && ea < tol) { return c; }
        if (fa * fc < 0.0) { b = c; fb = fc; } else { a = c; fa = fc; }
        prev = c; have_prev = 1;
    }
    return 0.5*(a+b);
}

double falsa_posicion(double a, double b, double tol, int maxit) {
    double fa = f(a), fb = f(b);
    if (fa * fb > 0.0) {
        printf("Falsa Posición: el intervalo inicial no encierra una raíz.\n");
        return NAN;
    }
    printf("\n>> Falsa Posición  [a=%.6f, b=%.6f]\n", a, b);
    print_header("it","a","b","x_r","f(x_r)","e_a (rel)");

    double prev = 0.0; int have_prev = 0;
    for (int k=1; k<=maxit; ++k) {
        double xr = b - fb*(b-a)/(fb-fa);
        double fr = f(xr);
        double ea = rel_error(xr, prev, have_prev);

        printf("%-4d | %-14.8f | %-14.8f | %-14.8f | %-14.8e | %-12s\n",
               k, a, b, xr, fr, (ea<0? "-" : (ea<1e-99? "0" : ""))), printf("");

        if (ea >= 0 && ea < tol) { return xr; }
        if (fa * fr < 0.0) { b = xr; fb = fr; } else { a = xr; fa = fr; }
        prev = xr; have_prev = 1;
    }
    return b - fb*(b-a)/(fb-fa);
}

double punto_fijo(double x0, double tol, int maxit) {
    /* g(x) = e^{-x}  =>  x = g(x) */
    printf("\n>> Punto Fijo  [x0=%.6f,  x_{n+1}=e^{-x_n}]\n", x0);
    print_header("it","x_n","x_{n+1}","f(x_{n+1})","-","e_a (rel)");

    double prev = 0.0; int have_prev = 0;
    double x = x0;
    for (int k=1; k<=maxit; ++k) {
        double x1 = exp(-x);
        double ea = rel_error(x1, prev, have_prev);

        printf("%-4d | %-14.8f | %-14.8f | %-14.8e | %-14s | %-12s\n",
               k, x, x1, f(x1), "", (ea<0? "-" : (ea<1e-99? "0" : ""))), printf("");

        if (ea >= 0 && ea < tol) { return x1; }
        prev = x1; have_prev = 1; x = x1;
    }
    return x;
}

double newton(double x0, double tol, int maxit) {
    printf("\n>> Newton-Raphson  [x0=%.6f]\n", x0);
    print_header("it","x_n","x_{n+1}","f(x_{n+1})","-","e_a (rel)");

    double prev = 0.0; int have_prev = 0;
    double x = x0;
    for (int k=1; k<=maxit; ++k) {
        double dfx = df(x);
        double x1  = x - f(x)/dfx;
        double ea  = rel_error(x1, prev, have_prev);

        printf("%-4d | %-14.8f | %-14.8f | %-14.8e | %-14s | %-12s\n",
               k, x, x1, f(x1), "", (ea<0? "-" : (ea<1e-99? "0" : ""))), printf("");

        if (ea >= 0 && ea < tol) { return x1; }
        prev = x1; have_prev = 1; x = x1;
    }
    return x;
}

double secante(double x0, double x1, double tol, int maxit) {
    printf("\n>> Secante  [x0=%.6f, x1=%.6f]\n", x0, x1);
    print_header("it","x_{n-1}","x_n","x_{n+1}","f(x_{n+1})","e_a (rel)");

    double prev = 0.0; int have_prev = 0;
    for (int k=1; k<=maxit; ++k) {
        double f0 = f(x0), f1 = f(x1);
        double x2 = x1 - f1*(x1-x0)/(f1-f0);
        double ea = rel_error(x2, prev, have_prev);

        printf("%-4d | %-14.8f | %-14.8f | %-14.8f | %-14.8e | %-12s\n",
               k, x0, x1, x2, f(x2), (ea<0? "-" : (ea<1e-99? "0" : ""))), printf("");

        if (ea >= 0 && ea < tol) { return x2; }
        x0 = x1; x1 = x2; prev = x2; have_prev = 1;
    }
    return x1;
}

/* ============================
   Main
   ============================ */
int main(void) {
    const double tol   = 1e-8;
    const int    maxit = 200;

    printf("\n==== Raíz de f(x) = exp(-x) - x (MISMA raíz con todos los métodos) ====\n");

    double xb = biseccion(0.0, 1.0, tol, maxit);
    printf("Raíz (Bisección)      ≈ %.12f   | f(x)≈ %.3e\n", xb, f(xb));

    double xr = falsa_posicion(0.0, 1.0, tol, maxit);
    printf("Raíz (Falsa Posición) ≈ %.12f   | f(x)≈ %.3e\n", xr, f(xr));

    double xg = punto_fijo(0.0, tol, 1000);
    printf("Raíz (Punto Fijo)     ≈ %.12f   | f(x)≈ %.3e\n", xg, f(xg));

    double xn = newton(0.5, tol, 50);
    printf("Raíz (Newton)         ≈ %.12f   | f(x)≈ %.3e\n", xn, f(xn));

    double xs = secante(0.0, 1.0, tol, 50);
    printf("Raíz (Secante)        ≈ %.12f   | f(x)≈ %.3e\n", xs, f(xs));

    printf("\n==== Resumen ====\n");
    printf("Bisección        : %.12f\n", xb);
    printf("Falsa Posición   : %.12f\n", xr);
    printf("Punto Fijo       : %.12f\n", xg);
    printf("Newton-Raphson   : %.12f\n", xn);
    printf("Secante          : %.12f\n", xs);
    return 0;
}
