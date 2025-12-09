#include <stdio.h>
#include <math.h>

double k = 0.45;
double Cmax = 1.0;

double f(double t, double C) {
    return k * C * (Cmax - C);
}

int main(void) {
    double dt   = 0.1;     // paso de tiempo
    double tmax = 25.0;    // tiempo final
    double t    = 0.0;

    // DOS variables: una para RK4 y otra para RK2
    double C4 = 0.02;      // C(0) para RK4
    double C2 = 0.02;      // C(0) para RK2

    printf("# t(min)\tC_RK4(mol/L)\tC_RK2(mol/L)\n");

    while (t <= tmax + 1e-9) {
        printf("%8.3f\t%10.6f\t%10.6f\n", t, C4, C2);

        /* === Paso RK4 para C4 === */
        double k1 = f(t, C4);
        double k2 = f(t + 0.5*dt, C4 + 0.5*dt*k1);
        double k3 = f(t + 0.5*dt, C4 + 0.5*dt*k2);
        double k4 = f(t + dt,     C4 + dt*k3);
        C4 = C4 + dt * (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0;

        /* === Paso RK2 (segundo orden) para C2 === */
        double k1_2 = f(t, C2);
        double k2_2 = f(t + dt, C2 + dt*k1_2);
        C2 = C2 + dt * (k1_2 + k2_2) / 2.0;

        t = t + dt;
    }

    return 0;
}

