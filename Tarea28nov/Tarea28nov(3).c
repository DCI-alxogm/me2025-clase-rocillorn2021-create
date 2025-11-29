#include <stdio.h>
#include <math.h>


double f(double t, double y) {
    (void)t;
    return -2.0 * y;
}

int main(void) {
    double t0 = 0.0;
    double y0 = 1.0;
    double tf = 2.0;
    int    N  = 20;
    double h  = (tf - t0) / N;

    double t = t0;
    double y = y0;

    printf("Metodo RK4 para y' = -2y, y(0)=1\n");
    printf("t\t\t y_num\n");
    printf("%g\t %g\n", t, y);

    for (int i = 0; i < N; ++i) {
        double k1 = f(t, y);
        double k2 = f(t + 0.5*h, y + 0.5*h*k1);
        double k3 = f(t + 0.5*h, y + 0.5*h*k2);
        double k4 = f(t + h,     y + h*k3);

        
        y = y + (h/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
        t = t + h;

        printf("%g\t %g\n", t, y);
    }

    return 0;
}