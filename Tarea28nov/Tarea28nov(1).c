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

    printf("Metodo de Euler para y' = -2y, y(0)=1\n");
    printf("t\t\t y_num\n");
    printf("%g\t %g\n", t, y);

    for (int i = 0; i < N; ++i) {
       
        y = y + h * f(t, y);
        t = t + h;

        printf("%g\t %g\n", t, y);
    }

    return 0;
}