#include <stdio.h>
#include <math.h>

// Declaración de la función f(x)
double f(double x) {
    // Ejemplo: f(x) = x^2 - 4 (puedes cambiar esta función según necesites)
    return x * x - 4;
}

int main() {
    double a, b, k, k_old;
    double epsilon_max, epsilon_n;
    int max_iter, iter = 0;
    
    // Solicitar valores iniciales al usuario
    printf("Metodo de Biseccion (Euler)\n");
    printf("Ingrese el valor de a: ");
    scanf("%lf", &a);
    printf("Ingrese el valor de b: ");
    scanf("%lf", &b);
    printf("Ingrese el error maximo (epsilon_max): ");
    scanf("%lf", &epsilon_max);
    printf("Ingrese el numero maximo de iteraciones: ");
    scanf("%d", &max_iter);
    
    // Verificar que f(a) y f(b) tienen signos opuestos
    if (f(a) * f(b) >= 0) {
        printf("Error: f(a) y f(b) deben tener signos opuestos.\n");
        printf("f(a) = %lf, f(b) = %lf\n", f(a), f(b));
        return 1;
    }
    
    k_old = a; // Inicializar k_old
    
    printf("\nIter\t  a\t\t  b\t\t  k\t\t  f(k)\t\t  Error\n");
    printf("------------------------------------------------------------------------\n");
    
    do {
        k_old = k; // Guardar el valor anterior de k
        k = (a + b) / 2.0; // Calcular el punto medio
        
        printf("%3d\t%10.6lf\t%10.6lf\t%10.6lf\t%10.6lf", iter, a, b, k, f(k));
        
        if (iter > 0) {
            epsilon_n = fabs((k - k_old) / k);
            printf("\t%10.6lf", epsilon_n);
        } else {
            printf("\t      -");
        }
        printf("\n");
        
        // Verificar si hemos encontrado la raíz exacta
        if (f(k) == 0.0) {
            break;
        }
        // Actualizar los intervalos
        else if (f(k) * f(a) < 0) {
            b = k;
        } else {
            a = k;
        }
        
        iter++;
        
    } while (iter < max_iter && (iter == 0 || fabs((k - k_old) / k) > epsilon_max));
    
    printf("\nRaiz aproximada: %.8lf\n", k);
    printf("f(raiz) = %.8lf\n", f(k));
    printf("Iteraciones realizadas: %d\n", iter);
    
    if (iter < max_iter) {
        printf("Convergencia alcanzada con error: %.8lf\n", epsilon_n);
    } else {
        printf("Numero maximo de iteraciones alcanzado.\n");
    }
    
    return 0;
}