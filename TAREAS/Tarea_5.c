#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define MAX_ITER 100
#define ERROR_RELATIVO 0.001

// Función f(x) = sin(10x) - cos(3x)
double f(double x) {
    return sin(10 * x) - cos(3 * x);
}

// Método de bisección para encontrar una raíz en el intervalo [a, b]
int bisection(double a, double b, double *raiz, double *error_relativo, int *iteraciones) {
    double fa, fb, fc, c, c_anterior;
    int iter = 0;
    
    fa = f(a);
    fb = f(b);
    
    // Verificar si hay cambio de signo en el intervalo
    if (fa * fb > 0) {
        return 0; // No hay raíz en este intervalo
    }
    
    c_anterior = a;
    
    for (iter = 0; iter < MAX_ITER; iter++) {
        c = (a + b) / 2.0;
        fc = f(c);
        
        // Calcular error relativo aproximado
        if (iter > 0) {
            *error_relativo = fabs((c - c_anterior) / c) * 100;
        }
        
        // Verificar si encontramos la raíz exacta
        if (fc == 0.0) {
            *raiz = c;
            *iteraciones = iter + 1;
            return 1;
        }
        
        // Determinar en qué subintervalo está la raíz
        if (fa * fc < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
        
        c_anterior = c;
        
        // Verificar criterio de convergencia
        if (iter > 0 && *error_relativo < ERROR_RELATIVO) {
            *raiz = c;
            *iteraciones = iter + 1;
            return 1;
        }
    }
    
    *raiz = c;
    *iteraciones = MAX_ITER;
    return 1;
}

// Función para encontrar todas las raíces en el intervalo [x_min, x_max]
void encontrar_raices(double x_min, double x_max, double paso) {
    double a, b;
    double raiz, error;
    int iteraciones, encontrada;
    int contador_raices = 0;
    
    printf("Buscando raices de f(x) = sin(10x) - cos(3x) en [%.2f, %.2f]\n", x_min, x_max);
    printf("Error relativo maximo: %.4f\n\n", ERROR_RELATIVO);
    
    // Dividir el intervalo en subintervalos más pequeños
    a = x_min;
    while (a < x_max) {
        b = a + paso;
        if (b > x_max) b = x_max;
        
        // Verificar si hay cambio de signo en el subintervalo
        if (f(a) * f(b) <= 0) {
            encontrada = bisection(a, b, &raiz, &error, &iteraciones);
            
            if (encontrada) {
                contador_raices++;
                printf("Raiz %d encontrada:\n", contador_raices);
                printf("  x = %.6f\n", raiz);
                printf("  f(x) = %.6e\n", f(raiz));
                printf("  Iteraciones: %d\n", iteraciones);
                printf("  Error relativo final: %.6f%%\n", error);
                
                // Verificar la relación Ea = δx / 2^n
                double delta_x = b - a;
                double ea_teorico = delta_x / pow(2, iteraciones);
                double ea_real = fabs(raiz - ((a + b) / 2.0));
                
                printf("  Ea teorico (δx/2^n): %.6e\n", ea_teorico);
                printf("  Ea real: %.6e\n", ea_real);
                printf("  Relacion satisfecha: %s\n\n", 
                       fabs(ea_real - ea_teorico) < 1e-10 ? "SI" : "NO");
            }
        }
        a = b;
    }
    
    if (contador_raices == 0) {
        printf("No se encontraron raices en el intervalo especificado.\n");
    } else {
        printf("Total de raices encontradas: %d\n", contador_raices);
    }
}

int main() {
    double x_min = 3.0;
    double x_max = 5.0;
    double paso = 0.1; // Tamaño del paso para buscar subintervalos
    
    printf("METODO DE BISECCION\n");
    printf("===================\n");
    
    encontrar_raices(x_min, x_max, paso);
    
    return 0;
}