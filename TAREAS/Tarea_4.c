/*Tarea 4 10 DE SEPTIEMBRE DEL 2025*/

#include <stdio.h>

int main() {
    double x0, x1, x2;
    double f0, f1, f2;
    double derivada_real;

    // Entrada de datos por teclado
    printf("=== Metodo de Diferencias Finitas ===\n");
    printf("Ingresa x0: "); scanf("%lf", &x0);
    printf("Ingresa f(x0): "); scanf("%lf", &f0);

    printf("Ingresa x1: "); scanf("%lf", &x1);
    printf("Ingresa f(x1): "); scanf("%lf", &f1);

    printf("Ingresa x2: "); scanf("%lf", &x2);
    printf("Ingresa f(x2): "); scanf("%lf", &f2);

    printf("Ingresa la derivada real en x1: "); scanf("%lf", &derivada_real);

    // Cálculos de derivadas aproximadas
    double adelante = (f2 - f1) / (x2 - x1);  // Hacia adelante
    double atras = (f1 - f0) / (x1 - x0);     // Hacia atrás
    double centrada = (f2 - f0) / (x2 - x0);  // Centrada

    // Cálculos de errores relativos
    double error_adelante = (derivada_real - adelante) / derivada_real;
    double error_atras = (derivada_real - atras) / derivada_real;
    double error_centrada = (derivada_real - centrada) / derivada_real;

    // Resultados
    printf("\n=== Resultados en x = %.2f ===\n", x1);
    printf("Hacia adelante: %.6f\n", adelante);
    printf("Hacia atras: %.6f\n", atras);
    printf("Centrada: %.6f\n", centrada);

    printf("\n=== Errores relativos ===\n");
    printf("Hacia adelante: %.6f\n", error_adelante);
    printf("Hacia atras: %.6f\n", error_atras);
    printf("Centrada: %.6f\n", error_centrada);

    return 0;
}

