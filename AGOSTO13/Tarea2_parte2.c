#include <stdio.h>
#include <math.h>
#include <ctype.h>

// Función para calcular raíces de ecuación cuadrática
void calcularRaices(double a, double b, double c) {
    double discriminante, real, imag, x1, x2;
    
    discriminante = b*b - 4*a*c;
    
    if (discriminante > 0) {
        x1 = (-b + sqrt(discriminante)) / (2*a);
        x2 = (-b - sqrt(discriminante)) / (2*a);
        printf("Raíces reales distintas:\n");
        printf("x₁ = %.3lf\nx₂ = %.3lf\n", x1, x2);
    } 
    else if (discriminante == 0) {
        x1 = -b / (2*a);
        printf("Raíz real única (doble):\n");
        printf("x = %.3lf\n", x1);
    } 
    else {
        real = -b / (2*a);
        imag = sqrt(-discriminante) / (2*a);
        printf("Raíces complejas conjugadas:\n");
        printf("x₁ = %.3lf + %.3lfi\n", real, imag);
        printf("x₂ = %.3lf - %.3lfi\n", real, imag);
    }
}

int main() {
    double coef_a, coef_b, coef_c;
    char opcion;
    
    printf("============================================\n");
    printf("    SOLUCIONADOR DE ECUACIONES CUADRÁTICAS\n");
    printf("============================================\n");
    
    do {
        // Entrada de coeficientes
        printf("\n--- Ingrese los coeficientes ---\n");
        printf("Coeficiente a: ");
        scanf("%lf", &coef_a);
        printf("Coeficiente b: ");
        scanf("%lf", &coef_b);
        printf("Coeficiente c: ");
        scanf("%lf", &coef_c);
        
        // Validación de casos especiales
        if (coef_a == 0 && coef_b == 0) {
            printf("\n Error: Ambos coeficientes a y b son cero.\n");
            printf("La ecuación no es válida.\n");
        } 
        else if (coef_a == 0) {
            // Ecuación lineal
            double solucion = -coef_c / coef_b;
            printf("\n Ecuación lineal (a = 0):\n");
            printf("Solución: x = %.3lf\n", solucion);
        } 
        else if (coef_b == 0) {
            // Ecuación cuadrática pura
            if (-coef_c/coef_a >= 0) {
                x1 = sqrt(-coef_c/coef_a);
                x2 = -sqrt(-coef_c/coef_a);
                printf("\n Ecuación cuadrática pura (b = 0):\n");
                printf("Soluciones: x₁ = %.3lf, x₂ = %.3lf\n", x1, x2);
            } else {
                real = 0;
                imag = sqrt(coef_c/coef_a);
                printf("\n Ecuación cuadrática pura (b = 0):\n");
                printf("Soluciones complejas: x₁ = %.3lfi, x₂ = -%.3lfi\n", imag, imag);
            }
        } 
        else {
            // Ecuación cuadrática general
            printf("\n Ecuación cuadrática completa:\n");
            printf("%.2lfx² + %.2lfx + %.2lf = 0\n", coef_a, coef_b, coef_c);
            calcularRaices(coef_a, coef_b, coef_c);
        }
        
        // Preguntar si desea continuar
        printf("\n¿Desea resolver otra ecuación? (S/N): ");
        scanf(" %c", &opcion);
        
    } while (tolower(opcion) == 's');
    
    printf("\n============================================\n");
    printf("   Programa finalizado!\n");
    printf("============================================\n");
    
    return 0;
}