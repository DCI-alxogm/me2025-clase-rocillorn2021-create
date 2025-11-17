#include <stdio.h>
#include <math.h>

/* Ajustes generales */
#define DIM 5 /* es una constante que reemplaza los 5 nodos, para no poner uno por uno cada variable */ 
#define LIMITE_ERR 1e-6
#define CICLOS_MAX 1000

/* Función por nodo) */
void reporteSolucion(double vector[DIM], const char *etiqueta) {
    int i;
    printf("\n%s:\n", etiqueta);
    for(i = 0; i < DIM; i++){
        printf("n%d = %lf\n", i+1, vector[i]);
    }
}

/* Matriz */
void copiarMatriz(double origen[DIM][DIM], double destino[DIM][DIM]) {
    int i, j;
    for(i = 0; i < DIM; i++){
        for(j = 0; j < DIM; j++){
            destino[i][j] = origen[i][j];
        }
    }
}

/*  Método Gauss–Jordan  */
void metodo_GaussJordan(double matriz[DIM][DIM], double vector_b[DIM], double solucion[DIM]) {
    double expandida[DIM][DIM+1];
    int i, j, p;

    /* Matriz aumentada */
    for(i = 0; i < DIM; i++){
        for(j = 0; j < DIM; j++){
            expandida[i][j] = matriz[i][j];
        }
        expandida[i][DIM] = vector_b[i];
    }

    /* Eliminación por filas */
    for(p = 0; p < DIM; p++){

        /* Normalización por pivote */
        double piv = expandida[p][p];
        for(j = 0; j <= DIM; j++){
            expandida[p][j] /= piv;
        }

        /* Anulación de pivoteo */
        for(i = 0; i < DIM; i++){
            if(i != p){
                double factor = expandida[i][p];
                for(j = 0; j <= DIM; j++){
                    expandida[i][j] -= factor * expandida[p][j];
                }
            }
        }
    }


    for(i = 0; i < DIM; i++){
        solucion[i] = expandida[i][DIM];
    }
}

/*  Método Gauss–Seidel */
void metodo_GaussSeidel(double matriz[DIM][DIM], double vector_b[DIM], double solucion[DIM]) {
    double previo[DIM] = {0.0};
    double ajuste, diferencia;
    int i, j, ciclo;

    /* >>> Apertura de archivo para registrar el error por iteración <<< */
    FILE *registroError = fopen("registro_error_seidel.txt", "w");
    if (registroError == NULL) {
        printf("No se pudo abrir el archivo de registro de error.\n");
        return;
    }

    /* Inicialización de nodos en 0 */
    for(i = 0; i < DIM; i++){
        solucion[i] = 0.0;
    }

    printf("\n========================================================================\n");
    printf(" Iter |          n1          n2          n3          n4          n5        |   Error\n");
    printf("========================================================================\n");

    for(ciclo = 1; ciclo <= CICLOS_MAX; ciclo++){

        /* Copia de valores  */
        for(i = 0; i < DIM; i++){
            previo[i] = solucion[i];
        }

        /* Iteración principal */
        for(i = 0; i < DIM; i++){
            ajuste = vector_b[i];

            for(j = 0; j < DIM; j++){
                if(j != i){
                    ajuste -= matriz[i][j] * solucion[j];
                }
            }
            solucion[i] = ajuste / matriz[i][i];
        }

        /* Cálculo de error máximo  */
        diferencia = fabs(solucion[0] - previo[0]);
        for(i = 1; i < DIM; i++){
            double diff = fabs(solucion[i] - previo[i]);
            if(diff > diferencia){
                diferencia = diff;
            }
        }

        /* >>> Guardar iteración y error en el archivo <<< */
        fprintf(registroError, "%d %e\n", ciclo, diferencia);

        /*  fila de iteración en pantalla */
        printf(" %4d | %12lf %12lf %12lf %12lf %12lf | %e\n",
               ciclo, solucion[0], solucion[1], solucion[2], solucion[3], solucion[4], diferencia);

        /* Criterio de parada */
        if(diferencia < LIMITE_ERR){
            printf("========================================================================\n");
            printf(">>> Gauss-Seidel convergió en %d iteraciones con error final = %e\n", ciclo, diferencia);
            
            /* >>> Cerrar archivo antes de salir <<< */
            fclose(registroError);
            return;
        }
    }

    printf("========================================================================\n");
    printf(">>> Gauss-Seidel NO convergió dentro del límite de %d iteraciones.\n", CICLOS_MAX);

    /* >>> Cerrar archivo si no convergió <<< */
    fclose(registroError);
}



/*   Método LU   */

void factorizacionLU(double matriz[DIM][DIM], double L[DIM][DIM], double U[DIM][DIM]){
    int i, j, k;

    for(i = 0; i < DIM; i++){
        for(j = 0; j < DIM; j++){
            L[i][j] = 0.0;
            U[i][j] = 0.0;
        }
    }

    for(k = 0; k < DIM; k++){
        
        for(j = k; j < DIM; j++){
            double suma = 0.0;
            for(i = 0; i < k; i++){
                suma += L[k][i] * U[i][j];
            }
            U[k][j] = matriz[k][j] - suma;
        }

  
        for(i = k; i < DIM; i++){
            if(i == k){
                L[i][k] = 1.0;
            } else {
                double suma = 0.0;
                for(j = 0; j < k; j++){
                    suma += L[i][j] * U[j][k];
                }
                L[i][k] = (matriz[i][k] - suma) / U[k][k];
            }
        }
    }
}

void resolverLU(double L[DIM][DIM], double U[DIM][DIM], double vector_b[DIM], double solucion[DIM]){
    double y[DIM];
    int i, j;

    /* Sustitución hacia adelante: L * y = b */
    for(i = 0; i < DIM; i++){
        double suma = 0.0;
        for(j = 0; j < i; j++){
            suma += L[i][j] * y[j];
        }
        y[i] = vector_b[i] - suma;
    }

    /* Sustitución hacia atrás: U * x = y */
    for(i = DIM-1; i >= 0; i--){
        double suma = 0.0;
        for(j = i+1; j < DIM; j++){
            suma += U[i][j] * solucion[j];
        }
        solucion[i] = (y[i] - suma) / U[i][i];
    }
}

int main(){
    /* Matriz de coeficientes y vector de términos independientes */
    double matrizCoef[DIM][DIM] = {
        { 11,  -3,  -2,   0,   0},
        {-10,  23,  -8,  -5,   0},
        {-35, -42, 142, -30, -35},
        {  0,  -7,  -8,  29, -14},
        {  0,   0, -10, -15,  31}
    };

    double vectorB[DIM] = {60, 0, 1050, 0, 0};

    double salidaJordan[DIM], salidaSeidel[DIM], salidaLU[DIM];
    double replica[DIM][DIM];
    double L[DIM][DIM], U[DIM][DIM];

    /* Gauss–Jordan */
    copiarMatriz(matrizCoef, replica);
    metodo_GaussJordan(replica, vectorB, salidaJordan);
    reporteSolucion(salidaJordan, "Resultados mediante Gauss-Jordan");

    /* Gauss–Seidel */
    copiarMatriz(matrizCoef, replica);
    metodo_GaussSeidel(replica, vectorB, salidaSeidel);
    reporteSolucion(salidaSeidel, "Resultados mediante Gauss-Seidel");

    /* Factorización LU */
    copiarMatriz(matrizCoef, replica);
    factorizacionLU(replica, L, U);
    resolverLU(L, U, vectorB, salidaLU);
    reporteSolucion(salidaLU, "Resultados mediante Factorizacion LU");

    return 0;
}
