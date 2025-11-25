#include <stdio.h>
#include <math.h>

#define N_MAX   200   
#define NP      4     
#define MAX_IT  5000 
#define STEP_INIT  0.1
#define STEP_MIN   1e-6

/* -------- DATOS EXPERIMENTALES  -------- */

int N;                    
double aw[N_MAX];          
double Xe[N_MAX];          

/* -------- MODELOS -------- */

// Modelo Peleg: Xe = b0*aw^b1 + b2*aw^b3
double peleg_model(double aw, const double *b) {
    return b[0]*pow(aw, b[1]) + b[2]*pow(aw, b[3]);
}

// Modelo DLP: Xe = b0 + b1*x + b2*x^2 + b3*x^3, x = ln(-ln aw)
double dlp_model(double aw, const double *b) {
    double x = log(-log(aw));
    return b[0] + b[1]*x + b[2]*x*x + b[3]*x*x*x;
}

/* -------- Cálculo de chi^2 para un modelo genérico -------- */

double chi2(double (*model)(double, const double *), const double *b) {
    double sum = 0.0;
    for (int i = 0; i < N; ++i) {
        double Xe_pred = model(aw[i], b);
        double r = Xe[i] - Xe_pred;
        sum += r * r;
    }
    return sum;
}

/* -------- Minimización sencilla de chi^2 (búsqueda por coordenadas) -------- */

void fit(double (*model)(double, const double *), double *b) {
    double step[NP];
    for (int j = 0; j < NP; ++j) step[j] = STEP_INIT;

    double chi_best = chi2(model, b);

    for (int iter = 0; iter < MAX_IT; ++iter) {
        int improved = 0;

        for (int j = 0; j < NP; ++j) {
            double old = b[j];

            
            b[j] = old + step[j];
            double chi_plus = chi2(model, b);

            if (chi_plus < chi_best) {
                chi_best = chi_plus;
                improved = 1;
                continue;
            }

            /
            b[j] = old - step[j];
            double chi_minus = chi2(model, b);

            if (chi_minus < chi_best) {
                chi_best = chi_minus;
                improved = 1;
                continue;
            }

            
            b[j] = old;
        }

        if (!improved) {
            
            double max_step = 0.0;
            for (int j = 0; j < NP; ++j) {
                step[j] *= 0.5;
                if (step[j] > max_step) max_step = step[j];
            }
            if (max_step < STEP_MIN) break;   
        }
    }
}



int main(void) {
    

    N = 5;  

    
    aw[0] = 0.10; Xe[0] = 2.10;
    aw[1] = 0.20; Xe[1] = 2.20;
    aw[2] = 0.30; Xe[2] = 2.30;
    aw[3] = 0.40; Xe[3] = 2.50;
    aw[4] = 0.50; Xe[4] = 3.00;

   
    double b_peleg[NP] = {1.0, 1.0, 1.0, 1.0};
    double b_dlp[NP]   = {1.0, 1.0, 1.0, 1.0};

   
    fit(peleg_model, b_peleg);
    double chiP = chi2(peleg_model, b_peleg);
    double chiP_red = chiP / (N - NP);

   
    fit(dlp_model, b_dlp);
    double chiD = chi2(dlp_model, b_dlp);
    double chiD_red = chiD / (N - NP);

   
    printf("=== Resultados para esta curva ===\n\n");

    printf("Modelo Peleg:\n");
    printf("  b0 = %g\n  b1 = %g\n  b2 = %g\n  b3 = %g\n",
           b_peleg[0], b_peleg[1], b_peleg[2], b_peleg[3]);
    printf("  chi2      = %g\n", chiP);
    printf("  chi2_red  = %g\n\n", chiP_red);

    printf("Modelo DLP:\n");
    printf("  b0 = %g\n  b1 = %g\n  b2 = %g\n  b3 = %g\n",
           b_dlp[0], b_dlp[1], b_dlp[2], b_dlp[3]);
    printf("  chi2      = %g\n", chiD);
    printf("  chi2_red  = %g\n\n", chiD_red);

    if (chiP_red < chiD_red)
        printf("=> El mejor modelo para esta curva es: Peleg\n");
    else
        printf("=> El mejor modelo para esta curva es: DLP\n");

    return 0;
}
