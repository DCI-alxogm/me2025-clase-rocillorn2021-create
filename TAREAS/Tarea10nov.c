// Tarea10nov.c — Gradiente ascendente (máximo) con derivadas finitas

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


double objective(int n, const double *x);
void   finite_diff_grad(int n, const double *x, double h, double *g);
double norm2(int n, const double *v);
void   gradient_ascent(int n, double *x, double alpha, double h, double tol, int max_iter);
int    read_input(FILE *f, int *n, double **x0, double *alpha, double *h, double *tol, int *max_iter);


double objective(int n, const double *x) {
    if (n == 1) {
        double t = x[0];
        return -(t - 1.5)*(t - 1.5);
    } else if (n == 2) {
        double x1 = x[0], x2 = x[1];
        return - (x1 - 2.0)*(x1 - 2.0) - (x2 + 1.0)*(x2 + 1.0);
    } else {
        double s=0;
        for(int i=0;i<n;++i){
            double c = (double)(i+1);
            double d = x[i] - c;
            s -= d*d;
        }
        return s;
    }
}


void finite_diff_grad(int n, const double *x, double h, double *g){
    double *xp = malloc(n*sizeof(double));
    double *xm = malloc(n*sizeof(double));
    for(int i=0;i<n;++i){ xp[i]=x[i]; xm[i]=x[i]; }

    for(int i=0;i<n;++i){
        xp[i] = x[i] + h;
        xm[i] = x[i] - h;
        double fp = objective(n, xp);
        double fm = objective(n, xm);
        g[i] = (fp - fm)/(2*h);
        xp[i] = x[i];
        xm[i] = x[i];
    }
    free(xp);
    free(xm);
}


double norm2(int n, const double *v){
    double s=0;
    for(int i=0;i<n;++i) s += v[i]*v[i];
    return sqrt(s);
}


void gradient_ascent(int n, double *x, double alpha, double h, double tol, int max_iter){
    double *g  = malloc(n*sizeof(double));
    double *xn = malloc(n*sizeof(double));

    double fx = objective(n, x);

    for(int it=1; it<=max_iter; ++it){
        finite_diff_grad(n, x, h, g);
        double ng = norm2(n, g);
        if(ng < tol){
            printf("Convergencia: ||grad||=%.3e iter=%d\n",ng,it);
            break;
        }

        for(int i=0;i<n;++i) xn[i] = x[i] + alpha*g[i];
        double fn = objective(n, xn);

        double a = alpha;
        int bt=0;
        while(fn < fx){
            a *= 0.5;
            if(a < 1e-16) break;
            for(int i=0;i<n;++i) xn[i] = x[i] + a*g[i];
            fn = objective(n, xn);
            if(++bt > 20) break;
        }

        if(fn > fx){
            for(int i=0;i<n;++i) x[i]=xn[i];
            fx = fn;
        } else {
            if(a < 1e-16){
                printf("Paso mínimo alcanzado\n");
                break;
            }
        }
    }
    printf("f(x*) = %.10f\n",fx);
    free(g);
    free(xn);
}

/*Leer entrada */
int read_input(FILE *f, int *n, double **x0, double *alpha, double *h, double *tol, int *max_iter){
    if(fscanf(f,"%d",n)!=1 || *n<=0) return -1;
    *x0 = malloc((*n)*sizeof(double));
    for(int i=0;i<*n;++i)
        if(fscanf(f,"%lf",&(*x0)[i])!=1) return -1;

    double A,H,T;
    int M;
    if(fscanf(f,"%lf %lf %lf %d",&A,&H,&T,&M)==4){
        *alpha=A; *h=H; *tol=T; *max_iter=M;
    } else {
        *alpha=0.1; *h=1e-6; *tol=1e-6; *max_iter=10000;
    }
    return 0;
}


int main(int argc, char **argv){
    if(argc<2){
        printf("Uso:\n  %s <archivo>\n  %s --stdin\n",argv[0],argv[0]);
        return 1;
    }

    FILE *f;
    int use_stdin = (strcmp(argv[1],"--stdin")==0);
    if(use_stdin) f=stdin;
    else {
        f=fopen(argv[1],"r");
        if(!f){ perror("No pude abrir archivo"); return 1; }
    }

    int n, max_iter;
    double *x0, alpha, h, tol;

    if(read_input(f,&n,&x0,&alpha,&h,&tol,&max_iter)){
        printf("Error leyendo entrada.\n");
        return 1;
    }
    if(!use_stdin) fclose(f);

    printf("n=%d alpha=%.2g h=%.1e tol=%.1e max_iter=%d\n",n,alpha,h,tol,max_iter);
    gradient_ascent(n,x0,alpha,h,tol,max_iter);

    printf("x* =");
    for(int i=0;i<n;++i) printf(" %.6f",x0[i]);
    printf("\n");

    free(x0);
    return 0;
}




