#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define EPS 1e-12


static double **mat(int n){
    double **M = (double**)malloc(n*sizeof(*M));
    if(!M) return NULL;
    for(int i=0;i<n;++i){
        M[i]=(double*)calloc(n,sizeof(**M));
        if(!M[i]){ for(int k=0;k<i;++k) free(M[k]); free(M); return NULL; }
    }
    return M;
}
static void freemat(double **M,int n){ if(!M) return; for(int i=0;i<n;++i) free(M[i]); free(M); }
static double *vec(int n){ return (double*)calloc(n,sizeof(double)); }


static int read_system(FILE *f, int *n, double ***A, double **b){
    if(fscanf(f,"%d",n)!=1 || *n<=0) return -1;
    int N=*n;
    double **AA=mat(N); if(!AA) return -1;
    double *bb=vec(N);  if(!bb){ freemat(AA,N); return -1; }

    for(int i=0;i<N;++i){
        for(int j=0;j<N;++j){
            if(fscanf(f,"%lf",&AA[i][j])!=1){ freemat(AA,N); free(bb); return -1; }
        }
    }
    for(int i=0;i<N;++i){
        if(fscanf(f,"%lf",&bb[i])!=1){ freemat(AA,N); free(bb); return -1; }
    }
    *A=AA; *b=bb; return 0;
}


static int lu_factor(double **A,int n,double ***Lout,double ***Uout,int *piv,double *det_out){
    double **L=mat(n), **U=mat(n);
    if(!L || !U){ freemat(L,n); freemat(U,n); return -1; }

    for(int i=0;i<n;++i){ for(int j=0;j<n;++j) U[i][j]=A[i][j]; L[i][i]=1.0; }

    int sign = 1;
    for(int k=0;k<n;++k){
        int p=k; double maxabs=fabs(U[k][k]);
        for(int i=k+1;i<n;++i){ double v=fabs(U[i][k]); if(v>maxabs){maxabs=v;p=i;} }
        piv[k]=p;
        if(p!=k){
            for(int j=0;j<n;++j){ double t=U[k][j]; U[k][j]=U[p][j]; U[p][j]=t; }
            for(int j=0;j<k;++j){ double t=L[k][j]; L[k][j]=L[p][j]; L[p][j]=t; }
            sign*=-1;
        }
        if(fabs(U[k][k])<EPS){ freemat(L,n); freemat(U,n); return -1; }

        double pivot = U[k][k];
        for(int i=k+1;i<n;++i){
            L[i][k] = U[i][k]/pivot;
            U[i][k] = 0.0;
            for(int j=k+1;j<n;++j) U[i][j] -= L[i][k]*U[k][j];
        }
    }
    double det=sign;
    for(int i=0;i<n;++i) det*=U[i][i];
    if(det_out) *det_out=det;
    *Lout=L; *Uout=U; return 0;
}


static void apply_pivots(int n,const int *piv,double *b){
    for(int k=0;k<n;++k){ int p=piv[k]; if(p!=k){ double t=b[k]; b[k]=b[p]; b[p]=t; } }
}
static void forward(int n,double **L,double *y){         /* L con diag=1 */
    for(int i=0;i<n;++i) for(int j=0;j<i;++j) y[i]-=L[i][j]*y[j];
}
static void back(int n,double **U,double *y,double *x){
    for(int i=n-1;i>=0;--i){
        double s=y[i];
        for(int j=i+1;j<n;++j) s-=U[i][j]*x[j];
        x[i]=s/U[i][i];
    }
}

/* ===== guardar / cargar LU (con verificación de fscanf) ===== */
static int save_lu(const char *path,int n,const int *piv,double **L,double **U){
    FILE *f=fopen(path,"w"); if(!f) return -1;
    fprintf(f,"%d\n",n);
    for(int i=0;i<n;++i) fprintf(f,"%d%c",piv[i], i==n-1?'\n':' ');
    for(int i=0;i<n;++i){
        for(int j=0;j<n;++j) fprintf(f,"%.17g%c",L[i][j], j==n-1?'\n':' ');
    }
    for(int i=0;i<n;++i){
        for(int j=0;j<n;++j) fprintf(f,"%.17g%c",U[i][j], j==n-1?'\n':' ');
    }
    fclose(f); return 0;
}

static int load_lu(const char *path,int *n,int **piv,double ***L,double ***U){
    FILE *f=fopen(path,"r"); if(!f) return -1;
    if(fscanf(f,"%d",n)!=1){ fclose(f); return -1; }
    int N=*n;

    int *pv=(int*)malloc(N*sizeof(int));
    if(!pv){ fclose(f); return -1; }

    for(int i=0;i<N;++i){
        if(fscanf(f,"%d",&pv[i])!=1){ free(pv); fclose(f); return -1; }
    }

    double **Ll=mat(N), **Uu=mat(N);
    if(!Ll || !Uu){ free(pv); freemat(Ll,N); freemat(Uu,N); fclose(f); return -1; }

    for(int i=0;i<N;++i){
        for(int j=0;j<N;++j){
            if(fscanf(f,"%lf",&Ll[i][j])!=1){ free(pv); freemat(Ll,N); freemat(Uu,N); fclose(f); return -1; }
        }
    }
    for(int i=0;i<N;++i){
        for(int j=0;j<N;++j){
            if(fscanf(f,"%lf",&Uu[i][j])!=1){ free(pv); freemat(Ll,N); freemat(Uu,N); fclose(f); return -1; }
        }
    }

    fclose(f);
    *piv=pv; *L=Ll; *U=Uu; return 0;
}

/* ===== leer solo b (con verificación) ===== */
static int load_b(const char *path,int *n,double **b){
    FILE *f=fopen(path,"r"); if(!f) return -1;
    if(fscanf(f,"%d",n)!=1){ fclose(f); return -1; }
    int N=*n;
    double *bb=vec(N); if(!bb){ fclose(f); return -1; }
    for(int i=0;i<N;++i){
        if(fscanf(f,"%lf",&bb[i])!=1){ free(bb); fclose(f); return -1; }
    }
    fclose(f);
    *b=bb; return 0;
}


static void uso(const char *p){
    fprintf(stderr,
      "Uso:\n"
      "  %s <archivo_in>                 # lee n,A,b; resuelve; imprime det y x\n"
      "  %s --stdin                      # lee n,A,b desde terminal\n"
      "  %s <archivo_in> --save lu.txt   # además guarda L,U y pivotes\n"
      "  %s --use-lu lu.txt --b b.txt    # reutiliza LU para nuevo b\n", p,p,p,p);
}

int main(int argc,char **argv){
    if(argc<2){ uso(argv[0]); return 1; }

    
    if(strcmp(argv[1],"--use-lu")==0){
        if(argc!=5 || strcmp(argv[3],"--b")){ uso(argv[0]); return 1; }
        int n, *piv=NULL; double **L=NULL, **U=NULL;
        if(load_lu(argv[2],&n,&piv,&L,&U)!=0){
            fprintf(stderr,"Error al cargar LU.\n"); return 1;
        }
        int nb; double *b=NULL;
        if(load_b(argv[4],&nb,&b)!=0 || nb!=n){
            fprintf(stderr,"Error al leer b o dimensiones incompatibles.\n");
            free(piv); freemat(L,n); freemat(U,n); if(b) free(b);
            return 1;
        }
        apply_pivots(n,piv,b);
        double *x=vec(n); if(!x){ free(piv); freemat(L,n); freemat(U,n); free(b); return 1; }
        forward(n,L,b); back(n,U,b,x);
        for(int i=0;i<n;++i) printf("x[%d] = %.10f\n",i+1,x[i]);
        free(x); free(b); free(piv); freemat(L,n); freemat(U,n);
        return 0;
    }

    FILE *f = NULL; int from_stdin = (strcmp(argv[1],"--stdin")==0);
    if(from_stdin){ f=stdin; }
    else {
        f=fopen(argv[1],"r");
        if(!f){ perror("No pude abrir archivo"); return 1; }
    }

    int n; double **A=NULL,*b=NULL;
    if(read_system(f,&n,&A,&b)!=0){
        fprintf(stderr,"Error leyendo datos (n, A, b).\n");
        if(!from_stdin) fclose(f);
        return 1;
    }
    if(!from_stdin) fclose(f);

    double **L=NULL, **U=NULL; int *piv=(int*)malloc(n*sizeof(int));
    if(!piv){ freemat(A,n); free(b); return 1; }

    double det;
    if(lu_factor(A,n,&L,&U,piv,&det)!=0){
        printf("Sistema no bien determinado (det≈0).\n");
        freemat(A,n); free(b); free(piv); freemat(L,n); freemat(U,n);
        return 1;
    }

    printf("det(A) = %.10e\n", det);
    if(fabs(det) < 1e-14) printf("[AVISO] det≈0 -> problema mal determinado.\n");

    apply_pivots(n,piv,b);
    double *x=vec(n); if(!x){ freemat(A,n); free(b); free(piv); freemat(L,n); freemat(U,n); return 1; }
    forward(n,L,b);
    back(n,U,b,x);

    printf("Solucion:\n");
    for(int i=0;i<n;++i) printf("x[%d] = %.10f\n", i+1, x[i]);

    if(argc==4 && strcmp(argv[2],"--save")==0){
        if(save_lu(argv[3],n,piv,L,U)==0)
            printf("LU guardado en %s\n", argv[3]);
        else
            fprintf(stderr,"No se pudo guardar LU en %s\n", argv[3]);
    }

    freemat(A,n); free(b); free(piv); freemat(L,n); freemat(U,n); free(x);
    return 0;
}

